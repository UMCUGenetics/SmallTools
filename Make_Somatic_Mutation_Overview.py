#!/usr/bin/python
import os
import re
import vcf
import glob
import numpy as np

#GENE FORMAT
##chr    start    stop    name
#3    178866311    178952497    PIK3CA

from optparse import OptionParser
# -------------------------------------------------
parser = OptionParser()
parser.add_option("--vcfdir",   dest="vcfdir",     help="Path to directory containing VCF files",  default=False)
parser.add_option("--outdir",   dest="outdir",     help="Path to directory to write output to",    default="./DriverProfile/")
parser.add_option("--genelist", dest="genelist",   help="File containing Genes to test/plot)",     default=False)

parser.add_option("--bgzip",    dest="bgzip",      help="Path to bgzip binary",                    default="bgzip")
parser.add_option("--tabix",    dest="tabix",      help="Path to tabix binary",                    default="tabix")

parser.add_option("--t",        dest="nrcpus",     help="Number of CPUs to use per sample",        default=2)

parser.add_option("--dp",       dest="mindepth",   help="Minumum read depth to consider reliable", default=10)
parser.add_option("--af",       dest="minvaf",     help="Minumum variant allele fraction",         default=0.25)
parser.add_option("--pf",       dest="popfreq",    help="Maximum popultaion frequency",            default=0.05)

parser.add_option("--debug",    dest="debug",      help="Flag for debug logging",                  default=False)
parser.add_option("--format",   dest="format",     help="VCF output format [GATK/FREEB/..]",       default="GATK")
(options, args) = parser.parse_args()
# -------------------------------------------------

vocabulary = {"None":-1, "clean":0, "sequence_feature":0, "synonymous_variant":0, "intron_variant":0, "3_prime_UTR_variant":0.5, "5_prime_UTR_variant":0.5, "non_coding_exon_variant":0.5, "TF_binding_site_variant":1.0, "splice_region_variant":1.1, "missense_variant":1.5, "splice_donor_variant":2, "splice_acceptor_variant":2, "inframe_deletion":2.1, "inframe_insertion":2.1, "disruptive_inframe_deletion":2.5, "disruptive_inframe_insertion":2.5, "5_prime_UTR_premature_start_codon_gain_variant":3, "stop_gained":4, "nonsense_mediated_decay":4, "frameshift_variant":5}
toselect = [k for k,v in vocabulary.items() if v >= 1.5]

mapping = {"missense_variant":"Missense_Mutation", "disruptive_inframe_deletion":"Frame_Shift_Del", "disruptive_inframe_insertion":"Frame_Shift_Ins", "5_prime_UTR_premature_start_codon_gain_variant":"Nonsense_Mutation", "stop_gained":"Nonsense_Mutation", "nonsense_mediated_decay":"Nonsense_Mutation", "frameshift_variant":"Frame_Shift_???"}
lolipop = ["Hugo_Symbol","Sample_ID","Protein_Change","Mutation_Type","Chromosome","Start_Position","End_Position","Reference_Allele","Variant_Allele"]
# -------------------------------------------------
debug = options.debug
DEPTH_KEY=""
VAF_KEY=""
# -------------------------------------------------
def check_arguments():
    global DEPTH_KEY
    global VAF_KEY

    if not os.path.exists(options.vcfdir):
        print("Invalid VCF folder %s"%(options.vcfdir))
        return False

    if not os.path.exists(options.outdir):
        print("Creating output folder %s"%(options.outdir))
        try:
            os.makedir(options.outdir)
        except OSError:
            print("Invalid / unable to create, output folder %s"%(options.outdir))
            return False

    if options.format == "GATK":
        DEPTH_KEY="AD"
        VAF_KEY="AD"

    if options.format == "FREEB":
        DEPTH_KEY = "DP"
        VAF_KEY = "DPR"


    print("Running with the following settings:")
    print("------------------------------------")
    print(options)
    print("DEPTH FIELD:"+DEPTH_KEY)
    print("ALLELE FIELD:"+VAF_KEY)
    print("------------------------------------")
    return True

# -------------------------------------------------

# Extract population frequency from VCF record
# Annoation assumed to be in SNPeff formatting
def find_popfreq(vcf_record):
    popfreq=[0.0]
    #print(vcf_record.INFO)
    freq_fields = ["dbNSFP_ExAC_AF", "dbNSFP_ExAC_Adj_AF", "GoNLv5_Freq"]

    for field in freq_fields:
        if field in vcf_record.INFO:
            #print vcf_record.INFO[field]
            popfreq.append([float(x) for x in vcf_record.INFO[field]])

    #print(popfreq)
    return(popfreq)

# Determine the most damaging effect of the variant
def find_effects(vcf_record):
    maxeffect="None"
    if debug: print(vcf_record.INFO)

    if "ANN" not in vcf_record.INFO:
        return maxeffect

    # TRAVERSE ALL ANNOTATIONS
    for pred in vcf_record.INFO["ANN"]:
        # SPLIT THE SEPERATE FIELDS WITHIN THE ANNOTATION
        items = pred.split("|")
        allele = items[0]
        effects = items[1].split("&")
        for effect in effects:
            if effect not in vocabulary:
                # A NEW MUTATION EFFECT WAS FOUND
                if debug:
                    print "NEW Mutation effect identified:"
                    print pred
                    print effect

            else:
                # STORE THE MOST DELETERIOUS EFFECT
                if vocabulary[effect] > vocabulary[maxeffect]:
                    maxeffect = effect
    if debug: print maxeffect
    return(maxeffect)

# ETRACT THE MOST DELETERIOUS MUTATIONS IN A GENE
def select_maximum_effect(effects):
    effectvalues = [vocabulary[eff] for eff in effects]
    if debug: print(effectvalues)
    indices = np.argmax(effectvalues)
    return(indices)

# CHECK AND GENERATE GZ AND TBI
def zip_and_index(vcffile):
    if not os.path.exists(vcffile+".gz"):
        os.system(options.bgzip+" -c "+vcffile+" > "+vcffile+".gz")
    if not os.path.exists(vcffile+".gz"+".tbi"):
        os.system(options.tabix+" "+vcffile+".gz")

# -------------------------------------------------

def main():
    global DEPTH_KEY
    global VAF_KEY

    file_list = glob.glob(os.path.join(options.vcfdir, "*.vcf"))
    for vcf_file in file_list:
        zip_and_index(vcf_file)

    genelist=open(options.genelist, 'r').read().split('\n')

    df = {}
    rdf= {}

    for vcf_file in file_list:
        if (debug):
            print "------"
            print vcf_file
        vcfread = vcf.Reader(open(vcf_file+".gz",'r'), compressed="gz")

        sample = False
        samplename = False
        if (debug): print options.format
        if options.format == "GATK":
            sample = vcfread.samples[0]
            samplename = sample
        elif options.format == "FREEB":
            if (debug): print("++ "+vcfread.samples[1])
            sample = vcfread.samples[1]
            samplename = vcf_file.split(".")[1].split("_")[1]
            if (debug): print("-- "+samplename)

        if (debug):
            print(vcfread.samples)
            print(sample)

        if not sample:
            print("Error, no sample found "+vcf_file)
            continue

        df[samplename] = {}
        rdf[samplename] = {}

        # FOR EACH GENE OF INTREST
        for gene in genelist:
            if len(gene)<=0:
                continue
            thisgene = dict(zip(["Chr","Start","Stop","SYMBOL"], gene.strip().split('\t')))

            if (debug):    print(thisgene)

            # FOR EACH TUMOR SAMPLE
            vcf_records=False
            try:
                vcf_records = vcfread.fetch(thisgene["Chr"], int(thisgene["Start"])-20, int(thisgene["Stop"])+20)
            except ValueError as e:
                df[samplename][thisgene["SYMBOL"]] = "None"
                continue

            effects = []
            records = []
            # FILTER NON-QC RECORDS
            for vcf_record in vcf_records:
                #CHEK IF AD FIELD PRESENT
                #if not DEPTH_KEY in vcf_record.genotype(sample):
                #    print("Error, key not found "+DEPTH_KEY)
                #    continue

                # CHECK TOTAL COVERAGE OF IDENTIFIED ALLELLES
                if isinstance(vcf_record.genotype(sample)[DEPTH_KEY], int):
                    # SKIP LOW DEPTH POSITIONS
                    if vcf_record.genotype(sample)[DEPTH_KEY] < int(options.mindepth):
                        continue
                    if debug: print sum(vcf_record.genotype(sample)[VAF_KEY][1:])*1.0/vcf_record.genotype(sample)[DEPTH_KEY]
                    # CHECK VAF
                    if (sum(vcf_record.genotype(sample)[VAF_KEY][1:])*1.0/vcf_record.genotype(sample)[DEPTH_KEY]) < float(options.minvaf):
                        continue

                else:
                    # SKIP LOW DEPTH POSITIONS
                    if sum(vcf_record.genotype(sample)[DEPTH_KEY]) < int(options.mindepth):
                        continue
                    if debug: print sum(vcf_record.genotype(sample)[VAF_KEY][1:])*1.0/sum(vcf_record.genotype(sample)[DEPTH_KEY])
                    # CHECK VAF
                    if (sum(vcf_record.genotype(sample)[VAF_KEY][1:])*1.0/sum(vcf_record.genotype(sample)[DEPTH_KEY])) < float(options.minvaf):
                        continue

                # CHECK POPULATION FREQUENCY
                if max(find_popfreq(vcf_record)) > float(options.popfreq):
                    continue

                effects.append(find_effects(vcf_record))
                records.append(vcf_record)

            if len(effects) <= 0:
                df[samplename][thisgene["SYMBOL"]] = "None"
            else:
                loc = select_maximum_effect(effects)
                eff = effects[loc]
                if eff in toselect:
                    df[samplename][thisgene["SYMBOL"]] = eff
                    if eff in mapping:
                        rdf[samplename][thisgene["SYMBOL"]] = {}
                        rdf[samplename][thisgene["SYMBOL"]]["REC"] = records[loc]
                        rdf[samplename][thisgene["SYMBOL"]]["EFF"] = eff
                else:
                    df[samplename][thisgene["SYMBOL"]] = "None"


        #print(sample, df[sample])
        if (debug): print "Sample\t"+'\t'.join(df[samplename].keys())

    outfile = open(options.outdir+"/"+"MutationOverview.txt",'w')
    outfile.write("Sample\t"+'\t'.join(df[samplename].keys())+"\n")
    #print "##############################"
    for sp in df:
        outfile.write(sp+'\t'+'\t'.join(df[sp].values())+"\n")
    #    print sp+'\t'+'\t'.join(df[sp].values())
    #print "##############################"
    outfile.close()

    outfile = open(options.outdir+"/"+"MutationChart.txt",'w')
    outfile.write("Sample\t"+'\t'.join(lolipop)+"\n")
    #print "##############################"
    for samplename in rdf:
        for gene in rdf[samplename]:
            thisrec = rdf[samplename][gene]["REC"]

            proteffect=None
            for pred in thisrec.INFO["ANN"]:
                if rdf[samplename][gene]["EFF"] in pred.split("|")[1].split("&"):
                    proteffect=pred.split("|")[10]
                    break
            if (debug):
                print(thisrec.INFO["ANN"])
                print(gene, samplename, proteffect, mapping[rdf[samplename][gene]["EFF"]], str(thisrec.CHROM), str(thisrec.POS), str(thisrec.POS+len(thisrec.ALT[0])), thisrec.REF, str(thisrec.ALT[0]))
            outfile.write("\t".join([gene, samplename, proteffect, mapping[rdf[samplename][gene]["EFF"]], str(thisrec.CHROM), str(thisrec.POS), str(thisrec.POS+len(thisrec.ALT[0])), thisrec.REF, str(thisrec.ALT[0])] )+"\n")
    #print "##############################"
    outfile.close()




# -------------------------------------------------

print("Starting Analysis")

if __name__ == '__main__':
    if check_arguments():
        main()
    else:
        print("Error in provided arguments")

print("DONE")

# -------------------------------------------------
