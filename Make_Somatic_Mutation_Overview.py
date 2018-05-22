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

parser.add_option("--dp",       dest="mindepth",   help="Minimum read depth to consider reliable", default=10)
parser.add_option("--af",       dest="minvaf",     help="Minimum variant allele fraction",         default=0.25)
parser.add_option("--pf",       dest="popfreq",    help="Maximum popultaion frequency",            default=0.05)
parser.add_option("--cf",       dest="cohfreq",    help="Maximum cohort frequency",                default=0.10)
parser.add_option("--me",       dest="mineff",     help="Minimum variant effect score",            default=1.50)

parser.add_option("--debug",    dest="debug",      help="Flag for debug logging",                  default=False)
parser.add_option("--format",   dest="format",     help="VCF output format [GATK/FREEB/..]",       default="GATK")
(options, args) = parser.parse_args()
# -------------------------------------------------

# -------------------------------------------------
vocabulary = {
    "None":-1, "clean":0,
    "sequence_feature":0, "intron_variant":0,
    "3_prime_UTR_variant":0, "5_prime_UTR_variant":0, "non_coding_exon_variant":0,
    "TF_binding_site_variant":0.5, "splice_region_variant":0.5,
    "synonymous_variant":1.0,
    "missense_variant":1.5,
    "splice_donor_variant":2, "splice_acceptor_variant":2,
    "inframe_deletion":2.1, "inframe_insertion":2.1,
    "disruptive_inframe_deletion":2.5, "disruptive_inframe_insertion":2.5,
    "5_prime_UTR_premature_start_codon_gain_variant":3,
    "stop_gained":4, "nonsense_mediated_decay":4, "frameshift_variant":4
}

# Mapping of SNEPeff effects to 'MAF' names for variation effects, enables later use in MAF tools
# https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v1.0
# https://bioconductor.org/packages/3.7/bioc/vignettes/maftools/inst/doc/maftools.html
mapping = {
    "synonymous_variant":"Silent", "missense_variant":"Missense_Mutation", "disruptive_inframe_deletion":"Frame_Shift_Del", "disruptive_inframe_insertion":"Frame_Shift_Ins",
    "5_prime_UTR_premature_start_codon_gain_variant":"Nonsense_Mutation", "stop_gained":"Nonsense_Mutation", "nonsense_mediated_decay":"Nonsense_Mutation", "frameshift_variant":"Frame_Shift_???"
}
# Data fields needed to make lollipop plots
lollipop = ["Hugo_Symbol","Sample_ID","Protein_Change","Mutation_Type","Chromosome","Start_Position","End_Position","Reference_Allele","Variant_Allele","VAF"]

# Known fields with information on population frequency
FREQ_FIELDS = ["dbNSFP_ExAC_AF", "dbNSFP_ExAC_Adj_AF", "GoNLv5_Freq", "GoNLv5_AF"]



# -------------------------------------------------
# DETERMINE which effects to report based on 'abribitrary' variant impact score
toselect = [k for k,v in vocabulary.items() if v >= float(options.mineff)]
# -------------------------------------------------


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
            os.mkdir(options.outdir)
        except OSError:
            print("Invalid / unable to create, output folder %s"%(options.outdir))
            return(False)

    if options.format == "GATK":
        DEPTH_KEY="AD"
        VAF_KEY="AD"

    if options.format == "FREEB":
        DEPTH_KEY="DP"
        VAF_KEY="DPR"


    print("Running with the following settings:")
    print("------------------------------------")
    print(options)
    print("DEPTH FIELD:"+DEPTH_KEY)
    print("ALLELE FIELD:"+VAF_KEY)
    print("------------------------------------")
    return(True)

# -------------------------------------------------

# Extract population frequency from VCF record
# Annoation assumed to be in SNPeff formatting
def find_popfreq(vcf_record):
    popfreq=[0.0]
    for field in FREQ_FIELDS:
        if field in vcf_record.INFO:
            #if debug: print(vcf_record.INFO[field])
            for x in vcf_record.INFO[field]:
                if x is None:
                    popfreq.append(0.0)
                else:
                    popfreq.append(float(x))
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
                    print("NEW Mutation effect identified:")
                    print(pred)
                    print(effect)

            else:
                # STORE THE MOST DELETERIOUS EFFECT
                if vocabulary[effect] > vocabulary[maxeffect]:
                    maxeffect = effect
    if debug: print(maxeffect)
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


# GENE FORMAT
# Gene name + location + variants or not
# VARIANT FORMAT
# Variant + DEPTH + POP FREQ + MLEAF + EFFECT

def check_ad(sample_vcf):
    try:
        ad_item = sample_vcf[DEPTH_KEY]
    except AttributeError as e:
        return(False)
    if sample_vcf[DEPTH_KEY] is None:
        return(False)
    return(True)

#sample_vcf == vcf_record.genotype(sample)
def check_depth(sample_vcf):
    #single depth field
    if isinstance(sample_vcf[DEPTH_KEY], int):
        # SKIP LOW DEPTH POSITIONS
        if sample_vcf[DEPTH_KEY] < int(options.mindepth):
            return(False)
    #multi depth field
    else:
        # SKIP LOW DEPTH POSITIONS
        if sum(sample_vcf[DEPTH_KEY]) < int(options.mindepth):
            return(False)
    return(True)

def check_vaf(sample_vcf):
    #single depth field
    if isinstance(sample_vcf[DEPTH_KEY], int):
        # CHECK VAF
        if (sum(sample_vcf[VAF_KEY][1:])*1.0/sample_vcf[DEPTH_KEY]) < float(options.minvaf):
            return(False)
    #multi depth field
    else:
        # CHECK VAF
        if (sum(sample_vcf[VAF_KEY][1:])*1.0/sum(sample_vcf[DEPTH_KEY])) < float(options.minvaf):
            return(False)
    return(True)

def condense_bed(genelist):
    newlist = {}
    for genebody in genelist:
        gene=genebody[0]
        if gene not in newlist:
            newlist[gene] = [int(genebody[1]),int(genebody[2])]
        else:
            newlist[gene][0] = min(newlist[gene][0],int(genebody[1]))
            newlist[gene][1] = min(newlist[gene][1],int(genebody[2]))

    return([[i,j[0],j[1]] for i,j in newlist.iteritems()])

# -------------------------------------------------

def main():
    global DEPTH_KEY
    global VAF_KEY

    file_list = glob.glob(os.path.join(options.vcfdir, "*.vcf"))
    for vcf_file in file_list:
        zip_and_index(vcf_file)

    genelist=open(options.genelist, 'r').read().split('\n')
    genelist=condense_bed(genelist)

    # DF to keep the mutation effcts per gene
    df = {}
    #VCF record df, for MAX effects only, used for lollipop data
    rdf= {}
    #Count data frame
    cdf = {}

    # FOR ALL VCF FILES
    for vcf_file in file_list:
        if (debug):
            print("------")
            print(vcf_file)
        vcfread = vcf.Reader(open(vcf_file+".gz",'r'), compressed="gz")

        if (debug): print(vcfread.samples)
        if (debug): print(options.format)

        # FOR EACH SAMPLE
        for i,sample in enumerate(vcfread.samples):
            samplename = False

            if options.format == "GATK":
                samplename = sample
            elif options.format == "FREEB":
                if (debug): print("++ "+vcfread.samples[1])
                samplename = vcfread.samples[i+1]
                #samplename = vcf_file.split(".")[1].split("_")[1]
            df[samplename] = {}
            rdf[samplename] = {}
            cdf[samplename] = {}

        if debug: print(df)

        # FOR EACH GENE OF INTREST
        for gene in genelist:
            nr_of_positions = 0
            if len(gene)<=0:
                continue
            thisgene = dict(zip(["Chr","Start","Stop","SYMBOL"], gene.strip().split('\t')))

            #if debug: print(")
            vcf_records=False
            try:
                vcf_records = vcfread.fetch(thisgene["Chr"], int(thisgene["Start"])-20, int(thisgene["Stop"])+20)
            except ValueError as e:
                if debug: print("-- {}\tNO RECORDS FOUND".format(thisgene))
                for samplename in df:
                    df[samplename][thisgene["SYMBOL"]] = "None"
                continue

            # Prep containers
            effects = {}
            records = {}
            for samplename in df:
                effects[samplename] = []
                records[samplename] = []

            # For each variant position within gene
            for vcf_record in vcf_records:
                nr_of_positions += 1
                # For each sample
                for samplename in df:
                    #CHECK IF SAMPLE GENOTYPE AVAILABLE
                    sgenot = None
                    try:
                        sgenot = vcf_record.genotype(samplename)
                    except AttributeError as e:
                        if debug: print("-- {}\t{}\tNO GT FOUND".format(thisgene, samplename))
                        continue

                    # FILTER NON-QC RECORDS
                    PASS = False
                    log = "++ {}\t{}\t{}".format(thisgene,samplename,vcf_record)
                    # CHEK IF AD FIELD PRESENT
                    if check_ad(sgenot):
                        log += "\tAD:PASS"
                        log += "\tDEPTH:{}".format(vcf_record.genotype(samplename)[DEPTH_KEY])
                        # CHECK TOTAL COVERAGE OF IDENTIFIED ALLELLES
                        if check_depth(sgenot):
                            log += ":PASS"
                            log += "\tVAF:{}".format(sum(vcf_record.genotype(samplename)[VAF_KEY][1:])*1.0/sum(vcf_record.genotype(samplename)[DEPTH_KEY]))

                            # add clean if sufficient depth is measured
                            effects[samplename].append("clean")
                            records[samplename].append(None)

                            # CHECK VARIANT ALLELE FREQUENCY
                            if check_vaf(sgenot):
                                log +=":PASS"
                                log +="\tPOP:{}".format([vcf_record.INFO[rf] for rf in FREQ_FIELDS if rf in vcf_record.INFO])
                                # CHECK POPULATION FREQUENCY
                                if max(find_popfreq(vcf_record)) <= float(options.popfreq):
                                    log += ":PASS"
                                    log += "\tMLEAF:{}".format(vcf_record.INFO["MLEAF"])
                                    # CHECK OCCURENCE IN TOTAL POOL
                                    if max(vcf_record.INFO["MLEAF"]) <= float(options.cohfreq):
                                        log +=":PASS"
                                        PASS = True

                    if debug: print(log)
                    if PASS:
                        effects[samplename].append(find_effects(vcf_record))
                        records[samplename].append(vcf_record)

            # ON GENE+SAMPLE LEVEL determine the number of mutations and the maximum mutation effect
            for samplename in df:
                # If no murtations/effects measured consider the gene as 'not assesed'
                if len(effects[samplename]) <= 0:
                    df[samplename][thisgene["SYMBOL"]] = "None"
                    cdf[samplename][thisgene["SYMBOL"]] = 0

                # Else determine the max effect
                else:
                    cdf[samplename][thisgene["SYMBOL"]] = sum([eff in toselect for eff in effects[samplename]])
                    #len(effects[samplename]) - effects[samplename].count("clean")
                    loc = select_maximum_effect(effects[samplename])
                    eff = effects[samplename][loc]

                    # If a 'strong enough' effect is detected report it in the summary
                    if eff in toselect:
                        df[samplename][thisgene["SYMBOL"]] = eff
                        if eff in mapping:
                            rdf[samplename][thisgene["SYMBOL"]] = {}
                            rdf[samplename][thisgene["SYMBOL"]]["REC"] = records[samplename][loc]
                            rdf[samplename][thisgene["SYMBOL"]]["EFF"] = eff

                    # Else check if gene was not observed 'None' or not mutated 'clean'
                    else:
                        # check number of 'clean' positions
                        # if 50% of positions passes DP metric count as clean
                        if effects[samplename].count("clean") >= (nr_of_positions/2):
                            df[samplename][thisgene["SYMBOL"]] = "clean"
                        else:
                            df[samplename][thisgene["SYMBOL"]] = "None"

                if debug: print("** {}\t{}\t{}\t{}\t{}".format(thisgene, samplename, df[samplename][thisgene["SYMBOL"]], cdf[samplename][thisgene["SYMBOL"]], ",".join(effects[samplename])))


    # Printing the mutation overview table
    outfile = open(options.outdir+"/"+"MutationOverview.txt",'w')
    # Print header with gene names
    firstsample = list(df.keys())[0]
    outfile.write("Sample\t{}\n".format('\t'.join(df[firstsample].keys()) ))
    if debug: print("##############################")
    # Loop all samples
    for sp in df:
        if debug: print("{}\t{}\n".format(sp, '\t'.join(df[sp].values()) ))
        outfile.write("{}\t{}\n".format(sp, '\t'.join(df[sp].values()) ))

    if debug: print("##############################")
    outfile.close()


    # Printing the mutation count table
    outfile = open(options.outdir+"/"+"MutationCounts.txt",'w')
    # Print header with gene names
    firstsample = list(cdf.keys())[0]
    outfile.write("Sample\t{}\tTotMutCount\n".format('\t'.join(cdf[firstsample].keys()) ))
    if debug: print("##############################")
    # Loop all samples
    for sp in cdf:
        if debug: print("{}\t{}\t{}\n".format(sp, '\t'.join([str(i) for i in cdf[sp].values()]), sum(cdf[sp].values()) ))
        outfile.write("{}\t{}\t{}\n".format(sp, '\t'.join([str(i) for i in cdf[sp].values()]), sum(cdf[sp].values()) ))

    if debug: print("##############################")
    outfile.close()


    # Printing the mutation details chart/table
    outfile = open(options.outdir+"/"+"MutationChart.txt",'w')
    # Printing annotations header
    outfile.write("{}\n".format('\t'.join(lollipop)))

    if debug: print("##############################")
    for samplename in rdf:
        for gene in rdf[samplename]:
            thisrec = rdf[samplename][gene]["REC"]
            vaf=(sum(thisrec.genotype(samplename)[VAF_KEY][1:])*1.0)/sum(thisrec.genotype(samplename)[DEPTH_KEY])

            proteffect=None
            for pred in thisrec.INFO["ANN"]:
                # Look for the first transcript with this effect
                if rdf[samplename][gene]["EFF"] in pred.split("|")[1].split("&"):
                    proteffect=pred.split("|")[10]
                    break

            if (debug):
                #print(thisrec.INFO["ANN"])
                print(gene, samplename, proteffect, mapping[rdf[samplename][gene]["EFF"]], str(thisrec.CHROM), str(thisrec.POS), str(thisrec.POS+len(thisrec.ALT[0])), thisrec.REF, str(thisrec.ALT[0]), vaf)

            outfile.write("\t".join([gene, samplename, proteffect, mapping[rdf[samplename][gene]["EFF"]], str(thisrec.CHROM), str(thisrec.POS), str(thisrec.POS+len(thisrec.ALT[0])), thisrec.REF, str(thisrec.ALT[0]), str(vaf)])+"\n")
    if debug: print("##############################")
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
