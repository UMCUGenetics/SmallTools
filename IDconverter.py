import requests, sys
import json
import os.path
import csv
# ------------------------------------------------------------------------
#               "human"                 "chicken"           "zebrafish"         "mouse"
speciesmap = {"homo_sapiens":9606, "gallus_gallus":9031, "danio_rerio":7955, "mus_musculus":10090}
# ------------------------------------------------------------------------
wdir="/Users/jdeligt/data/Laurent_FunctionalAnalysis/KEGGmapping/"
# ------------------------------------------------------------------------

def get_kegg_data(pathway):
    server="http://togows.org"
    ext="/entry/kegg-pathway/"+pathway+"/genes.json"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()[0]
    return(decoded)

def find_symbols(species, symbols):
    server = "https://rest.ensembl.org"
    ext = "/lookup/symbol/"+species
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}

    r = requests.post(server+ext, headers=headers, data='{ "symbols" : ["'+'", "'.join(symbols)+'" ] }')
    #print('{ "symbols" : ["'+'", "'.join(symbols)+'" ] }')
    if not r.ok:
      r.raise_for_status()
      sys.exit()

    decoded = r.json()
    return(decoded)

def retrieve_orthologs(ensid):
    server = "https://rest.ensembl.org"
    ext = "/homology/id/"+ensid+"?format=condensed;type=orthologues"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()['data'][0]
    return(decoded)


def parse_kegg_genes(species, jsondata):
    decoded = {}
    for kid in jsondata:
        hgcn = jsondata[kid].split(";")[0]
        decoded[hgcn] = {}
        decoded[hgcn]["KEGG"]=species+":"+kid

    mappings = find_symbols("homo_sapiens", decoded.keys())
    for gene in mappings:
        decoded[gene]["mapping"] = mappings[gene]
        decoded[gene]["orthologs"] = retrieve_orthologs(mappings[gene]['id'])

    del(mappings)
    return(decoded)
# ------------------------------------------------------------------------
def update_and_write_kegg_colors(genedata, fulldata, colors, pathwayid):

    # INIT colors
    keggcolors = {fulldata[i]["KEGG"]:{} for i in fulldata}

    for species in genedata:
        for condition in genedata[species]:
            # print(species+"_"+condition)
            # STORE DEFAULT VALUE
            for keggid in keggcolors:
                keggcolors[keggid][species+"_"+condition] = colors["nodata"]

            # Retrieve ID mapping data
            ofile = wdir+condition+"_mappings.json"
            mappingdata={}
            if not os.path.isfile(ofile):
                mappingdata = find_symbols(species, genedata[species][condition])
                with open(ofile, 'w') as outfile:
                    json.dump(mappingdata, outfile)
            else:
                with open(ofile, 'r') as data_file:
                    mappingdata = json.load(data_file)

            # CHECK FOR GENE MATCHES
            for agene in mappingdata:
                # IF FULL SYMBOL MATCH
                if agene in fulldata:
                    keggcolors[fulldata[agene]["KEGG"]][species+"_"+condition] = colors[condition]
                    continue

                for i in fulldata:
                    #print(fulldata[i]["orthologs"]["homologies"][0])
                    if not 'orthologs' in fulldata[i]:
                        continue
                    #print([fulldata[i]["orthologs"]["homologies"][j]["id"] for j in range(0,len( fulldata[i]["orthologs"]["homologies"]))])
                    if mappingdata[agene]['id'] in [fulldata[i]["orthologs"]["homologies"][j]["id"] for j in range(0,len( fulldata[i]["orthologs"]["homologies"]))]:
                        keggcolors[fulldata[i]["KEGG"]][species+"_"+condition] = colors[condition]

            # WRITE PER CODITION FILE
            # with open(wdir+pathwayid+"/"+species+"_"+condition+"_colors.csv", 'wb') as csv_file:
            #    writer = csv.writer(csv_file)
            #    for key, value in keggcolors.items():
            #       writer.writerow([key, value])

    conditions = list(keggcolors[list(keggcolors)[0]])
    #print(conditions)
    with open(wdir+pathwayid+"/"+"All_Codition_Colors.txt", 'w') as outfile:
        outfile.write("#hsa\t"+'\t'.join(conditions)+'\n')
        for keggid in keggcolors:
            #print(keggid)
            #print([j for i,j in keggcolors[keggid].items()])
            outfile.write(keggid.split(":")[1]+'\t'+'\t'.join([j for i,j in keggcolors[keggid].items()])+'\n')
# ------------------------------------------------------------------------
#pathwayid="hsa04640"
pathwayid="hsa05146"
# ------------------------------------------------------------------------
fulldata = {}
pathwayfile = pathwayid+"_mappings.json"

if not os.path.isfile(wdir+pathwayid+"/"+pathwayfile):
    # parse all info for genes in pathway:
    fulldata = parse_kegg_genes("hsa", get_kegg_data(pathwayid))
    # store in file
    with open(wdir+pathwayid+"/"+pathwayfile, 'w') as outfile:
        json.dump(fulldata, outfile)

else:
    # load data from pre-generated file
    with open(wdir+pathwayid+"/"+pathwayfile, 'r') as data_file:
        fulldata = json.load(data_file)

#print(fulldata)

# ------------------------------------------------------------------------
# make kegg color maps
colors = {"nodata":"#ffffff","pos":"#1f78b4","neg":"#a6cee3","test":"#b2df8a", "de":"#33a02c", "E105":"#fdbf6f","E115":"#ff7f00","28hpf_AP":"#fb9a99","28hpf_DV":"#e31a1c","40hpf_AP":"#cab2d6","40hpf_DV":"#6a3d9a"}

with open(wdir+"degene_data.json", 'r') as data_file:
    genedata = json.load(data_file)
    update_and_write_kegg_colors(genedata, fulldata, colors, pathwayid)












# --- LEGACY ----

def convert_and_map_human_gene_symbol(gene, target_taxon):
    server = "https://rest.ensembl.org"

    if "-" in gene:
        gene = gene.replace("-","%2D")
    ext = "/homology/symbol/human/"+gene+"?target_taxon="+str(target_taxon)+";format=condensed;type=orthologues"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()['data'][0]
    return(decoded)



def convert_symbol_file(fname, target_taxon):
    mappings = {}
    with open(fname) as f:
        header = f.readline()
        header = header + '\n'+ f.readline()
        for line in f:
            gene=line.strip()
            mappings[gene] = convert_and_map_human_gene_symbol(gene, target_taxon)

    return(mappings)

def output_relevant_ids(fname, mappings):
    outf = open(fname, 'w')
    outf.write('\t'.join(["SYMBOL","HUMAN","MOUSE"])+'\n')
    for gid in mappings:
        #print(gid)
        #print(mappings[gid])
        if len(mappings[gid]['homologies']) <= 0:
            outf.write('\t'.join([gid, mappings[gid]['id'], ""])+'\n')
        else:
            outf.write('\t'.join([gid, mappings[gid]['id'], mappings[gid]['homologies'][0]['id']])+'\n')
    outf.close()

#origin_taxon=9606
#target_taxon=10090
#base_folder="./"

# FUSION_DOWN
#mappings = convert_symbol_file(base_folder+"SETLUR_PROSTATE_CANCER_TMPRSS2_ERG_FUSION_DN.txt", target_taxon)
#output_relevant_ids(base_folder+"FUSION_DOWN.txt", mappings)

# FUSION_UP
#mappings = convert_symbol_file(base_folder+"SETLUR_PROSTATE_CANCER_TMPRSS2_ERG_FUSION_UP.txt", target_taxon)
#output_relevant_ids(base_folder+"FUSION_UP.txt", mappings)
