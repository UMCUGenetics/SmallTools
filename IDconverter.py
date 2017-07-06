import requests, sys
import json
import os.path
import csv
import argparse

import seaborn as sns
import struct

def hex2rgb(rgb):
    return struct.unpack('BBB', rgb.decode('hex'))

def rgb2hex(rgb):
    return struct.pack('BBB',*rgb).encode('hex')

# ------------------------------------------------------------------------
# TODO put into a sensible JSON file
#               "human"                 "chicken"           "zebrafish"         "mouse"
speciesmap = {"homo_sapiens":9606, "gallus_gallus":9031, "danio_rerio":7955, "mus_musculus":10090}
keggmap  = {"homo_sapiens":"hsa", "gallus_gallus":"gga", "danio_rerio":"dre", "mus_musculus":"mmu"}
# ------------------------------------------------------------------------
# TODO put into a JSON file
# make kegg color map
colors = {"nodata":"#ffffff","pos":"#1f78b4","neg":"#a6cee3","test":"#b2df8a", "de":"#33a02c", "E105":"#fdbf6f","E115":"#ff7f00","28hpf_AP":"#fb9a99","28hpf_DV":"#e31a1c","40hpf_AP":"#cab2d6","40hpf_DV":"#6a3d9a"}
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Convert and MAP IDs across species to one KEGG map')
parser.add_argument('-i', '--input',   type=str, help='input json data file')
parser.add_argument('-o', '--outdir',  type=str, help='output directory')
parser.add_argument('-p', '--pathway', type=str, help='pathway to map to')

args = parser.parse_args()
print("-"*60)
print("Running with these settings:")
print(args)
print("-"*60)
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# RESTfull functions
def generic_json_request_handler(server, ext):
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    return(r.json())

def get_kegg_genes(pathway):
    #if args.debug:
    #    print("Parsing pathway: "+pathway)
    server="http://togows.org"
    ext="/entry/kegg-pathway/"+pathway+"/genes.json"
    return(generic_json_request_handler(server, ext)[0])

def get_ens_orthologues(ensid):
    server = "https://rest.ensembl.org"
    ext = "/homology/id/"+ensid+"?format=condensed;type=orthologues"
    return(generic_json_request_handler(server, ext)['data'][0])

def map_ens_to_species(ensid, targettaxon):
    server = "https://rest.ensembl.org"
    ext = "/homology/id/"+ensid+"?format=condensed;type=orthologues;target_taxon="+targettaxon
    return(generic_json_request_handler(server, ext)['data'][0])

def get_sym_orthologues(symbol, species):
    server="https://rest.ensembl.org"
    ext="/homology/symbol/"+species+"/"+symbol+"?format=condensed;type=orthologues"
    return(generic_json_request_handler(server, ext)['data'][0])

def find_symbols(symbols, species):
    server = "https://rest.ensembl.org"
    ext = "/lookup/symbol/"+species
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}

    # requires slightly more formatting of the list
    r = requests.post(server+ext, headers=headers, data='{ "symbols" : ["'+'", "'.join(symbols)+'" ] }')
    #print('{ "symbols" : ["'+'", "'.join(symbols)+'" ] }')
    if not r.ok:
      r.raise_for_status()
      sys.exit()

    decoded = r.json()
    return(decoded)
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
def parse_kegg_genes(species, jsondata):
    decoded = {}
    for kid in jsondata:
        hgcn = jsondata[kid].split(";")[0]
        decoded[hgcn] = {}
        decoded[hgcn]["KEGG"]=keggmap[species]+":"+kid

    mappings = find_symbols(decoded.keys(), species)
    for gene in mappings:
        decoded[gene]["mapping"] = mappings[gene]
        decoded[gene]["orthologues"] = get_ens_orthologues(mappings[gene]['id'])

    for hgcn in decoded:
        # NOT found through lookup let's try something else
        if hgcn not in mappings:
            decoded[hgcn]["orthologues"] = get_sym_orthologues(hgcn, species)
            decoded[hgcn]['id'] = decoded[hgcn]["orthologues"]['id']

    # REPORT incomplete genes
    for i in decoded:
        if not 'orthologues' in decoded[i]:
            print('[WARN]'+'no orthologues found for: %s'%(i))

    return(decoded)
# ------------------------------------------------------------------------
def fill_kegg_colors(genedata, fulldata, colors):
    # INIT colors
    keggcolors = {fulldata[i]["KEGG"]:{} for i in fulldata}


    for species in genedata:
        for condition in genedata[species]:
            # print(species+"_"+condition)
            # STORE DEFAULT VALUE
            for keggid in keggcolors:
                keggcolors[keggid][species+"_"+condition] = colors["nodata"]

            # Retrieve ID mapping data
            ofile = args.outdir+condition+"_mappings.json"
            mappingdata={}
            if not os.path.isfile(ofile):
                #print(genedata[species][condition])
                mappingdata = find_symbols(genedata[species][condition], species)
                with open(ofile, 'w') as outfile:
                    json.dump(mappingdata, outfile)
            else:
                with open(ofile, 'r') as data_file:
                    mappingdata = json.load(data_file)

            print("[INFO] Finished %s"%(species+"_"+condition))

            # CHECK FOR GENE MATCHES
            for agene in mappingdata:
                # IF FULL SYMBOL MATCH
                if agene in fulldata:
                    keggcolors[fulldata[agene]["KEGG"]][species+"_"+condition] = colors[condition]
                    continue

                for i in fulldata:
                    if not 'orthologues' in fulldata[i]:
                        # SKIP IDs without ortholog information
                        continue

                    if mappingdata[agene]['id'] in [fulldata[i]["orthologues"]["homologies"][j]["id"] for j in range(0,len( fulldata[i]["orthologues"]["homologies"]))]:
                        keggcolors[fulldata[i]["KEGG"]][species+"_"+condition] = colors[condition]

    return(keggcolors)
# ------------------------------------------------------------------------
def write_kegg_colors(keggcolors, outputfile):

    pallette = sns.color_palette("Reds",9).as_hex()
    #['#fee5d8', '#fdcab5', '#fcab8f', '#fc8a6a', '#fb694a', '#f14432', '#d92523', '#bc141a', '#980c13']

    # ADD the all column
    for keggid in keggcolors:
        keggcolors[keggid]["all"] = sum([j!=colors["nodata"] for i,j in keggcolors[keggid].items()])

    for keggid in keggcolors:
        if keggcolors[keggid]["all"] == 0:
            keggcolors[keggid]["all"] = colors["nodata"]
        else:
            keggcolors[keggid]["all"] = pallette[keggcolors[keggid]["all"]-1]

    conditions = [i for i in keggcolors[list(keggcolors)[0]]]
    with open(outputfile, 'w') as outfile:
        outfile.write("#hsa\t"+'\t'.join(conditions)+'\n')
        for keggid in keggcolors:
            outfile.write(keggid.split(":")[1]+'\t'+'\t'.join([str(j) for i,j in keggcolors[keggid].items()])+'\n')
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------

fulldata = {}
pathwayfile = args.pathway+"_mappings.json"

if not os.path.isfile(os.path.join(os.sep, args.outdir, args.pathway, pathwayfile)):
    # parse all info for genes in pathway
    fulldata = parse_kegg_genes("homo_sapiens", get_kegg_genes(args.pathway))
    # store in file
    with open(os.path.join(os.sep, args.outdir, args.pathway, pathwayfile), 'w') as outfile:
        json.dump(fulldata, outfile)

else:
    # load data from pre-generated file
    with open(os.path.join(os.sep, args.outdir, args.pathway, pathwayfile), 'r') as data_file:
        fulldata = json.load(data_file)

# ------------------------------------------------------------------------

with open(args.input, 'r') as data_file:
    genedata = json.load(data_file)
    keggcolors = fill_kegg_colors(genedata, fulldata, colors)
    write_kegg_colors(keggcolors, os.path.join(os.sep, args.outdir, args.pathway, args.pathway+"_all_conditions_colors.txt"))

# ------------------------------------------------------------------------
