#!/opt/local/bin/python2.7

import glob, os

def shellquote(s):
    return "'" + s.replace("'", "'\\''").replace("\"","") + "'"


samples = {}
#os.chdir("/Users/test/data/CNVs/freec/new_analysis/")
controls = ["BLOOD","BULK"]

filenames = glob.glob("*_CNVs.txt")
filenames.sort()


for cnvfile in filenames:
	sample,condition,rest = cnvfile.split("_")

	if (sample not in samples):
		samples[sample] = {"control":False, "derivatives":[]}
	#print sample, condition

	if (condition in controls):
		samples[sample]["control"] = cnvfile

	else:
		samples[sample]["derivatives"].append(cnvfile)

#print(samples)



# filter CNVs
outfiles = {}
for sample in samples:
	if not samples[sample]["control"]:
		continue
	
	cont = samples[sample]["control"]
	for test in samples[sample]["derivatives"]:
		out = test.replace("_CNVs.txt", "_CNVs.filtered.txt")
		os.system("bedtools intersect -a %s -b overlapmap/ALL_Controls_multisample_merged_CNVs.bed -v | bedtools intersect -a stdin -b %s -v | awk '{if(($3-$2)>50000) print}' > %s" %(test, cont, out))
		outfiles[out] = sample

wind_size = 50000
regions = []

rats = [i.replace("_CNVs.txt","_ratio.txt") for i in filenames]
samplenames = [i.replace("_CNVs.txt","") for i in filenames]


cnvout = "ratios/All_CNVs_merged_ratios.txt"
writer = open(cnvout,'w')
writer.write("Chr\tPos\t"+'\t'.join(samplenames)+"\n")
writer.close()

# gather plot data
for out in outfiles:
	reader = open(out,'r')
	#sample = outfiles[out]

	for region in reader:
		chrom,start,stop,cn,typ = region.strip().split('\t')
		reg = chrom+":"+start+stop
		if reg in regions:
			continue

		regions.append(reg)

		start = max(0,str(int(start)-wind_size))
		stop  = str(int(stop) +wind_size)
		

		contfile = samples[sample]["control"].replace("_CNVs.txt","_ratio.txt")
		paster = "paste <(awk '{if($1==%s && $2>=%s && $2<=%s) print}' %s | cut -f 1-2,4)"%(chrom,start,stop, rats[0])

		for test in rats[1:]:
			scommand = " <(awk '{if($1==%s && $2>=%s && $2<=%s) print}' %s | cut -f 4)"%(chrom,start,stop, test)
			paster += scommand

		paster+=" >> "+cnvout
		
		paster = paster.replace("$","\\$")
		cmd = "/bin/bash -c \""+paster+"\""
		os.system(cmd)

os.system("mv *_CNVs.filtered.txt results/")






# Log for generation of reference map and some usefull use-cases
"""
cat *_BULK_CNVs.txt > ALL_Controls_CNVs.bed
cat *_BLOOD_CNVs.txt >> ALL_Controls_CNVs.bed

cat ALL_Controls_CNVs.bed | sort -k1,1 -k2,2n > ALL_Controls_sorted_CNVs.bed

#/data_fedor12/common_scripts/bedtools2_latest/bedtools2/bin/bedtools merge -i ALL_Controls_sorted_CNVs.bed -d 20000 -c 1 -o count > ALL_Controls_sorted_merged_CNVs.bed
#cat ALL_Controls_sorted_merged_CNVs.bed | awk '{if($4>=2) print}' > ALL_Controls_sorted_merged_multisample_CNVs.bed
#/data_fedor12/common_scripts/bedtools2_latest/bedtools2/bin/bedtools intersect -wa -wb -a ALL_Controls_sorted_merged_multisample_CNVs.bed -b controls/*_CNVs.txt -filenames > overlaps.txt

bedtools intersect -wa -wb -a ALL_Controls_sorted_CNVs.bed -b controls/*_CNVs.txt -filenames > overlaps2.txt

bedtools groupby -i overlaps2.txt -grp 1-3,6 -c 1 -o count > overlaps2.counts
bedtools groupby -i overlaps2.counts -grp 1-3 -c 1 -o count > overlaps2.counts.counts

cat overlaps2.counts.counts | awk '{if ($4>=2) print}' > ALL_Controls_CNVs_multisample.bed

# bedtools merge -i ALL_Controls_CNVs_multisample.bed -c 1 -o count > ALL_Controls_multisample_merged_CNVs.bed
# cat overlapmap/ALL_Controls_CNVs_multisample.bed | sort -k1,1 -k2,2n > UnmatchedControl_CNV_Regions.bed

bedtools intersect -a STE0071_C21_CNVs.txt -b overlapmap/ALL_Controls_multisample_merged_CNVs.bed -v | bedtools intersect -a stdin -b STE0071_BLOOD_CNVs.txt -v | awk '{if(($3-$2)>50000) print}'



cat STE0071_BLOOD_ratio.txt | awk '{if($1==1 && $2>=168300000 && $2 <= 168350000) print}'
cat STE0071_*_ratio.txt | awk '{if($1==1 && $2>=168300000 && $2<=168350000) print}'

cat STE0071_*_ratio.txt | awk '{if($1==1 && $2==168315001) print}'
cat *_ratio.txt | awk '{if($1==1 && $2==168315001) print}'


bedtools intersect -a STE0071_C21_CNVs.txt -b overlapmap/ALL_Controls_multisample_merged_CNVs.bed -v | bedtools intersect -a stdin -b STE0071_BLOOD_CNVs.txt -v | awk '{if(($3-$2)>50000) print}'
bedtools intersect -a STE0097_C_CNVs.txt -b overlapmap/ALL_Controls_multisample_merged_CNVs.bed -v | bedtools intersect -a stdin -b STE0097_BLOOD_CNVs.txt -v | awk '{if(($3-$2)>50000) print}'
bedtools intersect -a C30913D_SC109_CNVs.txt -b overlapmap/ALL_Controls_multisample_merged_CNVs.bed -v | bedtools intersect -a stdin -b C30913D_BULK_CNVs.txt -v | awk '{if(($3-$2)>50000) print}'
bedtools intersect -a STE0099_C_CNVs.txt -b overlapmap/ALL_Controls_multisample_merged_CNVs.bed -v | bedtools intersect -a stdin -b STE0099_BLOOD_CNVs.txt -v | awk '{if(($3-$2)>50000) print}'
bedtools intersect -a STE0100_C_CNVs.txt -b overlapmap/ALL_Controls_multisample_merged_CNVs.bed -v | bedtools intersect -a stdin -b STE0100_BLOOD_CNVs.txt -v | awk '{if(($3-$2)>50000) print}'


cat STE0100_C_ratio.txt | awk '{if($1==7 && $2>=62560000 && $2<=62570000) print}'

"""
