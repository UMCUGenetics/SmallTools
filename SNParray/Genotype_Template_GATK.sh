# SETTINGS
ref="/../../Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta"
snp="taqman_design.vcf"
gatk="/../../GenomeAnalysisTK-3.2-0/GenomeAnalysisTK.jar"
call="UnifiedGenotyper"
panel=${snp%_*}

# First argument is bamfolder
folder=$1
outfolder=$2

for bam in $folder/*.bam; do
	filename=${bam##*/}
	samplename=${filename%.*}
	outfile="${outfolder}${samplename}_${panel}.vcf"

	echo $samplename
	java -Xms2g -Xmx9g -jar ${gatk} -T ${call} -R ${ref} -L ${snp} --output_mode EMIT_ALL_SITES -I ${bam} -o ${outfile}
done
