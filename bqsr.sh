#/usr/bin/sh

# indel realignment script for ZH-03152016-1039
# jre 1.8.0_74
# GATK v. 3.4-46


samples=$1
use_dedup=$2
wd=`pwd`
tmp_dir=${wd}/tmp
app_dir="/ifs/labs/cccb/projects/cccb/apps"
lib_dir="/ifs/labs/cccb/projects/cccb/pipelines/GATK_variant_discovery"
hg19_dir="/ifs/labs/cccb/projects/db/gatk/hg19"
hg19="${hg19_dir}/ucsc.hg19.fasta"
dbsnp="${hg19_dir}/dbsnp_137.hg19.vcf"
indels="${hg19_dir}/Mills_and_1000G_gold_standard.indels.hg19.vcf"
gatk="${app_dir}/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar"
java="${app_dir}/jre1.8.0_74/bin/java"
samstat="${app_dir}/samstat-1.5/samstat"
parallel="${lib_dir}/parallel"
PATH=$PATH:${app_dir}/samtools-1.2-5/
if [ ! -d $tmp_dir ]; then
    mkdir $tmp_dir;
fi;


##############################
# 
# Skip -known; results in slower computation time though
##############################

# Create directories
# Create indel intervals
# Realign indels
for f in `awk '{print $1}' $samples`; do
    if [ $use_dedup == "true" ]; then
        merged_bam="${wd}/Sample_${f}/bwa_mem_align/${f}.sort.dedup.bam";
        indel_dir="${wd}/Sample_${f}/gatk_indel_realn";
        bqsr_dir="${wd}/Sample_${f}/gatk_bqsr";
        indel_bam="${indel_dir}/${f}.sort.dedup.indel_realn.bam";
        bqsr_bam="${bqsr_dir}/${f}.sort.dedup.indel_realn.bqsr.bam";
    else
        merged_bam="${wd}/Sample_${f}/bwa_mem_align/${f}.sort.bam";
        indel_dir="${wd}/Sample_${f}/gatk_indel_realn";
        bqsr_dir="${wd}/Sample_${f}/gatk_bqsr";
        indel_bam="${indel_dir}/${f}.sort.indel_realn.bam";
        bqsr_bam="${bqsr_dir}/${f}.sort.indel_realn.bqsr.bam";
    fi;
    if [ ! -d $bqsr_dir ]; then
        mkdir $bqsr_dir;
    fi;
    #if [ ! -f ${bqsr_dir}/recal_data.table ]; then
	#    $java -jar $gatk \
	#        -T BaseRecalibrator \
	#        -R $hg19 \
	#        -I $indel_bam \
	#        -knownSites $dbsnp \
	#        -knownSites $indels \
	#        -o ${bqsr_dir}/recal_data.table;
    #fi;
    #if [ ! -f ${bqsr_dir}/post_recal_data.table ]; then
	#    $java -jar $gatk \
	#        -T BaseRecalibrator \
	#        -R $hg19 \
	#        -I $indel_bam \
	#        -knownSites $dbsnp \
	#        -knownSites $indels \
	#        -BQSR ${bqsr_dir}/recal_data.table \
	#        -o ${bqsr_dir}/post_recal_data.table;
    #fi;
    #if [ ! -f ${bqsr_dir}/${f}.sort.d ]; then
	#    $java -jar $gatk \
	#        -T AnalyzeCovariates \
	#        -R $hg19 \
	#        -before ${bqsr_dir}/recal_data.table \
	#        -after ${bqsr_dir}/post_recal_data.table \
	#        -plots ${bqsr_dir}/recalibration_plots.pdf;
	#    $java -jar $gatk \
	#        -T PrintReads \
	#        -R $hg19 \
	#        -I $indel_bam \
	#        -BQSR ${bqsr_dir}/recal_data.table \
	#        -o $bqsr_bam;
    #    echo $bqsr_bam >> ${wd}/bqsrRecalBam.list;
    #    $samstat $bqsr_bam;
    #fi;
    #echo $bqsr_bam;
    echo $bqsr_bam >> ${wd}/bqsrRecalBam.list;
    echo "$java -jar $gatk -T BaseRecalibrator -R $hg19 -I $indel_bam -knownSites $dbsnp -knownSites $indels -o ${bqsr_dir}/recal_data.table; $java -jar $gatk -T BaseRecalibrator -R $hg19 -I $indel_bam -knownSites $dbsnp -knownSites $indels -BQSR ${bqsr_dir}/recal_data.table -o ${bqsr_dir}/post_recal_data.table; $java -jar $gatk -T AnalyzeCovariates -R $hg19 -before ${bqsr_dir}/recal_data.table -after ${bqsr_dir}/post_recal_data.table -plots ${bqsr_dir}/recalibration_plots.pdf; $java -jar $gatk -T PrintReads -R $hg19 -I $indel_bam -BQSR ${bqsr_dir}/recal_data.table -o $bqsr_bam;" >> cmds.txt;
    $parallel -j10 < cmds.txt;
done;
for f in `awk '{print $1}' $samples`; do
    if [ $use_dedup == "true" ]; then
        merged_bam="${wd}/Sample_${f}/bwa_mem_align/${f}.sort.dedup.bam";
        indel_dir="${wd}/Sample_${f}/gatk_indel_realn";
        bqsr_dir="${wd}/Sample_${f}/gatk_bqsr";
        indel_bam="${indel_dir}/${f}.sort.dedup.indel_realn.bam";
        bqsr_bam="${bqsr_dir}/${f}.sort.dedup.indel_realn.bqsr.bam";
    else
        merged_bam="${wd}/Sample_${f}/bwa_mem_align/${f}.sort.bam";
        indel_dir="${wd}/Sample_${f}/gatk_indel_realn";
        bqsr_dir="${wd}/Sample_${f}/gatk_bqsr";
        indel_bam="${indel_dir}/${f}.sort.indel_realn.bam";
        bqsr_bam="${bqsr_dir}/${f}.sort.indel_realn.bqsr.bam";
    fi;
    $samstat $bqsr_bam;
    echo "";
done;
