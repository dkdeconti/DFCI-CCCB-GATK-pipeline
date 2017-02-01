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
samtools="${app_dir}/samtools-1.2-5/samtools"
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
    header=${bqsr_bam%.bam}.header.txt;
    echo $bqsr_bam
    if [ ! -f ${bqsr_bam}.bak ]; then
        cp $bqsr_bam ${bqsr_bam}.bak;
    fi;
    $samtools view -H $bqsr_bam > $header;
    python ${lib_dir}/fix_headers.py $header > ${header}.fixed;
    $samtools view ${bqsr_bam}.bak | cat ${header}.fixed - | $samtools view -bt $hg19 -o $bqsr_bam -;
done;

