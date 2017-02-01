#/usr/bin/sh

# indel realignment script for ZH-03152016-1039
# jre 1.8.0_74
# GATK v. 


samples=$1
sample_lanes=$2
probe_list=$3 # as BED file
is_paired=$4 # true
use_dedup=$5 # true
wd=`pwd`
tmp_dir=${wd}/tmp
app_dir="/ifs/labs/cccb/projects/cccb/apps"
hg19_dir="/ifs/labs/cccb/projects/db/gatk/hg19"
hg19="${hg19_dir}/ucsc.hg19.fasta"
dbsnp="${hg19_dir}/dbsnp_137.hg19.vcf"
indels="${hg19_dir}/Mills_and_1000G_gold_standard.indels.hg19.vcf"
gatk="${app_dir}/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar"
java="${app_dir}/jre1.8.0_74/bin/java"

lib_dir="/ifs/labs/cccb/projects/cccb/pipelines/GATK_variant_discovery"

# RG assignment, alignment, and dedup
${lib_dir}/bwa_aln.sh $samples $sample_lanes $is_paired

# Realign indels
${lib_dir}/indel_realn.sh $samples $use_dedup

# Base quality score recalibration
${lib_dir}/bqsr.sh $samples $use_dedup

# Depth of coverage QC
${lib_dir}/depthofcoverageQC.sh $samples $probe_list $use_dedup
