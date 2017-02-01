#/usr/bin/sh

# indel realignment script for ZH-03152016-1039
# jre 1.8.0_74
# GATK v. 3.4-46


samples=$1
use_dedup=$2
wd=`pwd`
tmp_dir=${wd}/tmp
app_dir="/ifs/labs/cccb/projects/cccb/apps"
hg19="/ifs/labs/cccb/projects/db/gatk/hg19/ucsc.hg19.fasta"
gatk="${app_dir}/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar"
java="${app_dir}/jre1.8.0_74/bin/java"
lib_dir="/ifs/labs/cccb/projects/cccb/pipelines/GATK_variant_discovery"
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
        merged_bam="${wd}/Sample_${f}/bwa_mem_align/${f}.sort.dedup.bam"
        indel_dir="${wd}/Sample_${f}/gatk_indel_realn"
        intervals="${f}_indel.intervals"
        if [ ! -d $indel_dir ]; then
            mkdir $indel_dir;
        fi;
        if [ ! -f ${indel_dir}/${f%.bam}_indel.intervals ]; then
            $java -jar $gatk \
              -T RealignerTargetCreator \
              -R $hg19 \
              -I $merged_bam \
              -nt 10 \
              -o ${indel_dir}/${intervals};
        fi;
    else
        merged_bam="${wd}/Sample_${f}/bwa_mem_align/${f}.sort.bam"
        indel_dir="${wd}/Sample_${f}/gatk_indel_realn"
        intervals="${f}_indel.intervals"
        if [ ! -d $indel_dir ]; then
            mkdir $indel_dir;
        fi;
        if [ ! -f ${indel_dir}/${f%.bam}_indel.intervals ]; then
            $java -jar $gatk \
              -T RealignerTargetCreator \
              -R $hg19 \
              -I $merged_bam \
              -nt 10 \
              -o ${indel_dir}/${intervals};
        fi;
    fi;
done;

for f in `awk '{print $1}' $samples`; do
    if [ $use_dedup == "true" ]; then
        merged_bam="${wd}/Sample_${f}/bwa_mem_align/${f}.sort.dedup.bam";
        indel_dir="${wd}/Sample_${f}/gatk_indel_realn";
        indel_bam="${indel_dir}/${f}.sort.dedup.indel_realn.bam";
        intervals="${f}_indel.intervals";
    else
        merged_bam="${wd}/Sample_${f}/bwa_mem_align/${f}.sort.bam";
        indel_dir="${wd}/Sample_${f}/gatk_indel_realn";
        indel_bam="${indel_dir}/${f}.sort.indel_realn.bam";
        intervals="${f}_indel.intervals";
    fi;
    if [ ! -f $indel_bam ]; then
        echo "$java -jar $gatk \
          -T IndelRealigner \
          -R $hg19 \
          -targetIntervals ${indel_dir}/${intervals} \
          -I $merged_bam \
          -o $indel_bam" >> tmp_cmds.txt;
    fi;
done;
${lib_dir}/parallel -j10 < tmp_cmds.txt;
rm tmp_cmds.txt;

