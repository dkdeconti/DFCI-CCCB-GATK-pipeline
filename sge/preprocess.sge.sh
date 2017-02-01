#!/usr/bin/bash

f=$1;
sample_lanes=$2;
wd=$3; # not used
use_dedup=$4; # true

tmp_dir=${wd}/tmp
app_dir="/ifs/labs/cccb/projects/cccb/apps"
hg19_dir="/ifs/labs/cccb/projects/db/gatk/hg19"
hg19="${hg19_dir}/ucsc.hg19.fasta"
dbsnp="${hg19_dir}/dbsnp_137.hg19.vcf"
indels="${hg19_dir}/Mills_and_1000G_gold_standard.indels.hg19.vcf"
bwa="${app_dir}/bwa-0.7.12/bwa"
gatk="${app_dir}/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar"
java="${app_dir}/jre1.8.0_74/bin/java"
picard="/ifs/labs/cccb/projects/cccb/apps/picard-tools-2.1.1/picard.jar"
samstat="${app_dir}/samstat-1.5/samstat"
PATH=$PATH:${app_dir}/samtools-1.2-5/
if [ ! -d $tmp_dir ]; then
    mkdir $tmp_dir;
fi;

# Note: requires hard coded path to working directory until fix is created
wd="/ifs/labs/cccb/projects/cccb/projects/2016/5/ZH_05172016_1085"



##############################
# initial alignment
##############################

# Create directories
if [ ! -d ${wd}/Sample_${f}/bwa_mem_align ]; then
    mkdir ${wd}/Sample_${f}/bwa_mem_align;
fi;
if [ ! -d ${wd}/Sample_${f}/lane_specific_fastq ]; then
    mkdir ${wd}/Sample_${f}/lane_specific_fastq;
fi;
if [ ! -d ${wd}/Sample_${f}/lane_specific_align ]; then
    mkdir ${wd}/Sample_${f}/lane_specific_align;
fi;
# symlink the lane specific fastq (for GATK RG purposes)
for line in `grep $f $sample_lanes`; do
    sample_name=`python -c "s='${line}'; print s.strip('\n').split()[0]"`;
    lane_fastq=`python -c "s='${line}'; print s.strip('\n').split()[1]"`;
    target_fastq=${lane_fastq##*/};
    target_fastq=${target_fastq//-/_};
    target_dir=${wd}/Sample_${sample_name}/lane_specific_fastq;
    target_fastq=${target_dir}/${target_fastq};
    if [ ! -f $target_fastq ]; then
        echo $lane_fastq $target_fastq;
        ln -s $lane_fastq $target_fastq;
    fi;
done;
# align with bwa
# added RG group via bwa
for read1 in ${wd}/Sample_${f}/lane_specific_fastq/*R1*; do
    read2="${read1%_R1_*}_R2_${read1##*_R1_}";
    sam_file=${read1%_R1_*};
    sam_file="./Sample_${f}/lane_specific_align/${sam_file##*/}.sam";
    sort_file="${sam_file%.sam}.sort.bam"
    SAMPLE=`python -c \
      "s='${read1##*/}'; print '_'.join(s.strip('\n').split('_')[0:-4])"`;
    ID=`python -c \
      "s='${read1##*/}'; print '_'.join(s.strip('\n').split('_')[0:-1])"`;
    RG='@RG\tID:'$ID'\tPL:ILLUMINA\tLB:'$SAMPLE'\tSM:'$SAMPLE'\tCN:CCCB';
    if [ ! -f $sam_file ] && [ ! -f $sort_file ]; then
        if [ $is_paired == "true" ]; then
            $bwa mem -t 10 -R "$RG" $hg19 $read1 $read2 > $sam_file;
        else
            $bwa mem -t 10 -R "$RG" $hg19 $read1 > $sam_file;
        fi;
    fi;
done;

##############################
# convert to BAM and merge and sort
##############################

# convert to BAM with picardtools
for read1 in ${wd}/Sample_${f}/lane_specific_fastq/*R1*; do
    read2="${read1%_R1_*}_R2_${read1##*_R1_}";
    sam_file=${read1%_R1_*};
    sam_file="${wd}/Sample_${f}/lane_specific_align/${sam_file##*/}.sam";
    bam_file="${sam_file%.sam}.bam";
    sort_file="${sam_file%.sam}.sort.bam";
    SAMPLE=`python -c \
      "s='${read1##*/}'; print '_'.join(s.strip('\n').split('_')[0:3])"`;
    if [ ! -f $bam_file ] && [ ! -f $sort_file ]; then
        $java -jar $picard SamFormatConverter \
          INPUT=${sam_file} OUTPUT=${bam_file};
    fi;
done;


# sort Bam with picardtools
for read1 in ${wd}/Sample_${f}/lane_specific_fastq/*R1*; do
    read2="${read1%_R1_*}_R2_${read1##*_R1_}";
    sam_file=${read1%_R1_*};
    sam_file="${wd}/Sample_${f}/lane_specific_align/${sam_file##*/}.sam";
    bam_file="${sam_file%.sam}.bam";
    sort_file="${sam_file%.sam}.sort.bam";
    SAMPLE=`python -c \
      "s='${read1##*/}'; print '_'.join(s.strip('\n').split('_')[0:3])"`;
    if [ ! -f $sort_file ]; then
        $java -jar $picard SortSam \
          INPUT=${bam_file} OUTPUT=${sort_file} SORT_ORDER=coordinate \
          TMP_DIR=$tmp_dir;
    fi;
done;

# Index sorted bam with picardtools
for read1 in ${wd}/Sample_${f}/lane_specific_fastq/*R1*; do
    read2="${read1%_R1_*}_R2_${read1##*_R1_}";
    sam_file=${read1%_R1_*};
    sam_file="${wd}/Sample_${f}/lane_specific_align/${sam_file##*/}.sam";
    bam_file="${sam_file%.sam}.bam";
    sort_file="${sam_file%.sam}.sort.bam";
    SAMPLE=`python -c \
      "s='${read1##*/}'; print '_'.join(s.strip('\n').split('_')[0:3])"`;
    if [ ! -f $sort_file ]; then
        $java -jar $picard BuildBamIndex \
          INPUT=${sort_file} TMP_DIR=$tmp_dir;
    fi;
done;

# merge BAM with picardtools
merged_bam="${wd}/Sample_${f}/bwa_mem_align/${f}.sort.bam"
INPUTS="";
sort_list=`ls ${wd}/Sample_${f}/lane_specific_align/*.sort.bam`
for sort_file in $sort_list; do
    INPUTS="INPUT=${sort_file} $INPUTS";
done;
# sort built into merge
if [ ! -f $merged_bam ]; then
    $java -jar $picard MergeSamFiles \
      ${INPUTS} OUTPUT=${merged_bam} SORT_ORDER=coordinate \
      USE_THREADING=true TMP_DIR=$tmp_dir;
    $java -jar $picard BuildBamIndex INPUT=${merged_bam} TMP_DIR=$tmp_dir;
fi;


##############################
# dedup (mark duplicates) with picardtools
##############################

merged_bam="${wd}/Sample_${f}/bwa_mem_align/${f}.sort.bam"
dedup_bam="${wd}/Sample_${f}/bwa_mem_align/${f}.sort.dedup.bam"
metrics_file="${wd}/Sample_${f}/bwa_mem_align/${f}.dedup_metrics.txt"
if [ ! -f $dedup_bam ]; then
    echo "$java -jar $picard MarkDuplicates \
      INPUT=${merged_bam} \
      OUTPUT=${dedup_bam} \
      METRICS_FILE=${metrics_file} \
      CREATE_INDEX=true \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      TMP_DIR=$tmp_dir" >> tmp_cmds.txt;
fi;


##############################
# Samstat summary statistic report for Dedup bam
##############################
merged_bam="${wd}/Sample_${f}/bwa_mem_align/${f}.sort.bam"
dedup_bam="${wd}/Sample_${f}/bwa_mem_align/${f}.sort.dedup.bam"
metrics_file="${wd}/Sample_${f}/bwa_mem_align/${f}.dedup_metrics.txt"
if [ -f $dedup_bam ] && [ ! -f ${dedup_bam}.samstat.html ]; then
    echo "$samstat $dedup_bam" >> tmp_cmds.txt;
fi;



##############################
# Indel realignment
##############################


# Create directories
# Create indel intervals
# Realign indels


if [ $use_dedup == "true" ]; then
    merged_bam="${wd}/Sample_${f}/bwa_mem_align/${f}.sort.dedup.bam"
    indel_dir="${wd}/Sample_${f}/gatk_indel_realn"
    intervals="${f}_indel.intervals"
    if [ ! -d $indel_dir ]; then
        mkdir $indel_dir;
    fi;
else
    merged_bam="${wd}/Sample_${f}/bwa_mem_align/${f}.sort.bam"
    indel_dir="${wd}/Sample_${f}/gatk_indel_realn"
    intervals="${f}_indel.intervals"
    if [ ! -d $indel_dir ]; then
        mkdir $indel_dir;
    fi;
fi;
$java -jar $gatk \
    -T RealignerTargetCreator \
    -R $hg19 \
    -I $merged_bam \
    -nt 10 \
    -o ${indel_dir}/${intervals};


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
$java -jar $gatk \
    -T IndelRealigner \
    -R $hg19 \
    -targetIntervals ${indel_dir}/${intervals} \
    -I $merged_bam \
    -o $indel_bam;


##############################
# BQSR
# Skip -known; results in slower computation time though
##############################



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
mkdir $bqsr_dir;
$java -jar $gatk \
    -T BaseRecalibrator \
    -R $hg19 \
    -I $indel_bam \
    -knownSites $dbsnp \
    -knownSites $indels \
    -o ${bqsr_dir}/recal_data.table;
$java -jar $gatk \
    -T BaseRecalibrator \
    -R $hg19 \
    -I $indel_bam \
    -knownSites $dbsnp \
    -knownSites $indels \
    -BQSR ${bqsr_dir}/recal_data.table \
    -o ${bqsr_dir}/post_recal_data.table;
$java -jar $gatk \
    -T AnalyzeCovariates \
    -R $hg19 \
    -before ${bqsr_dir}/recal_data.table \
    -after ${bqsr_dir}/post_recal_data.table \
    -plots ${bqsr_dir}/recalibration_plots.pdf;
$java -jar $gatk \
    -T PrintReads \
    -R $hg19 \
    -I $indel_bam \
    -BQSR ${bqsr_dir}/recal_data.table \
    -o $bqsr_bam;
$samstat $bqsr_bam;
