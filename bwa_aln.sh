#/usr/bin/sh

# alignment script for ZH-03152016-1039
# BWA MEM v. 0.7.12
# jre 1.8.0_74
# picardtools 2.1.1


samples=$1
sample_lanes=$2
is_paired=$3
wd=`pwd`
tmp_dir=${wd}/tmp
hg19="/ifs/labs/cccb/projects/db/gatk/hg19/ucsc.hg19.fasta"
bwa="/ifs/labs/cccb/projects/cccb/apps/bwa-0.7.12/bwa"
java="/ifs/labs/cccb/projects/cccb/apps/jre1.8.0_74/bin/java"
picard="/ifs/labs/cccb/projects/cccb/apps/picard-tools-2.1.1/picard.jar"
samstat="/ifs/labs/cccb/projects/cccb/apps/samstat-1.5/samstat"
lib_dir="/ifs/labs/cccb/projects/cccb/pipelines/GATK_variant_discovery"
PATH=$PATH:/ifs/labs/cccb/projects/cccb/apps/samtools-1.2-5/
if [ ! -d $tmp_dir ]; then
    mkdir $tmp_dir;
fi;

##############################
# initial alignment
##############################

# Create directories
for f in `awk '{print $1}' $samples`; do
    if [ ! -d ${wd}/Sample_${f}/bwa_mem_align ]; then
        mkdir ${wd}/Sample_${f}/bwa_mem_align;
    fi;
    if [ ! -d ${wd}/Sample_${f}/lane_specific_fastq ]; then
        mkdir ${wd}/Sample_${f}/lane_specific_fastq;
    fi;
    if [ ! -d ${wd}/Sample_${f}/lane_specific_align ]; then
        mkdir ${wd}/Sample_${f}/lane_specific_align;
    fi;
done;
# symlink the lane specific fastq (for GATK RG purposes)
while read line; do
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
done < $sample_lanes;
# align with bwa
# added RG group via bwa
for sample in `awk '{print $1}' $samples`; do
    for read1 in ${wd}/Sample_${sample}/lane_specific_fastq/*R1*; do
        read2="${read1%_R1_*}_R2_${read1##*_R1_}";
        sam_file=${read1%_R1_*};
        sam_file="./Sample_${sample}/lane_specific_align/${sam_file##*/}.sam";
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
done;

##############################
# convert to BAM and merge and sort
##############################

# convert to BAM with picardtools
for sample in `awk '{print $1}' $samples`; do
    for read1 in ${wd}/Sample_${sample}/lane_specific_fastq/*R1*; do
        read2="${read1%_R1_*}_R2_${read1##*_R1_}";
        sam_file=${read1%_R1_*};
        sam_file="${wd}/Sample_${sample}/lane_specific_align/${sam_file##*/}.sam";
        bam_file="${sam_file%.sam}.bam";
        sort_file="${sam_file%.sam}.sort.bam";
        SAMPLE=`python -c \
          "s='${read1##*/}'; print '_'.join(s.strip('\n').split('_')[0:3])"`;
        if [ ! -f $bam_file ] && [ ! -f $sort_file ]; then
            echo "$java -jar $picard SamFormatConverter \
              INPUT=${sam_file} OUTPUT=${bam_file}" >> tmp_cmds.txt;
            echo "rm $sam_file" >> rm_sam.sh;
        fi;
    done;
done;
${lib_dir}/parallel -j10 < tmp_cmds.txt;
rm tmp_cmds.txt;
#sh rm_sam.sh; rm rm_sam.sh;


# sort Bam with picardtools
for sample in `awk '{print $1}' $samples`; do
    for read1 in ${wd}/Sample_${sample}/lane_specific_fastq/*R1*; do
        read2="${read1%_R1_*}_R2_${read1##*_R1_}";
        sam_file=${read1%_R1_*};
        sam_file="${wd}/Sample_${sample}/lane_specific_align/${sam_file##*/}.sam";
        bam_file="${sam_file%.sam}.bam";
        sort_file="${sam_file%.sam}.sort.bam";
        SAMPLE=`python -c \
          "s='${read1##*/}'; print '_'.join(s.strip('\n').split('_')[0:3])"`;
        if [ ! -f $sort_file ]; then
            echo "$java -jar $picard SortSam \
              INPUT=${bam_file} OUTPUT=${sort_file} SORT_ORDER=coordinate \
              TMP_DIR=$tmp_dir" >> tmp_cmds.txt;
        fi;
    done;
done;
${lib_dir}/parallel -j10 < tmp_cmds.txt;
rm tmp_cmds.txt;

# Index sorted bam with picardtools
for sample in `awk '{print $1}' $samples`; do
    for read1 in ${wd}/Sample_${sample}/lane_specific_fastq/*R1*; do
        read2="${read1%_R1_*}_R2_${read1##*_R1_}";
        sam_file=${read1%_R1_*};
        sam_file="${wd}/Sample_${sample}/lane_specific_align/${sam_file##*/}.sam";
        bam_file="${sam_file%.sam}.bam";
        sort_file="${sam_file%.sam}.sort.bam";
        SAMPLE=`python -c \
          "s='${read1##*/}'; print '_'.join(s.strip('\n').split('_')[0:3])"`;
        if [ ! -f $sort_file ]; then
            echo "$java -jar $picard BuildBamIndex \
              INPUT=${sort_file} TMP_DIR=$tmp_dir" >> tmp_cmds.txt;
            echo "rm $bam_file" >> rm_bam.sh;
        fi;
        #rm $sam_file
    done;
done;
${lib_dir}/parallel -j10 < tmp_cmds.txt;
rm tmp_cmds.txt;
#sh rm_bam.sh; rm rm_bam.sh;

# merge BAM with picardtools
for sample in `awk '{print $1}' $samples`; do
    merged_bam="${wd}/Sample_${sample}/bwa_mem_align/${sample}.sort.bam"
    INPUTS="";
    sort_list=`ls ${wd}/Sample_${sample}/lane_specific_align/*.sort.bam`
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
done;


##############################
# dedup (mark duplicates) with picardtools
##############################
for sample in `awk '{print $1}' $samples`; do
    merged_bam="${wd}/Sample_${sample}/bwa_mem_align/${sample}.sort.bam"
    dedup_bam="${wd}/Sample_${sample}/bwa_mem_align/${sample}.sort.dedup.bam"
    metrics_file="${wd}/Sample_${sample}/bwa_mem_align/${sample}.dedup_metrics.txt"
    if [ ! -f $dedup_bam ]; then
        echo "$java -jar $picard MarkDuplicates \
          INPUT=${merged_bam} \
          OUTPUT=${dedup_bam} \
          METRICS_FILE=${metrics_file} \
          CREATE_INDEX=true \
          OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
          TMP_DIR=$tmp_dir" >> tmp_cmds.txt;
    fi;
done;
if [ -f tmp_cmds.txt ]; then
    ${lib_dir}/parallel -j10 < tmp_cmds.txt;
fi;
rm tmp_cmds.txt;

##############################
# Samstat summary statistic report for Dedup bam
##############################
for sample in `awk '{print $1}' $samples`; do
    merged_bam="${wd}/Sample_${sample}/bwa_mem_align/${sample}.sort.bam"
    dedup_bam="${wd}/Sample_${sample}/bwa_mem_align/${sample}.sort.dedup.bam"
    metrics_file="${wd}/Sample_${sample}/bwa_mem_align/${sample}.dedup_metrics.txt"
    if [ -f $dedup_bam ] && [ ! -f ${dedup_bam}.samstat.html ]; then
        echo "$samstat $dedup_bam" >> tmp_cmds.txt;
    fi;
done;
${lib_dir}/parallel -j10 < tmp_cmds.txt;
rm tmp_cmds.txt;
#rm -r $tmp_dir;

