#/usr/bin/sh

# indel realignment script for ZH-03152016-1039
# jre 1.8.0_74
# GATK v. 3.4-46


samples=$1
probe_list=$2
use_dedup=$3 # true
wd=`pwd`
tmp_dir=${wd}/tmp
if [ $use_dedup == "true" ]; then
    coverage_dir=${wd}/coverage;
else
    coverage_dir=${wd}/coverage_raw;
fi;
app_dir="/ifs/labs/cccb/projects/cccb/apps"
hg19_dir="/ifs/labs/cccb/projects/db/gatk/hg19"
hg19="${hg19_dir}/ucsc.hg19.fasta"
dbsnp="${hg19_dir}/dbsnp_137.hg19.vcf"
indels="${hg19_dir}/Mills_and_1000G_gold_standard.indels.hg19.vcf"
gatk="${app_dir}/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar"
java="${app_dir}/jre1.8.0_74/bin/java"
lib_dir="/ifs/labs/cccb/projects/cccb/pipelines/GATK_variant_discovery"
if [ ! -d $tmp_dir ]; then
    mkdir $tmp_dir;
fi;
if [ ! -d $coverage_dir ]; then
    mkdir $coverage_dir;
fi;

# checks if bqsr recal bam list is present in main directory
if [ ! -d ${wd}/recal_bam.list ]; then
    rm ${wd}/recal_bam.list;
fi;
#bam_list="${wd}/bqsrRecalBam.list"
for sample in `awk '{print $1}' $samples`; do
    bam_dir="${wd}/Sample_${sample}/bwa_mem_align"
    if [ $use_dedup == "true" ]; then
        merged_bam="${bam_dir}/${sample}.sort.dedup.indel_realn.bqsr.bam";
    else
        merged_bam="${bam_dir}/${sample}.sort.indel_realn.bqsr.bam";
    fi;
    echo $merged_bam >> ${wd}/recal_bam.list
done;
bam_list="${wd}/recal_bam.list"

# Runs GATK DepthOfCoverage module
if [ ! -f ${coverage_dir}/coverage.sample_statistics ]; then
    $java -jar $gatk \
      -T DepthOfCoverage \
      -R $hg19 \
      -L $probe_list \
      -I $bam_list \
      --partitionType library \
      -nt 10 \
      --omitDepthOutputAtEachBase --omitLocusTable --omitIntervalStatistics \
      --outputFormat rtable \
      --out ${coverage_dir}/coverage;
fi;

# Runs QC
Rscript ${lib_dir}/depthofcoverageQC.R $coverage_dir

