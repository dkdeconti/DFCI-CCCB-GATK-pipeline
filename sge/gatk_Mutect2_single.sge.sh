#/usr/bin/sh

# jre 1.8.0_74
# GATK v. 3.4-46


f=$1
wd=$2
use_dedup=$3 # true
tmp_dir=${wd}/tmp
app_dir="/ifs/labs/cccb/projects/cccb/apps"
lib_dir="/ifs/labs/cccb/projects/cccb/pipelines/GATK_variant_discovery"
hg19_dir="/ifs/labs/cccb/projects/db/gatk/hg19"
hg19="${hg19_dir}/ucsc.hg19.fasta"
cosmic="${hg19_dir}/cosmic_combined.vcf"
dbsnp="${hg19_dir}/dbsnp_137.hg19.vcf"
indels="${hg19_dir}/Mills_and_1000G_gold_standard.indels.hg19.vcf"
gatk="${app_dir}/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar"
java="${app_dir}/jre1.8.0_74/bin/java"

bqsr_dir="${wd}/Sample_${f}/gatk_bqsr";
vcf_dir="${wd}/Sample_${f}/gatk_vcf";
if [ ! -d $vcf_dir ]; then
    mkdir $vcf_dir;
fi;
if [ $use_dedup == "true" ]; then
    bqsr_bam="${bqsr_dir}/${f}.sort.dedup.indel_realn.bqsr.bam";
    vcf="${vcf_dir}/${f}.dedup.indel_realn.bqsr.MuTect2Single.vcf";
else
    bqsr_bam="${bqsr_dir}/${f}.sort.indel_realn.bqsr.bam"
    vcf="${vcf_dir}/${f}.indel_realn.bqsr.MuTect2Single.vcf"
fi;
$java -jar $gatk \
    -T MuTect2 \
    -R $hg19 \
    -I:tumor $bqsr_bam \
    --genotyping_mode DISCOVERY \
    -stand_emit_conf 10 \
    -stand_call_conf 30 \
    --dbsnp $dbsnp \
    --cosmic $cosmic \
    --emitRefConfidence GVCF \
    --variant_index_type LINEAR \
    --variant_index_parameter 128000 \
    -o $vcf
