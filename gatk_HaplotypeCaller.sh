#/usr/bin/sh

# jre 1.8.0_74
# GATK v. 3.4-46


samples=$1
use_dedup=$2 # true
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
parallel="${lib_dir}/parallel"

if [ -f ${wd}/cmds.txt ]; then
    rm ${wd}/cmds.txt;
fi;

for f in `awk '{print $1}' $samples`; do
    bqsr_dir="${wd}/Sample_${f}/gatk_bqsr";
    vcf_dir="${wd}/Sample_${f}/gatk_vcf";
    if [ ! -d $vcf_dir ]; then
        mkdir $vcf_dir;
    fi;
    if [ $use_dedup == "true" ]; then
        bqsr_bam="${bqsr_dir}/${f}.sort.dedup.indel_realn.bqsr.bam";
        vcf="${vcf_dir}/${f}.dedup.indel_realn.bqsr.vcf";
    else
        bqsr_bam="${bqsr_dir}/${f}.sort.indel_realn.bqsr.bam"
        vcf="${vcf_dir}/${f}.indel_realn.bqsr.vcf"
    fi;
    echo "$java -jar $gatk -T HaplotypeCaller -R $hg19 -I $bqsr_bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 --dbsnp $dbsnp --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o $vcf" >> ${wd}/cmds.txt;
done;
$parallel -j5 < ${wd}/cmds.txt;
