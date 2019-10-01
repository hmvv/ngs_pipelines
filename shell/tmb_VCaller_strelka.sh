#!/usr/bin/env bash

#################################################
# Strelka
#################################################

SAMPLE=$1
REF=$2
NORMAL_BAM=$3
TUMOR_BAM=$4
OUT_DIR=$5
ENV=$6
DEPTH=$7
NALF=$8
TALF=$9

echo " starting strelka:
REF : $REF
NORMAL_BAM : $NORMAL_BAM
TUMOR_BAM : $TUMOR_BAM
OUT_DIR : $OUT_DIR
"

# ###configuration
/opt/strelka/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam $NORMAL_BAM \
    --tumorBam $TUMOR_BAM \
    --ref $REF \
    --runDir $OUT_DIR \
    --exome

# execution on a single local machine with 8 parallel jobs
$OUT_DIR/runWorkflow.py -m local -j 8

gunzip < ${OUT_DIR}results/variants/somatic.snvs.vcf.gz > ${OUT_DIR}somatic.snvs.vcf.txt

gunzip < ${OUT_DIR}results/variants/somatic.indels.vcf.gz > ${OUT_DIR}somatic.indels.vcf.txt

java -jar /opt/GATK4/GenomeAnalysisTK.jar \
     -R $REF \
     -T VariantsToTable \
     -V ${OUT_DIR}somatic.snvs.vcf.txt \
     -F CHROM -F POS -F REF -F ALT -F FILTER \
     -F QSS -F TQSS -F NT -F QSS_NT -F TQSS_NT -F SGT -F SOMATIC -F DP -F MQ -F MQ0 -F  ReadPosRankSum -F SNVSB -F PNOISE -F PNOISE2 -F SomaticEVS \
     -GF DP -GF FDP -GF SDP -GF SUBDP -GF AU -GF CU -GF GU -GF TU \
     -o ${OUT_DIR}output.snvs.vcf.filter.txt \
     --showFiltered

java -jar /opt/GATK4/GenomeAnalysisTK.jar \
     -R $REF \
     -T VariantsToTable \
     -V ${OUT_DIR}somatic.indels.vcf.txt \
     -F CHROM -F POS -F REF -F ALT -F FILTER \
     -F QSI -F TQSI -F NT -F QSI_NT -F TQSI_NT -F SGT -F RU -F RC -F IC -F IHP -F MQ -F MQ0 -F SOMATIC -F OVERLAP -F SomaticEVS \
     -GF DP -GF DP2 -GF TAR -GF TIR -GF TOR -GF DP50 -GF FDP50 -GF SUBDP50 -GF BCN50 \
     -o ${OUT_DIR}output.indels.vcf.filter.txt \
     --showFiltered
