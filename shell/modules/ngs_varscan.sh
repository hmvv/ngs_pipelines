#!/bin/bash

heme_varscan()
{
  bam=$1
  outDir=$2
  fileName=${bam##*/}
  sample=${fileName%%.*}

  ref="/home/doc/ref/ref_genome/ucsc.hg19.fasta"

  ##Varscan call SNPs##
  /opt/samtools/samtools mpileup -f "$ref" "$bam" | \
  java -jar /opt/varscan/VarScan.v2.3.9.jar pileup2snp > "$outDir"/"$sample".snp.txt

  ##Varscan call indels##
  /opt/samtools/samtools mpileup -f "$ref" "$bam" | \
  java -jar /opt/varscan/VarScan.v2.3.9.jar pileup2indel > "$outDir"/"$sample".indel.txt
}

germline_varscan()
{
  bam=$1
  outDir=$2
  fileName=${bam##*/}
  sample=${fileName%%.*}

  ref="/home/doc/ref/ref_genome/ucsc.hg19.fasta"

  ##Varscan call SNPs##
  /opt/samtools/samtools mpileup -f "$ref" "$bam" | \
  java -jar /opt/varscan/VarScan.v2.3.9.jar pileup2snp  --p-value 0.01 --min-avg-qual 30 > "$outDir"/"$sample".snp.txt

  ##Varscan call indels##
  /opt/samtools/samtools mpileup -f "$ref" "$bam" | \
  java -jar /opt/varscan/VarScan.v2.3.9.jar pileup2indel --p-value 0.01 --min-avg-qual 30 > "$outDir"/"$sample".indel.txt
}

tmb_varscan()
{
  SAMPLE=$1
  REF=$2
  NORMAL_BAM=$3
  TUMOR_BAM=$4
  OUT_DIR=$5
  ENV=$6
  DEPTH=$7
  NALF=$8
  TALF=$9


  echo " starting varscan:
  REF : $REF
  NORMAL_BAM :$NORMAL_BAM
  TUMOR_BAM :$TUMOR_BAM
  OUT_DIR :$OUT_DIR
  "

  java -jar /opt/varscan/VarScan.v2.3.9.jar somatic \
            <(/opt/samtools19/bin/samtools  mpileup  -f $REF  $NORMAL_BAM ) \
            <(/opt/samtools19/bin/samtools  mpileup  -f $REF  $TUMOR_BAM) \
            ${OUT_DIR}${SAMPLE}.varscan.output \
            --min-avg-qual 30 \
            --min-var-freq 0.1

}
