#!/bin/bash

vep_83()
{
  /opt/perl/bin/perl /opt/vep/ensembl-tools-release-83/scripts/variant_effect_predictor/variant_effect_predictor.pl \
  -i $1 \
  -o $2 \
  --offline \
  --dir_cache /opt/vep/ensembl-tools-release-83/cache/ \
  --sift p \
  --polyphen p \
  --hgvs \
  --symbol \
  --vcf \
  --pubmed \
  --fasta /home/doc/ref/ref_genome/Homo_sapiens.GRCh37.73.dna_sm.primary_assembly.fa \
  --pick_allele \
  --check_alleles \
  --force_overwrite \
  --gmaf \
  --maf_1kg
}

vep_94_panel()
{

  /opt/vep_94/ensembl-tools-release-94/vep_94/ensembl-vep/vep \
  -i $1 \
  -o $2 \
  --offline \
  --dir_cache /opt/vep_94/ensembl-tools-release-94/cache \
  --sift p \
  --polyphen p \
  --hgvs \
  --symbol \
  --vcf \
  --pubmed \
  --fasta /home/doc/ref/ref_genome/Homo_sapiens.GRCh37.73.dna_sm.primary_assembly.fa \
  --pick_allele \
  --force_overwrite \
  --af \
  --af_1kg

}

vep_94_tmb()
{
  /opt/vep_94/ensembl-tools-release-94/vep_94/ensembl-vep/vep \
	-i $1 \
	-o $2 \
	--offline \
	--dir_cache /opt/vep_94/ensembl-tools-release-94/cache \
	--vcf \
	--refseq \
	--pick_allele \
	--sift p \
	--polyphen p \
	--hgvs \
	--symbol \
	--vcf \
	--pubmed \
  --force_overwrite \
	--fasta /home/doc/ref/ref_genome/Homo_sapiens.GRCh37.73.dna_sm.primary_assembly.fa
}
