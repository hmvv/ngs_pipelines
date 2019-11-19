#!/bin/bash

# ##############################################################################
# # functions
# ##############################################################################
display_usage()
{
cat <<EOF >> /dev/stderr

 USAGE: $0

 OPTIONS:
 	-s SAMPLE NAME
 	-f FASTQ Directory
	-o Output Directory

EOF
}

parse_options()
{
    IMPORT=0

		while getopts "hz:s:f:o:e:q:u:p:" opt; do

				case $opt in
					h)
					display_usage
					exit 1
					;;
	        z)
					IMPORT=$OPTARG
					;;
					s)
					SAMPLE=$OPTARG
					;;
					f)
					FASTQ_DIR=$OPTARG
					;;
					o)
					OUTPUT_DIR=$OPTARG
					;;
          e)
          ENVIRONMENT=$OPTARG
          ;;
          q)
          QUEUEID=$OPTARG
          ;;
          u)
          USER=$OPTARG
          ;;
          p)
          PASSWORD=$OPTARG
          ;;
					:)
					echo "Option -$OPTARG requires an argument."
					;;
					\?)
					echo "Invalid option: -$OPTARG"
		   esac
    done

    if [ $IMPORT -gt 0 ] ; then
        return 0
    fi

    return 1
}

load_modules()
{
	source /home/pipelines/ngs_${ENVIRONMENT}/shell/modules/ngs_utils.sh
  source /home/pipelines/ngs_${ENVIRONMENT}/shell/modules/ngs_align.sh

}

run_trimmomatic()
{
	log_info "Running Trimmomatic: Removing sequences < Q20 sample- $SAMPLE"
	trimmomatic="java -jar /opt/trimmomatic/Trimmomatic-0.33/trimmomatic-0.33.jar PE -phred33 -threads 8 \
	${FASTQ_DIR}${SAMPLE}_S*_R1_001.fastq.gz \
	${FASTQ_DIR}${SAMPLE}_S*_R2_001.fastq.gz \
	${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}_filt_paired_R1_001.fastq.gz \
	${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}_filt_unpaired_R1_001.fastq.gz \
	${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}_filt_paired_R2_001.fastq.gz \
	${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}_filt_unpaired_R2_001.fastq.gz \
	TRAILING:20 \
	AVGQUAL:20 \
	SLIDINGWINDOW:10:20 \
	MINLEN:30"
	($trimmomatic) 2>&1 | tee ${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.trimmomatic.summary.txt
}

run_alignment()
{
	log_info "Running bwa mem aligner: $SAMPLE"

  tmb_bwaAlign $SAMPLE  \
	$REF_GENOME  \
	$MAP_QUALITY \
	${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}_filt_paired_R1_001.fastq.gz \
	${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}_filt_paired_R2_001.fastq.gz  \
	$OUTPUT_DIR_SAMPLE_ALIGNMENT  \
	$LOG_FILE
}

sort_bam()
{
log_info "Generating sort bam by coordinate sample- $SAMPLE"
java -jar /opt/picard2/picard.jar SortSam \
          I=${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.bam  \
          O=${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.sort.bam  \
          SORT_ORDER=coordinate
}

generate_alignStat()
{
log_info "Generating alignment stat sample- $SAMPLE"
java -jar /opt/picard2/picard.jar CollectAlignmentSummaryMetrics \
          R=$REF_GENOME \
          I=${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.sort.bam \
          O=${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.sort.bam.alignmentMetrics.txt
}

remove_duplicates()
{
log_info "Removing duplicates sample- $SAMPLE "
java -jar /opt/picard2/picard.jar MarkDuplicates \
          REMOVE_DUPLICATES=true \
          INPUT=${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.sort.bam \
          OUTPUT=${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.sort.rmdups.bam \
          METRICS_FILE=${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.sort.rmdups.bam.metrics.txt
}

generate_bamIndex()
{
log_info "Generating bam index sample- $SAMPLE"
java -jar /opt/picard2/picard.jar BuildBamIndex \
          I=${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.sort.rmdups.bam
}

calculate_targetCoverage()
{
log_info "Generating CalculateHsMetrics sample- $SAMPLE "
# # java -jar /opt/picard/picard-tools-1.134/picard.jar BedToIntervalList  I=cre_v1_design.bed O=/home/hhadmin/exome_pipeline/01_bamQC/cre_v1_design_bed.interval_list SD=/doc/ref/ref_genome/ucsc.hg19.dict
java -jar /opt/picard/picard-tools-1.134/picard.jar CalculateHsMetrics \
          I=${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.sort.rmdups.bam  \
          O=${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.output_hs_metrics.txt \
          R=$REF_GENOME \
          BAIT_INTERVALS=   /home/environments/ngs_${ENVIRONMENT}/assayCommonFiles/tmbAssay/cre_design_bed.interval_list \
          TARGET_INTERVALS=   /home/environments/ngs_${ENVIRONMENT}/assayCommonFiles/tmbAssay/cre_design_bed.interval_list
}

calculate_breadthCoverage()
{
/opt/samtools19/bin/samtools depth ${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.sort.rmdups.bam  > ${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.depth.bed
awk '{ if ( ( $1 !~ "chrM") && ( $1 !~ /\_/  ) && ( $3 >= 10) ) {print $1 "\t" $2 "\t" $2+1 "\t" $3} }' ${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.depth.bed > ${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.depth.filter2.bed
/opt/bedtools2/bin/bedtools intersect -a ${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.depth.filter2.bed -b /home/environments/ngs_${ENVIRONMENT}/assayCommonFiles/tmbAssay/cre_design_ucsc_exon.txt_filter.csv -wa -wb > ${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.depth.filter2.exon_intersect.bed
}

calculate_uniformity()
{
log_info "#####generating uniformity calculation sample- $SAMPLE "
/opt/bedtools2/bin/bedtools coverage -a /home/environments/ngs_${ENVIRONMENT}/assayCommonFiles/tmbAssay/cre_design.bed -b ${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.filter.bed -mean > ${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.CREcoverage.mean.bed
# /opt/python3/bin/python3 /home/hhadmin/exome_pipeline/01_bamQC/06_uniformityPlots.py  ${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.CREcoverage.mean.bed

/opt/samtools19/bin/samtools view ${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.sort.rmdups.bam | awk '{ n=length($10); print gsub(/[AaTt]/,"",$10)/n;}' > ${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.ATCount.txt
}

# ##############################################################################
# main
# ##############################################################################
main()
{
	parse_options $*

	if [ $? -eq 0 ] ; then
			log_error "Import flag non-zero. Aborting $0"
			exit 1
	fi

	DIR_SCRIPT="/home/pipelines/ngs_${ENVIRONMENT}/"
	REF_GENOME="/home/doc/ref/ref_genome/ucsc.hg19.fasta"
	MAP_QUALITY="20"

  load_modules

	OUTPUT_DIR_SAMPLE="${OUTPUT_DIR}${SAMPLE}/"
	create_dir $OUTPUT_DIR_SAMPLE

	LOG_FILE="${OUTPUT_DIR}${SAMPLE}.log"


  exec >  >(tee -a ${LOG_FILE})
  exec 2> >(tee -a ${LOG_FILE} >&2)

  log_info "Starting alignment"

	OUTPUT_DIR_SAMPLE_ALIGNMENT="${OUTPUT_DIR_SAMPLE}Alignment/"
	create_dir $OUTPUT_DIR_SAMPLE_ALIGNMENT


	log_info "
	Starting - runAlignment_Interface.sh
	SAMPLE - $SAMPLE
	FASTQ_DIR - $FASTQ_DIR
	OUTPUT_DIR - $OUTPUT_DIR
	OUTPUT_DIR_SAMPLE - $OUTPUT_DIR_SAMPLE
	OUTPUT_DIR_SAMPLE_ALIGNMENT - $OUTPUT_DIR_SAMPLE_ALIGNMENT
	LOG_FILE - $LOG_FILE "

	run_trimmomatic

  run_alignment

	sort_bam

	generate_alignStat

  remove_duplicates

	generate_bamIndex

	calculate_targetCoverage

	calculate_breadthCoverage

}

# ##############################################################################
# run main
# ##############################################################################
main $*
