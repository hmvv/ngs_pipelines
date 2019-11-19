#!/bin/bash
#===============================================================================
#
# FILE: heme_interface.sh
#
#DESCRIPTION: This script is run to generate required variant files
#             for heme .
# REQUIREMENTS: illuminaPipelineInterface.sh
# COMPANY:Houston Methodist Hospital, Molecular Diagnostic Laboratory
#===============================================================================

# ##############################################################################
# # functions
# ##############################################################################
display_usage()
{
cat <<EOF >> /dev/stderr

 USAGE: $0

 OPTIONS:
 r - runID
 s - sampleName
 a - assay
 i - instrument
 e - environment
 q - queueID
 u - user
 p - password

EOF
}

parse_options()
{
    IMPORT=0

		while getopts "hz:r:s:a:i:e:q:u:p:" opt ; do
				case $opt in
					h)
						 display_usage
						 exit 1
						 ;;
					z)
						IMPORT=$OPTARG
						;;
					r)
						RUNID=$OPTARG
						;;
					s)
					  SAMPLENAME=$OPTARG
						;;
					a)
						ASSAY=$OPTARG
						;;
					i)
				    INSTRUMENT=$OPTARG
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
      source /home/pipelines/ngs_${ENVIRONMENT}/shell/modules/ngs_varscan.sh
}

create_rundate()
{
	##get runDate information##
	dateString=${RUNNAME%%_*}
  echo $dateString
	year=20${dateString:0:2}
	month=${dateString:2:2}
	day=${dateString:4:2}
	dateString=$year-$month-$day

	echo $dateString > ${HOME_ANALYSIS}runDate.txt

}

heme_run_alignment()
{

  update_status "$QUEUEID" "Trimming" "$DB" "$USER"  "$PASSWORD"

  log_info "Running Trimmomatic: Removing sequences < Q20 sample- $SAMPLENAME"
	trimmomatic="java -jar /opt/trimmomatic/Trimmomatic-0.33/trimmomatic-0.33.jar PE -phred33 -threads 8 \
	$FASTQ_R1 \
	$FASTQ_R2 \
	${HOME_ANALYSIS}variantCaller/${SAMPLENAME}_filt_paired_R1_001.fastq.gz \
	${HOME_ANALYSIS}variantCaller/${SAMPLENAME}_filt_unpaired_R1_001.fastq.gz \
	${HOME_ANALYSIS}variantCaller/${SAMPLENAME}_filt_paired_R2_001.fastq.gz \
	${HOME_ANALYSIS}variantCaller/${SAMPLENAME}_filt_unpaired_R2_001.fastq.gz \
	TRAILING:20 \
	AVGQUAL:20 \
	SLIDINGWINDOW:10:20 \
	MINLEN:30 "
	($trimmomatic) 2>&1 | tee ${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.trimmomatic.summary.txt


  ### alignment ########
  update_status "$QUEUEID" "Alignment" "$DB" "$USER"  "$PASSWORD"
  log_info "Running bwa mem aligner: $SAMPLENAME"
  tmb_bwaAlign $SAMPLENAME  \
	$REF_GENOME  \
	$MAP_QUALITY \
	${HOME_ANALYSIS}variantCaller/${SAMPLENAME}_filt_paired_R1_001.fastq.gz \
	${HOME_ANALYSIS}variantCaller/${SAMPLENAME}_filt_paired_R2_001.fastq.gz \
	${HOME_ANALYSIS}variantCaller/  \
	${HOME_ANALYSIS}process.log

  log_info "Generating sort bam by coordinate sample- $SAMPLENAME"
  java -jar /opt/picard2/picard.jar SortSam \
            I=${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.bam  \
            O=${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.sort.bam  \
            SORT_ORDER=coordinate

  log_info "Generating alignment stat sample- $SAMPLENAME"
  java -jar /opt/picard2/picard.jar CollectAlignmentSummaryMetrics \
            R=$REF_GENOME \
            I=${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.sort.bam \
            O=${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.sort.bam.alignmentMetrics.txt

  log_info "Generating bam index sample- $SAMPLENAME"
  java -jar /opt/picard2/picard.jar BuildBamIndex \
            I=${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.sort.bam

  log_info "Generating CalculateHsMetrics sample- $SAMPLENAME "
  java -jar /opt/picard/picard-tools-1.134/picard.jar CalculateHsMetrics \
            I=${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.sort.bam  \
            O=${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.output_hs_metrics.txt \
            R=$REF_GENOME \
            BAIT_INTERVALS= /home/environments/ngs_${ENVIRONMENT}/assayCommonFiles/hemeAssay/myeloid_design.interval_list \
            TARGET_INTERVALS= /home/environments/ngs_${ENVIRONMENT}/assayCommonFiles/hemeAssay/myeloid_design.interval_list


  if [ ! -f "${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.sort.bam" ] ; then

    log_error "ERROR:Alignment output file not found"
    update_status "$QUEUEID" "ERROR:Alignment" "$DB" "$USER"  "$PASSWORD"
    exit

  else

    log_info "Completed Alignment"

  fi
}


heme_run_variantCaller()
{

  log_info "Running variant calling for - $SAMPLENAME"

	VCS=( "varscan" "freebayes" "mutect")

  update_status "$QUEUEID" "VariantCaller" "$DB" "$USER"  "$PASSWORD"

	for v_caller in "${VCS[@]}";do

    if [[ "$v_caller" == "varscan" ]];then

      log_info "Starting :  $v_caller "

      heme_varscan ${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.sort.bam   ${HOME_ANALYSIS}variantCaller

      bash ${HOME_SHELLDIR}heme_parseVarScan.sh \
              -s ${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.snp.varscan.output \
              -i ${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.indel.varscan.output \
              -o ${HOME_ANALYSIS}variantCaller/  \
              -f ${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.varscan.vcf  \
              -e $ENVIRONMENT

    elif [[ "$v_caller" == "freebayes" ]];then

      log_info "Starting :  $v_caller "

      /opt/freebayes2/freebayes/bin/freebayes -f $REF_GENOME   ${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.sort.bam  > ${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.freebayes.output

      java -jar /opt/GATK4/GenomeAnalysisTK.jar \
           -R $REF_GENOME \
           -T VariantsToTable \
           --showFiltered \
           -V ${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.freebayes.output \
           -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER \
           -GF DP -GF RO  -GF AO  \
           -o ${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.freebayes.vcf


    elif [[ "$v_caller" == "mutect" ]];then

      log_info "Starting :  $v_caller "

      java  -jar /opt/GATK4/GenomeAnalysisTK.jar \
            -T MuTect2 \
            -R $REF_GENOME \
            -I:tumor ${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.sort.bam   \
            -o ${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.mutect.output


      java  -jar /opt/GATK4/GenomeAnalysisTK.jar \
            -R $REF_GENOME \
            -T VariantsToTable \
            --showFiltered \
            -V ${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.mutect.output \
            -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER \
            -GF GT -GF AD -GF AF  \
            -o ${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.mutect.vcf

    fi

	done
}

heme_selects_variants()
{
  /opt/python3/bin/python3   ${HOME_PYTHONDIR}heme_selectHQVariants.py \
        ${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.varscan.vcf \
        ${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.mutect.vcf \
        ${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.freebayes.vcf \
        ${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.combine.vcf \

}

# ##############################################################################
# main
# ##############################################################################
main()
{
    parse_options $*

    if [ $? -eq 0 ];then
			  echo "Error in previous step. Aborting $0"
        exit 0
    fi

		############################################################################
    # initialize variables
		############################################################################

    load_modules

    RUNFOLDER=$(ls -d /home/$INSTRUMENT/*_"$RUNID"_*)
    RUNNAME=${RUNFOLDER##/home/$INSTRUMENT/}

    HOME_CODE="/home/pipelines/ngs_${ENVIRONMENT}/"
    HOME_RUNDIR="${HOME_CODE}run_files/"
    HOME_SHELLDIR="${HOME_CODE}shell/"
    HOME_PYTHONDIR="${HOME_CODE}python/"
    DB="ngs_${ENVIRONMENT}"
    DB_HOST="hhplabngsp01"

    HOME_ANALYSIS="/home/environments/ngs_${ENVIRONMENT}/${INSTRUMENT}Analysis/${ASSAY}Assay/${RUNNAME}/${SAMPLENAME}/"

    FASTQ_R1=$RUNFOLDER/out1/"$SAMPLENAME"*_R1_001.fastq.gz
    FASTQ_R2=${FASTQ_R1/_R1_/_R2_}      #replace "R1" with "R2"
    REF_GENOME="/home/doc/ref/ref_genome/ucsc.hg19.fasta"
    MAP_QUALITY="20"

    log_info " Running heme interface for :
		assay : $ASSAY
		instrument : $INSTRUMENT
		runID : $RUNID
		sampleName : $SAMPLENAME
		environment : $ENVIRONMENT
		queueID : $QUEUEID
    fastq1 : $FASTQ_R1
    fastq2 : $FASTQ_R2
    ref: $REF_GENOME
    mapQ: $MAP_QUALITY "

    heme_run_alignment

    heme_run_variantCaller

    heme_selects_variants
}

# ##############################################################################
# run main
# ##############################################################################
main $*
