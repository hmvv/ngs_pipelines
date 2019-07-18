#!/bin/bash
#===============================================================================
#
# FILE: illuminaPipelineInterface.sh
#
#DESCRIPTION: This script is run to generate variant file and amplicon files
#             for heme assay.
# REQUIREMENTS: runPipelines.sh
# COMPANY:Houston Methodist Hospital, Molecular Diagnostic Laboratory
#===============================================================================

#PBS -l mem=16gb,nodes=1:ppn=4
#PBS -l walltime=10:00:00
#PBS -q default
#PBS Error_Path=${PBS_JOBNAME}.err
#PBS Output_Path=${PBS_JOBNAME}.out
#PBS -k eo

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

	echo $dateString > ${HOME}runDate.txt

}

generate_variantFile()
{
  ### alignment ########
  update_status "$QUEUEID" "Alignment" "$DB" "$USER"  "$PASSWORD"
  heme_bwaAlign $FASTQ_R1 $FASTQ_R2 ${HOME}variantCaller


  if [ ! -f "${HOME}variantCaller/${SAMPLENAME}.sort.bam" ] ; then

    log_error "ERROR:Alignment output file not found"
    update_status "$QUEUEID" "ERROR:Alignment" "$DB" "$USER"  "$PASSWORD"
    exit

  else

    log_info "Completed Alignment"

  fi


  ### variant caller ########
  update_status "$QUEUEID" "VariantCaller" "$DB" "$USER"  "$PASSWORD"
  heme_varscan ${HOME}variantCaller/${SAMPLENAME}.sort.bam  ${HOME}variantCaller

  if [ ! -f "${HOME}variantCaller/${SAMPLENAME}.snp.txt" ] ; then

    log_error "ERROR:VariantCaller output file not found"
    update_status "$QUEUEID" "ERROR:VariantCaller" "$DB" "$USER"  "$PASSWORD"
    exit

  else

    log_info "Completed VariantCaller"

  fi


  bash ${HOME_SHELLDIR}hemeParseVarScan.sh \
          -s ${HOME}variantCaller/${SAMPLENAME}.snp.txt \
          -i ${HOME}variantCaller/${SAMPLENAME}.indel.txt \
          -o ${HOME}variantCaller/  \
          -e $ENVIRONMENT

   VARIANT_VCF=$(ls ${HOME}variantCaller/${SAMPLENAME}.comb.vcf)
}

generate_ampliconFile(){

  HEME_EXCLUDED_DESIGN="/home/doc/ref/Heme/trusight-myeloid-amplicon-track.excluded.bed"

  /opt/samtools-1.4/samtools-1.4/samtools bedcov  $HEME_EXCLUDED_DESIGN  ${HOME}variantCaller/${SAMPLENAME}.sort.bam | awk ' {print $4,"\t",int($13/($8-$7))} ' > ${HOME}variantAnalysis/${SAMPLENAME}.samtools.coverageDepth

  AMPLICON_FILE=$(ls ${HOME}variantAnalysis/${SAMPLENAME}.samtools.coverageDepth)

  if [ ! -f "${HOME}variantAnalysis/${SAMPLENAME}.samtools.coverageDepth" ] ; then

    log_error "ERROR:AmpliconFile output file not found"
    update_status "$QUEUEID" "ERROR:AmpliconFile" "$DB" "$USER"  "$PASSWORD"
    exit

  else
    log_info "Completed AmpliconFile"

  fi

}

prep_heme()
{

  HEME_EXCLUDED_DESIGN="/home/doc/ref/Heme/trusight-myeloid-amplicon-track.excluded.bed"
  #bedtools -u flag to write original A entry once if any overlaps found in B
  /opt/software/bedtools-2.17.0/bin/bedtools intersect -u -a $VARIANT_VCF -b $HEME_EXCLUDED_DESIGN > ${HOME}variantAnalysis/${SAMPLENAME}.filter.vcf


  HEME_EXCLUDED_AMPLICON="/home/doc/ref/Heme/excludedAmplicons.txt"
  # -v means "invert the match" in grep, in other words, return all non matching lines.
  grep -v -f $HEME_EXCLUDED_AMPLICON $AMPLICON_FILE > ${HOME}variantAnalysis/${SAMPLENAME}.amplicon.filter.txt

}

run_illuminaPipeline()
{
  bash ${HOME_SHELLDIR}illuminaPipeline.sh -d $HOME \
        -s $SAMPLENAME -e $ENVIRONMENT -q $QUEUEID -u $USER -p $PASSWORD
}

# ##############################################################################
# main
# ##############################################################################
main()
{
    parse_options $*

    if [ $? -eq 0 ]
    then
			  echo "Error in previous step. Aborting $0"
        exit 0
    fi

		############################################################################
		# initialize variables
		############################################################################

    load_modules

    RUNFOLDER=$(ls -d /home/$INSTRUMENT/*_"$RUNID"_*)
    RUNNAME=${RUNFOLDER##/home/$INSTRUMENT/}

    HOME="/home/environments/ngs_${ENVIRONMENT}/${INSTRUMENT}Analysis/${RUNNAME}/${SAMPLENAME}/"
		HOME_SHELLDIR="/home/pipelines/ngs_${ENVIRONMENT}/shell/"
    DB="ngs_${ENVIRONMENT}"

    update_status "$QUEUEID" "started" "$DB" "$USER"  "$PASSWORD"

    exec >  >(tee -a ${HOME}process.log)
    exec 2> >(tee -a ${HOME}process.log >&2)

    ##Aligning fastq files
    FASTQ_R1=$RUNFOLDER/out1/"$SAMPLENAME"*_R1_001.fastq.gz
    FASTQ_R2=${FASTQ_R1/_R1_/_R2_}      #replace "R1" with "R2"

    VARIANT_VCF=""
    AMPLICON_FILE=""

		create_dir ${HOME}variantCaller
		create_dir ${HOME}variantAnalysis

    show_pbsinfo

    log_info " Running illumina Pipeline Interface BY qsub for :
		assay : $ASSAY
		instrument : $INSTRUMENT
		runID : $RUNID
		sampleName : $SAMPLENAME
		environment : $ENVIRONMENT
		queueID : $QUEUEID
    fastq1 : $FASTQ_R1
    fastq2 : $FASTQ_R2 "

    create_rundate

    if [ $ASSAY == "heme" ] ; then

       generate_variantFile

       generate_ampliconFile

       log_info " variant file is - $VARIANT_VCF"
       log_info " amplicon file is - $AMPLICON_FILE"

       prep_heme

       run_illuminaPipeline


    fi


}

# ##############################################################################
# run main
# ##############################################################################
main $*
