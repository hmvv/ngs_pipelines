#!/bin/bash
#===============================================================================
#
# FILE: proton_annotationPipeline.sh
#
#DESCRIPTION: This script runs variant annotation for gene50 and neuro from
#             proton instrument,filters amplicon and variant files.
#             Then call script to load amplicon and variant files into the NGS database.
# REQUIREMENTS: protonPipelineInterface.sh
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
 d - sample directory under analysis environment
 s - sampleName
 c - coverageID
 v - callerID
 a - assay
 e - environment
 q - queueID
 u - user
 p - password

EOF
}

parse_options()
{
    IMPORT=0

		while getopts "hz:d:s:c:v:a:e:q:u:p:" opt ; do
				case $opt in
					h)
						 display_usage
						 exit 1
						 ;;
					z)
						IMPORT=$OPTARG
						;;
					d)
						HOME=$OPTARG
						;;
					s)
						SAMPLENAME=$OPTARG
						;;
					c)
						COVERAGEID=$OPTARG
						;;
					v)
						CALLERID=$OPTARG
						;;
          a)
            ASSAY=$OPTARG
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
      source /home/pipelines/ngs_${ENVIRONMENT}/shell/modules/ngs_vep.sh
}

split_protonVCF()
{
	python ${HOME_PYTHON}splitVcf.py \
				-I ${HOME}${CALLERID}/TSVC_variants.filter.vcf \
				-f FAO,FDP,AF \
				-o ${HOME}${CALLERID}/TSVC_variants.filter.split.vcf
}

run_vep()
{
	log_info "Running VEP"

  # vep_83 "${HOME}${CALLERID}/TSVC_variants.filter.split.vcf"  "${HOME}${CALLERID}/TSVC_variants.filter.split.vep.vcf"
  vep_94_panel "${HOME}${CALLERID}/TSVC_variants.filter.split.vcf"  "${HOME}${CALLERID}/TSVC_variants.filter.split.vep.vcf"

  if [ ! -f "${HOME}${CALLERID}/TSVC_variants.filter.split.vep.vcf" ] ; then

    log_error "ERROR:VEP output file not found"
    update_status "$QUEUEID" "ERROR:VEP" "$DB" "$USER"  "$PASSWORD"
    exit

  else

    log_info "Completed VEP"

  fi
}

parse_vep()
{
  log_info "Parsing VEP results"

	python ${HOME_PYTHON}parseVEP_v2.py \
					parseIonNewVarView \
					-I ${HOME}${CALLERID}/TSVC_variants.filter.split.vep.vcf \
					-o ${HOME}${CALLERID}/TSVC_variants.filter.split.vep.parse.txt
}

filter_vep()
{
	log_info "Filtering VEP variants : assay specific filtering"


  /opt/python3/bin/python3 ${HOME_PYTHON}filterVEP.py \
					-i ${HOME}${CALLERID}/TSVC_variants.filter.split.vep.parse.txt \
          -o ${HOME}${CALLERID}/TSVC_variants.filter.split.vep.parse.filter.txt \
          -a $ASSAY

  sample_ID=$(grep "^#CHROM" ${HOME}${CALLERID}/TSVC_variants.filter.vcf |cut -f 10)


  log_info "Filtering VEP variants (High/Moderate/>100)"

  shopt -s nocasematch
  if [[ $sample_ID =~ horizon ]]
	then
		awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 1 && $10 != "null" && $11 >=100 && $11 != "null") print}' ${HOME}${CALLERID}/TSVC_variants.filter.split.vep.parse.filter.txt \
		> ${HOME}${CALLERID}/TSVC_variants.filter.split.vep.parse.filter2.txt
	else
		awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 10 && $10 != "null" && $11 >=100 && $11 != "null") print}' ${HOME}${CALLERID}/TSVC_variants.filter.split.vep.parse.filter.txt \
		> ${HOME}${CALLERID}/TSVC_variants.filter.split.vep.parse.filter2.txt
	fi

}

filter_amplicon()
{
  log_info "Filtering <100 depth amplicons"

	awk -F'\t' -v OFS='\t' '{if($10 < 100) print $4,$10}' ${HOME}${COVERAGEID}/amplicon.filter.txt \
	> ${HOME}${COVERAGEID}/amplicon.lessThan100.txt
}

update_db()
{

  log_info "Updating Database"

  bash ${HOME_SHELL}heme_proton_common_DBUpdate.sh -d $HOME -s $SAMPLENAME -c $COVERAGEID -v $CALLERID -q $QUEUEID -e $ENVIRONMENT -u $USER -p $PASSWORD

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

		HOME_PYTHON="/home/pipelines/ngs_${ENVIRONMENT}/python/"
		HOME_SHELL="/home/pipelines/ngs_${ENVIRONMENT}/shell/"
		DB="ngs_${ENVIRONMENT}"

		log_info " Running proton pipeline for :
		directory : $HOME
		sampleName : $SAMPLENAME
		coverageID : $COVERAGEID
		callerID : $CALLERID
		environment : $ENVIRONMENT
		queueID : $QUEUEID "

    ## normalize variants generated by the ion proton instrument before annotation
		split_protonVCF

    update_status "$QUEUEID" "RunningVEP" "$DB" "$USER"  "$PASSWORD"

    run_vep

    update_status "$QUEUEID" "CompletedVEP" "$DB" "$USER"  "$PASSWORD"

		parse_vep

		filter_vep

		filter_amplicon

    update_status "$QUEUEID" "UpdatingDatabase" "$DB" "$USER"  "$PASSWORD"

		update_db

}


# ##############################################################################
# run main
# ##############################################################################
main $*
