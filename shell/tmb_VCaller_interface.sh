#!/bin/bash

# ##############################################################################
# # functions
# ##############################################################################

display_usage()
{
cat <<EOF >> /dev/stderr

 USAGE: $0

 OPTIONS:
 	-n Normal Sample
 	-t Tumor Sample
	-v Variant Callers
	-o Output Directory

EOF
}

parse_options()
{
    IMPORT=0

		while getopts "hz:n:t:v:o:e:i:q:u:p:" opt; do

				case $opt in
					h)
					display_usage
					exit 1
					;;
	        z)
					IMPORT=$OPTARG
					;;
					n)
					NORMAL_BAM=$OPTARG
					;;
					t)
					TUMOR_BAM=$OPTARG
					;;
					v)
					VARIANT_CALLERS=$OPTARG
					;;
					o)
					OUT_DIR=$OPTARG
					;;
          e)
          ENVIRONMENT=$OPTARG
          ;;
          i)
          SAMPLEID=$OPTARG
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
	source /home/pipelines/ngs_${ENVIRONMENT}/shell/modules/ngs_varscan.sh
  source /home/pipelines/ngs_${ENVIRONMENT}/shell/modules/ngs_vep.sh
}

run_UNpairedVC()
{
for v_caller in $VCS;do

	create_dir ${SAMPLE_DIR}/${v_caller}

	if [[ "$v_caller" == "varscan" ]];then

		log_info "Starting : " $v_caller

		germline_varscan $NORMAL_BAM   ${SAMPLE_DIR}/${v_caller}/

    bash /home/pipelines/ngs_test/shell/heme_parseVarScan.sh \
            -s ${SAMPLE_DIR}/${v_caller}/${SAMPLE}.snp.txt \
            -i ${SAMPLE_DIR}/${v_caller}/${SAMPLE}.indel.txt \
            -o ${SAMPLE_DIR}/${v_caller}/  \
            -e $ENVIRONMENT

	fi

done
}

run_pairedVC()
{

	DEPTH="10"
	NALF="10"
	TALF="10"

	for v_caller in $VCS;do

		create_dir ${SAMPLE_DIR}/${v_caller}

		log_info "Starting : " $v_caller

    ### run vc
    bash ${DIR_SCRIPT}shell/tmb_VCaller_${v_caller}.sh $SAMPLE $REF_GENOME_1  $NORMAL_BAM  $TUMOR_BAM  ${SAMPLE_DIR}/${v_caller}/  $ENVIRONMENT  $DEPTH  $NALF $TALF

    ## parse vc output
    /opt/python3/bin/python3 ${DIR_SCRIPT}python/tmb_parse_${v_caller}.py  "$SAMPLE"  "${SAMPLE_DIR}/${v_caller}/"  "$ENV" "$DEPTH" "$NALF" "$TALF"

	done

  # compare all three variants files
	/opt/python3/bin/python3 ${DIR_SCRIPT}python/tmb_compare_samples_3Venn.py   \
	${SAMPLE_DIR}/varscan/${SAMPLE}.varscan.${DEPTH}_${NALF}_${TALF}  \
	${SAMPLE_DIR}/strelka/${SAMPLE}.strelka.${DEPTH}_${NALF}_${TALF} \
	${SAMPLE_DIR}/mutect/${SAMPLE}.mutect.${DEPTH}_${NALF}_${TALF} \
	${SAMPLE_DIR}/ "${DEPTH}_${NALF}_${TALF}"


	# # #### combine vcfs from EACH variant caller
	/opt/python3/bin/python3 ${DIR_SCRIPT}python/tmb_combineVCFs.py  "$SAMPLE"  "$SAMPLE_DIR"  "$ENVIRONMENT"  "${DEPTH}_${NALF}_${TALF}"


  ###### run VEP
  update_status "$QUEUEID" "RunningVEP" "$DB" "$USER"  "$PASSWORD"

	tail -n +2 "${SAMPLE_DIR}/${SAMPLE}.variantcallers.combine.${DEPTH}_${NALF}_${TALF}" > ${SAMPLE_DIR}/${SAMPLE}.variantcallers.combinev2.${DEPTH}_${NALF}_${TALF}

  log_info  "running VEP"

  vep_94_tmb ${SAMPLE_DIR}/${SAMPLE}.variantcallers.combinev2.${DEPTH}_${NALF}_${TALF}   ${SAMPLE_DIR}/${SAMPLE}.variantcallers.combinev2.${DEPTH}_${NALF}_${TALF}.vep

  update_status "$QUEUEID" "CompletedVEP" "$DB" "$USER"  "$PASSWORD"

	##### FILTER vcfs from all variant callers
	/opt/python3/bin/python3 ${DIR_SCRIPT}python/tmb_parse_vep.py  "$SAMPLEID" "$SAMPLE"  "$SAMPLE_DIR"  "$ENVIRONMENT"  "${DEPTH}_${NALF}_${TALF}"

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
  DB="ngs_${ENVIRONMENT}"

	REF_GENOME_1="/home/doc/ref/ref_genome/ucsc.hg19.fasta"
	REF_GENOME_2="/home/doc/ref/ref_genome/Homo_sapiens.GRCh37.73.dna_sm.primary_assembly.fa"

  load_modules

	if [[ $TUMOR_BAM == "NONE" ]];then
		SAMPLE=$(echo "${NORMAL_BAM##*/}"| tr "." "\n" | head -1 )
	else
		SAMPLE=$(echo "${TUMOR_BAM##*/}"| tr "." "\n" | head -1 )
		SAMPLE+="_"$(echo "${NORMAL_BAM##*/}"| tr "." "\n" | head -1 )
	fi


	SAMPLE_DIR=${OUT_DIR}${SAMPLE}
	create_dir $SAMPLE_DIR

  LOG_FILE="${SAMPLE_DIR}${SAMPLE}.log"

  exec >  >(tee -a ${LOG_FILE})
  exec 2> >(tee -a ${LOG_FILE} >&2)

  log_info "Starting exome variant Callers.."

  #get all variant callers in list
	VCS=$(echo $VARIANT_CALLERS | tr "-" "\n")

	if [[ $TUMOR_BAM == "NONE" ]];then
		run_UNpairedVC
	else
		run_pairedVC
	fi

}

# ##############################################################################
# run main
# ##############################################################################
main $*
