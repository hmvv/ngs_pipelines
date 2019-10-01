#!/bin/bash
#===============================================================================
#
# FILE: runPipelines.sh
#
#DESCRIPTION: This script is run by cronjob.It checks the queue table in NGS
#             database and submits a sample analysis qsub based
#             on instrument and assay.
# REQUIREMENTS: cron
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
 	-e ENVIRONMENT
 	-u USER
	-p PASSWORD

EOF
}

parse_options()
{
    IMPORT=0

		while getopts "hz:e:u:p:" opt; do

				case $opt in
					h)
					display_usage
					exit 1
					;;
	        z)
					IMPORT=$OPTARG
					;;
					e)
					ENVIRONMENT=$OPTARG
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
      source /home/pipelines/ngs_${ENVIRONMENT}/shell/modules/runBcl2fastqPipeline.sh
}

update_status_host()
{
user=$4
password=$5
database=$3
insertstatement="INSERT INTO pipelineStatus (queueID, plStatus, timeUpdated) VALUES ('$1','$2',now());"
mysql  --user="$user" --password="$password" --database="$database" --execute="$insertstatement"
}

# ##############################################################################
# workflow:check_db_queue
# ##############################################################################
check_db_queue()
{
	default_Status=0
	status_query_statement="select pipelineQueue.queueID, samples.runID, samples.sampleName, \
	samples.coverageID, samples.callerID, assays.assayName, instruments.instrumentName, \
	pipelineQueue.status from pipelineQueue \
	join samples on samples.sampleID=pipelineQueue.sampleID \
	join assays on assays.assayID = samples.assayID \
	join instruments on instruments.instrumentID = samples.instrumentID \
	where pipelineQueue.status='$default_Status';"

  while read -r queueID runID  sampleName coverageID callerID assay instrument status; do

			log_info "Processing the following job queued in the NGS database:
			ENVIRONMENT - $ENVIRONMENT
			INSTRUMENT - $instrument
			ASSAY - $assay
			QUEUEID - $queueID
			QUEUE_STATUS - $status
			RUNID - $runID
			SAMPLE - $sampleName
			COVERAGE - $coverageID
			CALLER - $callerID"

	   ##update sampleAnalysisQueue table and set status of this queue to 1 i.e started processing
	   updatestatement="update pipelineQueue SET status=1 where queueID = $queueID;"
		 mysql  --user="$USER" --password="$PASSWORD" --database="$DB" --execute="$updatestatement"

     if [ "$instrument" == "proton" ] ; then

       if [ ! -f ${HOME_RUN}${CURRENTDT}_ProtonQueued.samples ]; then

         create_file "$HOME_RUN"  "${CURRENTDT}_ProtonQueued.samples"
         echo "##queueID;runID;sampleName;coverageID;callerID;assay;instrument;ENVIRONMENT##RUNDATE-$CURRENTDT"  >> ${HOME_RUN}${CURRENTDT}_ProtonQueued.samples

       fi

       echo "$queueID;$runID;$sampleName;$coverageID;$callerID;$assay;$instrument;$ENVIRONMENT" >> ${HOME_RUN}${CURRENTDT}_ProtonQueued.samples

     elif  [ $(echo ${instrument:0:7}) == "nextseq" ] ; then

       if [ ! -f ${HOME_RUN}${CURRENTDT}_IlluminaQueued.samples ]; then

         create_file "$HOME_RUN"  "${CURRENTDT}_IlluminaQueued.samples"
         echo "##queueID;runID;sampleName;coverageID;callerID;assay;instrument;ENVIRONMENT##RUNDATE-$CURRENTDT"  >> ${HOME_RUN}${CURRENTDT}_IlluminaQueued.samples

       fi

       echo "$queueID;$runID;$sampleName;$coverageID;$callerID;$assay;$instrument;$ENVIRONMENT" >> ${HOME_RUN}${CURRENTDT}_IlluminaQueued.samples

     fi


  done < <(mysql --user="$USER" --password="$PASSWORD" --database="$DB" --execute="$status_query_statement" -N)

}

# ##############################################################################
# workflow:submit job
# ##############################################################################
submit_jobs_proton()
{

  if [ ! -f  ${HOME_RUN}${CURRENTDT}_ProtonQueued.samples  ] ; then
    log_info "No samples queued from proton."
    return 1
  fi

  log_info "Samples queued from proton."

  tail -n +2 ${HOME_RUN}${CURRENTDT}_ProtonQueued.samples | while IFS=';' read -ra line; do

    queueID="${line[0]}"
  	runID="${line[1]}"
  	sampleName="${line[2]}"
  	coverageID="${line[3]}"
  	callerID="${line[4]}"
    assay="${line[5]}"
    instrument="${line[6]}"


    log_info "current sample - $queueID ; $runID ; $instrument ; $assay ; $sampleName"

    home_analysis_instrument="${HOME_ANALYSIS}${instrument}Analysis/"

    runFolder=$(ls -d /home/${instrument}/*$runID)
    runName=${runFolder##*/}

    create_dir ${home_analysis_instrument}$runName
    create_dir ${home_analysis_instrument}${runName}/$sampleName

    working_dir="${home_analysis_instrument}${runName}/${sampleName}/"

    /opt/torque/bin/qsub -d ${working_dir}  \
           -F "-r$runID -s$sampleName -c$coverageID -v$callerID -a$assay -i$instrument -e$ENVIRONMENT -q$queueID -u$USER -p$PASSWORD" \
           ${HOME_SHELL}protonPipelineInterface.sh

  done
}

submit_sample_illumina()
{

  queueID=$1
  runID=$2
  instrument=$3
  assay=$4
  sampleName=$5

  runFolder=$(ls -d /home/$instrument/*_"$runID"_*)
  runName=${runFolder##/home/$instrument/}
  home_analysis_instrument="${HOME_ANALYSIS}${instrument}Analysis/${assay}Assay/"

  create_dir "${home_analysis_instrument}$runName"
  create_dir "${home_analysis_instrument}${runName}/$sampleName"

  working_dir="${home_analysis_instrument}${runName}/${sampleName}/"

  log_info "Submitting illuminaPipelineInterface ; current sample - $queueID ; $runID ; $instrument ; $assay ; $sampleName"

  if [ $assay == "heme" ] ; then

    /opt/torque/bin/qsub -d ${working_dir} -l "walltime=10:00:00,nodes=1:ppn=4"   \
         -F "-r$runID -s$sampleName -a$assay -i$instrument -e$ENVIRONMENT -q$queueID -u$USER -p$PASSWORD" \
         ${HOME_SHELL}illuminaPipelineInterface.sh

  elif [ $assay == "tmb" ] ; then

    /opt/torque/bin/qsub -d ${working_dir} -l  "walltime=50:00:00,nodes=1:ppn=12"  \
         -F "-r$runID -s$sampleName -a$assay -i$instrument -e$ENVIRONMENT -q$queueID -u$USER -p$PASSWORD" \
         ${HOME_SHELL}illuminaPipelineInterface.sh

  fi

    # /opt/torque/bin/qsub -d ${working_dir}  \
    #      -F "-r$runID -s$sampleName -a$assay -i$instrument -e$ENVIRONMENT -q$queueID -u$USER -p$PASSWORD" \
    #      ${HOME_SHELL}illuminaPipelineInterface.sh

}

process_sample_llumina()
{
  line=$1
  queueID="${line[0]}"
  runID="${line[1]}"
  sampleName="${line[2]}"
  coverageID="${line[3]}"
  callerID="${line[4]}"
  assay="${line[5]}"
  instrument="${line[6]}"

  log_info "current sample - $queueID ; $runID ; $instrument ; $assay ; $sampleName"

  stausCheck="select pipelineStatusBcl2Fastq.status from pipelineStatusBcl2Fastq join instruments on pipelineStatusBcl2Fastq.instrumentID = instruments.instrumentID where pipelineStatusBcl2Fastq.runID='$runID' and instruments.instrumentName='$instrument'"
  fastqStatus=$(mysql --user="$USER" --password="$PASSWORD" --database="$DB" -se "$stausCheck")

  log_info "current fastqStatus - $fastqStatus, $queueID ; $runID ; $instrument ; $assay ; $sampleName"

  if [ "$fastqStatus" == 0 ] ; then

    lockdir=${HOME_RUN}nextseq_${runID}.lock

    ## mkdir is atomic, file or file variable is not
    if mkdir "$lockdir" ; then

      log_info "Acquired lock for run-$runID , sample- $sampleName , lockfile- $lockdir"
      log_info "bcl2fastq_starting_now run-$runID , sample- $sampleName "
      update_status_host "$queueID" "bcl2fastq_starting_now" "$DB" "$USER"  "$PASSWORD"

      run_bcl2fastq $USER  $PASSWORD  $DB  $DB_HOST  $runID $instrument $ENVIRONMENT

      ### sleep for an hour and keep checking database if bcl2fastq for this run is completed
      fastqStatus_r2="0"
      while [ "$fastqStatus_r2" == "0" ]; do

          update_status_host "$queueID" "bcl2fastq_running..." "$DB" "$USER"  "$PASSWORD"

          sleep 1h

          fastqStatus_r2=$(mysql --user="$USER" --password="$PASSWORD" --database="$DB" -se "$stausCheck")

          log_info "checking bcl2fastq status -- $fastqStatus_r2, run-$runID , sample- $sampleName"


      done

      log_info  "bcl2fastq_completed_now,  run-$runID , sample- $sampleName"
      update_status_host "$queueID" "bcl2fastq_completed_now" "$DB" "$USER"  "$PASSWORD"

      submit_sample_illumina $queueID  $runID $instrument $assay $sampleName

      rm -rf "$lockdir"

    else
      updatestatement="UPDATE pipelineQueue SET pipelineQueue.status=0 WHERE queueID = $queueID;"
      mysql --user="$USER" --password="$PASSWORD" --database="$DB" --execute="$updatestatement"

      log_info  "bcl2fastq_wait,  run-$runID , sample- $sampleName"
      update_status_host "$queueID" "bcl2fastq_wait" "$DB" "$USER"  "$PASSWORD"

    fi

  else

    log_info "bcl2fastq_completed_past,  run-$runID , sample- $sampleName"
    update_status_host "$queueID" "bcl2fastq_completed_past" "$DB" "$USER"  "$PASSWORD"
    submit_sample_illumina $queueID  $runID $instrument $assay $sampleName
  fi

}

submit_jobs_illumina()
{

  if [ ! -f  ${HOME_RUN}${CURRENTDT}_IlluminaQueued.samples  ] ; then
    log_info "No samples queued on illumina."
    return 1
  fi

  log_info "Samples queued on illumina."

  tail -n +2 ${HOME_RUN}${CURRENTDT}_IlluminaQueued.samples | while IFS=';' read -ra line; do

       process_sample_llumina $line &

  done
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

		# ##########################################################################
		# initialize Variables
		# ##########################################################################

    load_modules

		HOME="/home/pipelines/ngs_${ENVIRONMENT}/"
		HOME_RUN="${HOME}run_files/"
		HOME_SHELL="${HOME}shell/"
    HOME_ANALYSIS="/home/environments/ngs_${ENVIRONMENT}/"
		DB="ngs_${ENVIRONMENT}"
    DB_HOST="hhplabngsp01"
		CURRENTDT=`date '+%Y_%m_%d_%H_%M_%S'`


    # ##########################################################################
    # workflows
    # ##########################################################################
    check_db_queue

		submit_jobs_proton

    submit_jobs_illumina

    sleep 5s

    log_info "Completed cronjob for $0."
}

# ##############################################################################
# run main
# ##############################################################################
main $*
