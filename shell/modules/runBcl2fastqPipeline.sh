#!/bin/bash

load_modules()
{
      source /home/pipelines/ngs_${ENVIRONMENT}/shell/modules/ngs_utils.sh
}


run_bcl2fastq()
{
USER=$1
PASSWORD=$2
DB=$3
DB_HOST=$4
RUNID=$5
INSTRUMENT=$6
ENVIRONMENT=$7

load_modules

log_info "submiting bcl2fastq from runPipelines, runID is $RUNID"

/opt/torque/bin/qsub -k eo -F "-u$USER -p$PASSWORD -d$DB  -b$DB_HOST  -r$RUNID -i$INSTRUMENT" /home/pipelines/ngs_${ENVIRONMENT}/shell/modules/runBcl2fastq_qsubPipeline.sh

}
