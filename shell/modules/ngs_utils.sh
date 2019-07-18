#!/bin/bash

log_info()
{
   log " [INFO] : $1"
}

log_error()
{
   log " [ERROR] : $1"
}

log()
{
 MESSAGE=$1
 TIMESTAMP=`date "+%Y-%m-%d %H:%M:%S"`
 SCRIPT=$( basename $0 )
 echo " [ $TIMESTAMP ] [ $SCRIPT ] $MESSAGE "
}

create_file()
{
  path=$1
  file_name=$2
	touch  ${path}${file_name}
	chmod 775  ${path}${file_name}
}

create_dir()
{
	new_dir=$1
	if [ ! -d $new_dir ] ; then
		mkdir $new_dir
	fi
	chmod 775 $new_dir
}

update_status()
{
user=$4
password=$5
database=$3
insertstatement="INSERT INTO pipelineStatus (queueID, plStatus, timeUpdated) VALUES ('$1','$2',now());"
mysql --host="hhplabngsp01" --user="$user" --password="$password" --database="$database" --execute="$insertstatement"
}

show_pbsinfo()
{
  echo ------------------------------------------------------
  echo -n 'Job is running on node '; cat $PBS_NODEFILE
  echo ------------------------------------------------------
  echo PBS: qsub is running on $PBS_O_HOST
  echo PBS: originating queue is $PBS_O_QUEUE
  echo PBS: executing queue is $PBS_QUEUE
  echo PBS: working directory is $PBS_O_WORKDIR
  echo PBS: execution mode is $PBS_ENVIRONMENT
  echo PBS: job identifier is $PBS_JOBID
  echo PBS: job name is $PBS_JOBNAME
  echo PBS: node file is $PBS_NODEFILE
  echo PBS: current home directory is $PBS_O_HOME
  echo PBS: PATH = $PBS_O_PATH
  echo ------------------------------------------------------

}
