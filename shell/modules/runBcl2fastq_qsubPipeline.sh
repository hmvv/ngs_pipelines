#!/bin/bash
#PBS -l mem=32gb,nodes=1
#PBS -l walltime=5:00:00
#PBS -q default
# ##############################################################################
# # functions
# ##############################################################################

display_usuage()
{
cat <<EOF >> /dev/stderr

 USAGE: $0

 OPTIONS:
 u - user
 p - password
 d - database
 b - database host
 r - runID
EOF
}

parse_options()
{
    IMPORT=0

		while getopts "hz:u:p:d:b:r:i:e:" opt ; do
				case $opt in
					h)
						 display_usuage
						 exit 1
						 ;;
					z)
						IMPORT=$OPTARG
						;;
					u)
					  USER=$OPTARG
						;;
					p)
						PASSWORD=$OPTARG
						;;
				  d)
				    DB=$OPTARG
					  ;;
				  b)
						DB_HOST=$OPTARG
						;;
					r)
						RUNID=$OPTARG
						;;
          i)
            INSTRUMENT=$OPTARG
            ;;
          e)
  					ENVIRONMENT=$OPTARG
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

  DIR=$(ls -d /home/$INSTRUMENT/*_"$RUNID"_*)

  /usr/local/bin/bcl2fastq --no-lane-splitting --runfolder-dir $DIR --output-dir "${DIR}/out1" -d 10 -p 10
  wait $!

  updatestatement="update pipelineStatusBcl2Fastq join instruments on pipelineStatusBcl2Fastq.instrumentID = instruments.instrumentID set status=1 where pipelineStatusBcl2Fastq.runID='$RUNID' and instruments.instrumentName='$INSTRUMENT'"
  mysql --host=$DB_HOST --user="$USER" --password="$PASSWORD" --database="$DB" --execute="$updatestatement"

}

main $*
