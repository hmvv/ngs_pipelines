#!/usr/bin/env bash
export SHELL=/usr/bin/bash

if [ $# -eq 0 ]
then
	echo "Usage: createLocalBam.sh"
	echo "-f file name with parameters and coordinates"
	exit
fi

if test $# -gt 0
	then
	while getopts :f: opt
	do
	case $opt in
  f)
     fileName=$OPTARG
	   ;;
  :)
	echo "Option -$OPTARG requires an argument."
	;;
  \?)
	echo "Invalid option: -$OPTARG"
esac
done
   shift $((OPTIND-1))
fi


tempfilePath='/home/scratch/hmvv3/igv/'
IFS=';' read -ra params <<< "$(head -n 1 "$tempfilePath$fileName")"
environment="${params[0]}"
instrument="${params[1]}"
runID="${params[2]}"
assay="${params[3]}"
sampleName="${params[4]}"
callerID="${params[5]}"
coverageID="${params[6]}"



echo "
environment : $environment
instrument : $instrument
runID : $runID
assay : $assay
sampleName : $sampleName
callerID : $callerID
coverageID : $coverageID"

#remove previous bam/bai files
# rm -f $(find $tempfilePath -name "*.bam" -type f)
# rm -f $(find $tempfilePath -name "*.bai" -type f)

bamFile=""
if [ $instrument == "proton" ]
then

	bamFile=$(ls -f /home/$instrument/*$runID/plugin_out/$callerID/"$sampleName"_*.bam)

elif [ $instrument == "nextseq" ]
then

	bamFile=$(ls -f /home/environments/$environment/"$instrument"Analysis/*_"$runID"_*/"$sampleName"/variantCaller/"$sampleName".sort.bam)
fi

echo "BAM file path is - $bamFile"


COUNTER=0
CurrentCoordinate=""

{
	read
	while IFS='' read -r line; do
		coordinate="$(echo -e "${line}" | tr -d '[:space:]')"
		/opt/samtools/samtools view -b $bamFile $coordinate > "$tempfilePath$fileName$coordinate".bam
		COUNTER=$[$COUNTER +1]
		CurrentCoordinate=$coordinate

done
}< "$tempfilePath$fileName"


num=1
if [ "$COUNTER" -gt "$num" ]
then
	/opt/samtools/samtools merge "$tempfilePath$fileName".bam  $tempfilePath$fileName*.bam
else
	mv "$tempfilePath$fileName$CurrentCoordinate".bam "$tempfilePath$fileName".bam
fi

/opt/samtools/samtools index "$tempfilePath$fileName".bam
#
# rm -f $(find $tempfilePath -name "$fileName" -type f)
# rm -f $(find $tempfilePath -name  "$fileNamechr*"  -type f)
