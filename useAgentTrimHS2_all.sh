#!/bin/bash

#selfName=`basename "$0"`

#if [ $# -lt 3 ] || [ $# -gt 3 ]; then
#    echo "Usage: $selfName R1.fastq.gz R2.fastq.gz targetFolder"
#    echo "ie. $selfName HumanRef_50ng_2_S29_R1_001.fastq.gz HumanRef_50ng_2_S29_R2_001.fastq.gz _UMItrimmed"
#    echo "";
#    exit;
#fi
#echo $1
#cd $1

ls rawFASTQ | grep .fastq.gz | grep -v '_R2_' |while read id;
do
	echo $id
	aa=$(echo $id |sed  "s/_R1_001.fastq.gz//g")
	if [ ! -e "rawFASTQ/${aa}_R1_001.fastq.gz" ]; then  
    		echo "$R1 NOT found."
		exit;
	fi  
	if [ ! -e "rawFASTQ/${aa}_R2_001.fastq.gz" ]; then  
    	echo "$R2 NOT found."
		exit;
	fi 
	echo "/usr/local/genome/_UMItools/AGeNT_2.0.5/agent/agent.sh trim -v2 -fq1 rawFASTQ/${aa}_R1_001.fastq.gz -fq2 rawFASTQ/${aa}_R2_001.fastq.gz -out_loc trimUMIFASTQ"
	/usr/local/genome/_UMItools/AGeNT_2.0.5/agent/agent.sh trim -v2 -fq1 rawFASTQ/${aa}_R1_001.fastq.gz -fq2 rawFASTQ/${aa}_R2_001.fastq.gz -out_loc trimUMIFASTQ
done


# trim
# Mandatory Parameters:
# -fq1 filename	Read1 FASTQ file (Multiple files can be provided comma separated)
# -fq2 filename	Read2 FASTQ file (Multiple files can be provided comma separated)
# note: program outputs results of multiple comma seperated files into a single file

# At least one of these switch is mandatory:
# -halo	for a Haloplex sample
# -hs 	for a Haloplex HS sample
# -xt 	for a SureSelect XT sample
# -v2 	for a SureSelect XT HS V2 sample
# -qxt	for a SureSelect QXT sample

# Other optional parameters:
# -minFractionRead    Minimum read length as a fraction of original read length after trimming. Default is 30. Range is 0-99.
# -IDEE_FIXE	(to fix old style IDs in pairs read_id/1 read_id/2)
# -out_loc	comma seperated list of directories for output files.


