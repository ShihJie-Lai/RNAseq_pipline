#!/bin/bash

if [ $# != 4 ]; then
echo "Usage: $0 sourceFolder targetFolder R1_FASTQname baseNumber ";
echo "i.e. $0 rawFASTQ_ori rawFASTQ F8W_NTS_S76_R1_001.fastq.gz 112343232";
echo "i.e. $0 rawFASTQ_ori rawFASTQ WKY_PBMC_S77_R1_001.fastq.gz 234324324 ";
exit;
fi

sR1=$3;
if [ -e "$1/$sR1" ] ; then # Check whether file exists.

	sR2=$sR1;
	sR2=$(echo $sR1 | sed 's/_R1_001/_R2_001/')
else
	echo "$1/$sR1 not found!!!"
	exit;
fi


if [ ! -e "$1/$sR2" ] ; then # Check whether file exists.

	sR2=$sR1;
	sR2=$(echo $sR1 | sed 's/_R1\.fastq/_R2.fastq/')

fi

if [ ! -e "$1/$sR2" ] ; then # Check whether file exists.

        echo "$1/$sR2 not found!!!"
        exit;
fi

mkdir -p $2

tR1="$2/$sR1"
tR2="$2/$sR2"
Bnum=$4

echo "reformat.sh sampleseed=11 in1=$1/$sR1 in2=$1/$sR2 out1=$tR1 out2=$tR2 samplebasestarget=$Bnum"
reformat.sh sampleseed=11 in1=$1/$sR1 in2=$1/$sR2 out1=$tR1 out2=$tR2 samplebasestarget=$Bnum &> $sR1.extract.log 


