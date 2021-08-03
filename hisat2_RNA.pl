sub runHiSat2
{
	my ($FQ1, $FQ2, $sampleName, $hisat2PATH,$HISAT2_index,$samtoolsPATH,$hisat2StrandPE,$hisat2StrandSE,$fhRUN) = @_;
	my $Sam=$sampleName.".sam";
	my $Bam=$sampleName.".bam";
	my $SBam=$sampleName.".sorted.bam";
	my $NSFile=$sampleName.".nsfile";
	my $mapping_log=$sampleName."_mapping.log";
	my $Threads=24;
	
	print " Mapping on sample: $sampleName\n";
	print $fhRUN " Mapping on sample: $sampleName\n";
	my $cmdSTr;
	if ($FQ2=~/^SE$/) #SE
	{
		if ($hisat2StrandSE eq "") # default unstranded
		{
			$cmdStr="$hisat2PATH/hisat2 -x $HISAT2_index -U $FQ1 -S $Sam --phred33 --dta --novel-splicesite-outfile $NSFile -p $Threads > $mapping_log 2>&1";
		}
		else
		{
			$cmdStr="$hisat2PATH/hisat2 -x $HISAT2_index -U $FQ1 -S $Sam --phred33 --dta --rna-strandness $hisat2StrandSE --dta --novel-splicesite-outfile $NSFile -p $Threads > $mapping_log 2>&1";
		}
	}
	else 
	{
		if ($hisat2StrandPE eq "") # default unstranded
		{
			$cmdStr="$hisat2PATH/hisat2 -x $HISAT2_index -1 $FQ1 -2 $FQ2 -S $Sam --phred33 --dta --novel-splicesite-outfile $NSFile -p $Threads > $mapping_log 2>&1";
		}
		else
		{
			$cmdStr="$hisat2PATH/hisat2 -x $HISAT2_index -1 $FQ1 -2 $FQ2 -S $Sam --phred33 --rna-strandness $hisat2StrandPE --dta --novel-splicesite-outfile $NSFile -p $Threads > $mapping_log 2>&1";
		}
	}
	print $fhRUN " $cmdStr\n";
	print " $cmdStr\n";
	if (!-e $Bam && !-e $SBam){
		system ($cmdStr);
	}
	&runSam2sortedBam($Sam, $Bam, $SBam,$mapping_log,$samtoolsPATH, $fhRUN);

}
sub runSam2sortedBam
{
	my ($Sam, $Bam, $SBam,$mapping_log,$samtoolsPATH, $fhRUN) = @_;
	my $Threads=24;
	my $cmdStr="$samtoolsPATH/samtools view -bS $Sam >$Bam -\@ $Threads";
	print $fhRUN " $cmdStr\n";
	print " $cmdStr\n";
	my $cmdStr1="$samtoolsPATH/samtools sort -\@ $Threads  $Bam  -o $SBam";
	print $fhRUN " $cmdStr1\n";
	print " $cmdStr1\n";

	if (!-e $SBam){
		if (!-e $Bam){	
			system ("$cmdStr");	
		}
		if(-e $Bam){
			system ("$cmdStr1");
		}
		print $fhRUN " rm -rf $Sam && rm -rf $Bam\n";
		system ("rm -rf $Sam");system ("rm -rf $Bam");
		my $cmdStr2="$samtoolsPATH/samtools index $SBam";
		system ("$cmdStr2");
	}
	else{
		print " Sorted Bam file found: $SBam\n";
		print $fhRUN " Sorted Bam file found: $SBam\n";
	}
	
	
	if (! -e $mapping_log){	
		my $flagstatStr=`samtools flagstat $SBam`;
		my ($totalReadNum,$pairMateMappedNum,$singletonNum,$mappingRate)=();
		if ($flagstatStr=~/([0-9]+) \+ [0-9]* paired in sequencing/){$totalReadNum=$1;}
		if ($flagstatStr=~/([0-9]+) \+ [0-9]* with itself and mate mapped/){$pairMateMappedNum=$1;}
		if ($flagstatStr=~/([0-9]+) \+ [0-9]* singletons/){$singletonNum=$1;}
		$mappingRate=sprintf "%.2f%% overall alignment rate",($pairMateMappedNum+$singletonNum)/$totalReadNum *100;
		`echo $mappingRate > $mapping_log`;			
		print " Perform flagstat => $mappingRate\n";
		print $fhRUN " Perform flagstat => $mappingRate\n"; 
#48505970 + 0 paired in sequencing
#24252985 + 0 read1
#24252985 + 0 read2
#39422 + 0 properly paired (0.08%:-nan%)
#45558228 + 0 with itself and mate mapped
#1444121 + 0 singletons (2.98%:-nan%)
	}
}

return 1;
