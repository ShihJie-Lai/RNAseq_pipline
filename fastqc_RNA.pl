sub fastq_qual{
	my ($fastqcPATH,$fastqcrevPATH,$Paired_Check,@Sample_Name)= @_;
	if (!-d "fastqc"){system("mkdir -p fastqc");}
	my @html=<fastqc/*.html>;
	my @revhtml=<fastqc/rev_*.html>;
	my $fastq_Num;
	my $Threads=24;
	my $hostname=`hostname`; 
	if ($Paired_Check==0){
		$fastq_Num =scalar(@Sample_Name);
	}else{
		$fastq_Num =scalar(@Sample_Name)*2;
	}
	print "$fastq_Num\n";
	if ($fastq_Num == scalar(@html) && scalar(@revhtml)==0 )
	{
		chdir ("./fastqc");
		if ($hostname eq "SuperMicro"){system("$fastqcrevPATH/fastqc_rev");} #是不是要統一
		else {system("perl $fastqcrevPATH/fastqc_rev.pl");}
		
		chdir ("../");
	}
	elsif ($fastq_Num != (scalar(@html) - scalar(@revhtml)))
	{
		system("rm -rf fastqc && mkdir -p fastqc");
		$cmdStr="$fastqcPATH/fastqc -t $Threads -o fastqc *.fastq*";
		print "$cmdStr\n";
		system("$cmdStr");
	
		chdir ("./fastqc");
		if ($hostname eq "SuperMicro"){system("$fastqcrevPATH/fastqc_rev");}
		else {system("perl $fastqcrevPATH/fastqc_rev.pl");}
		
		chdir ("../");
	}
	@revhtml=<fastqc/rev_*.html>;
	if ($fastq_Num != scalar(@revhtml))
	{
		system("rm -rf fastqc && mkdir -p fastqc");
		$cmdStr="$fastqcPATH/fastqc -t $Threads *.fastq*";
		print "$cmdStr\n";
		system("$cmdStr");
		
		if ($hostname eq "SuperMicro"){system("$fastqcrevPATH/fastqc_rev");}
		else {system("perl $fastqcrevPATH/fastqc_rev.pl");}
		system ("mv *.html fastqc");
		system ("mv *.zip fastqc");
	}
	else
	{
		foreach my $revhtml (@revhtml)
		{
			my $tailLine=`tail -n 1 $revhtml`;
			if ($tailLine!~ /<\/html>/) {die "!! Something wrong with '$revhtml' since the last line is not </html> tag.\n\n";}
		}	
	}
	return $Paired_Check;
}

sub baseCount{
	my ($reformatPATH)= @_;
	my $Paired_Check=0;
	my @FASTQ_All = <*.fastq*>;
	#if ($#FASTQ_All ==0){exit;}
	foreach (@FASTQ_All){
		if ($_=~/_R2.fastq/) {
			$Paired_Check=1;
			last;
		}    
	}
	my $fastq_Num=scalar(@FASTQ_All);
	if ($Paired_Check==1) 
	{
		print "A pair-ends run.\n";
	}
	else
	{
		if(@FASTQ_All){
			print "A single-ends run.\n";
		}else{
			print "Bam file only.\n";
		}
	}
	
	my %rawRBstatFiles;
	my %trimRBstatFiles;
	my %rawReadBaseCount;
	my %trimmedReadBaseCount;
	my @tempFQ=<../Raw*/*.fastq* ../raw*/*.fastq* ../*rawFASTQ/*.fastq*>;
	if(@tempFQ){
		for(my $i=0;$i<=$#tempFQ;$i++){
			my $sampleName="";
			my $readType="";
			if ($tempFQ[$i]=~/^(.+)_S[0-9]+(_L00[1-8])?_(R[12])(_[0-9]+)?\.fastq/){
				$sampleName=$1; $readType=$3; $sampleName=~ s/^\.\.\/.+\///; 
			}
			elsif ($tempFQ[$i]=~/^(.+)_(R[12])(_[0-9]+)?\.fastq/){
				$sampleName=$1; $readType=$2; $sampleName=~ s/^\.\.\/.+\///;
			}
			elsif ($tempFQ[$i]=~/^(.+)_([12])(_[0-9]+)?\.fastq/){
				$sampleName=$1; $readType=$2; $sampleName=~ s/^\.\.\/.+\///;
			}
			elsif ($tempFQ[$i]=~/^(.+)\.fastq/){
				$sampleName=$1; $sampleName=~ s/^\.\.\/.+\///;
			}
			if ($readType !~ /2/){
				$rawRBstatFiles{$sampleName}="$sampleName.rawReadBases.txt";
				if ($Paired_Check==1) {$cmdStr="$reformatPATH/reformat.sh in1=$tempFQ[$i] in2=$tempFQ[$i+1] &> $rawRBstatFiles{$sampleName}";}
				else {$cmdStr="$reformatPATH/reformat.sh in=$tempFQ[$i] &> $rawRBstatFiles{$sampleName}";}
				if (!-e $rawRBstatFiles{$sampleName})
				{
					print "$cmdStr\n";
					system("$cmdStr");
				}
			}
		}
	}
	my @Sample_Name=();
	my @Group_Name=();
	my %sample2fastq;
	for(my $i=0;$i<=$#FASTQ_All;$i++){
		if ($FASTQ_All[$i]=~/^(.+)_(R[12])\.fastq/){
			my $sampleName=$1;
			$sample2fastq{$sampleName}{$2}=$FASTQ_All[$i];
			print "$sample2fastq{$sampleName}{$2}\t$sampleName\t$2\n";
			push (@Sample_Name, $sampleName) if ($2 eq "R1");
		}elsif ($FASTQ_All[$i]=~/^(.+)\.fastq/){
			my $sampleName=$1;
			$sample2fastq{$sampleName}{"R1"}=$FASTQ_All[$i];
			push (@Sample_Name, $sampleName);
		}
	}
	for (my $i=0;$i<=$#Sample_Name;$i++){
		my $RBstatFile="$Sample_Name[$i].rbStat.txt";
		if(-e "$Sample_Name[$i].rawReadBases.txt"){
			open rawReadBaseCount,"<","$Sample_Name[$i].rawReadBases.txt";
			my @tmp=<rawReadBaseCount>;
			close rawReadBaseCount;
			foreach (@tmp){
				if ($_=~/Input:\s+([0-9]+)\sreads\s+([0-9]+) bases/){
					$rawReadBaseCount{$Sample_Name[$i]}{'reads'}=$1;
					$rawReadBaseCount{$Sample_Name[$i]}{'bases'}=$2;
					#print "$1\n";
					print "$Sample_Name[$i]: $rawReadBaseCount{$Sample_Name[$i]}{'reads'} reads, ".$rawReadBaseCount{$Sample_Name[$i]}{'bases'}." bases.\n";
				}
			}
		}else{
			$rawReadBaseCount{$Sample_Name[$i]}{'reads'}="-";
			$rawReadBaseCount{$Sample_Name[$i]}{'bases'}="-";
		}
		my $Tr_RBstatFile="$Sample_Name[$i].trimmedReadBases.txt";
		if (!-e $Tr_RBstatFile ){
			if ($Paired_Check==1) {$cmdStr="$reformatPATH/reformat.sh in1=$sample2fastq{$Sample_Name[$i]}{R1} in2=$sample2fastq{$Sample_Name[$i]}{R2} &> $Tr_RBstatFile";}
			else {$cmdStr="$reformatPATH/reformat.sh in=$sample2fastq{$Sample_Name[$i]}{R1} &> $Tr_RBstatFile";}
			print "$cmdStr\n";
			system("$cmdStr");
		}
		if(!-e $RBstatFile){
			open trimmedReadBaseCount,"<","$Tr_RBstatFile";
			my @tmp=<trimmedReadBaseCount>;
			close trimmedReadBaseCount;
			foreach (@tmp){
				if ($_=~/Input:\s+([0-9]+) reads\s+([0-9]+) bases/){
					$trimmedReadBaseCount{$Sample_Name[$i]}{'reads'}=$1;
					$trimmedReadBaseCount{$Sample_Name[$i]}{'bases'}=$2;
					print "$Sample_Name[$i]: $trimmedReadBaseCount{$Sample_Name[$i]}{'reads'} reads, $trimmedReadBaseCount{$Sample_Name[$i]}{'bases'} bases.\n";
					open RBstatFile,">","$RBstatFile";
					print RBstatFile "$Sample_Name[$i]\t$rawReadBaseCount{$Sample_Name[$i]}{'reads'}\t$rawReadBaseCount{$Sample_Name[$i]}{'bases'}\t$trimmedReadBaseCount{$Sample_Name[$i]}{'reads'}\t$trimmedReadBaseCount{$Sample_Name[$i]}{'bases'}\n";
					close RBstatFile;
				}
			}
		}
	}
	
	return ($Paired_Check,\%sample2fastq,\@Sample_Name);
}

return 1;
