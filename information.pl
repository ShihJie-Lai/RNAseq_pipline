sub infor_RNA{
	my ($loc,$localPATH)= @_;
	my $ver="1.0.0";
	my $verDate="2020-9-1";
	my $scriptTitle="HISAT2-RNASeq-Pipeline v.$ver";
	my $verinfo="Written by SeanChang; Modified by Keh-Ming Wu, Ver $ver, $verDate";
	my $Usage =	
		"$scriptTitle\n".
		"=======================================\n".
		"$verinfo\n".
		"---------------------------------------\n".
		"* Something Wrong with the parameters *\n".
		"---------------------------------------\n".
		"Usage:\n\n".
		"ulimit -u 9999 -n 2048\n".
		"==From raw FASTQ.gz in ./ folder:\n".
		"mkdir 1_rawFASTQ 2_trimmedFASTQ 3_Run\n".
		"Trimmomatic_PE.pl RNAseq\n".
		"mv *_001.fastq.gz 1_rawFASTQ/.\n".
		"mv *.fastq.gz 2_trimmedFASTQ/.\n".
		"cd 3_Run\n".
		"ln -s ../2_trimmedFASTQ/*.fastq.gz .\n".
		"vi groupCompareInfo.txt\n".
		"[option] vi sample2GroupInfo.txt\n\n".
		"==Check Information of the POnumber==\n".
		"$loc -pc POnumber\n".
		"==List all POs with status execute or closed==\n".
		"$loc -pl execute|closed\n".
		"==Perform RNAseq Analysis==\n".
		"$loc -s species -l lib-strandness -n novelAssembly -t threadNum -f FoldChangeCutOff -p pValueCutOff -e expressonvalue [1:FPKM|0:TPM] -r re cal value -pc poinforcheck:Ponum. -gf geneFusion:Y|N -as alternativeSplicing:Y|N -sP seperateReport:Y|N\ -af arribageneFusionTag:Y|N -sn GSTKsnpTag:Y|N -cr circRNA_CIRCexplorer2:Y|N -gfF geneFusion_FFPE:Y|N \n".
		"$loc -s [*human] -l [*fr-firststrand] -n [*N] -t [*24] -f [*2] -p [*0.05] -lk [*1] -gf [*N] -as [*N] -sP [*N]\n".
		"Species:\n".
		"(1) H. sapiens (2) M. musculus (3) R. norvegicus (4) Arabidopsis (5) C.elgean (6) ZebraFish (7) Fruit Fly (8) Rice (9) Custom\n".
		"Lib-strandness: \n".
		"(1) fr-firststrand	 (2) fr-secondstrand  (3) fr-unstranded\n".
		"FoldChange>=2 or 1.5\n".
		"P-value: 0.05\n".
		"LibraryKit: \n".
		"(1) Agilent's SureSelect Strand-Specific RNA Library Preparation Kit (2) Illumina's TruSeq Stranded Total RNA Library Prep Gold Kit\n".
		"i.e.\n".
		"For Human:\n".
		"$loc -s 1 -l 1 -n N -f 2 -p 0.05 -e 0 -r 0 [-gf Y -as Y -sP Y]\n".
		"$loc -s 1 -l 1 -n N -f 2 -p 0.05 -e 0 -r 0 -gf Y -as Y -sP Y\n".
		"For Custom reference:\n".
		"$loc -s 9 -l 1 -n N -f 2  -p 0.05 -e 0 -r 0\n".
		"> Please write sample to group info. in 'sample2GroupInfo.txt' (SampleName\\tGroup).\n".
		"> Please write sample comparison info. in 'groupCompareInfo.txt' (-Treat-\\t-Ctrl-).\n".
		"> For permanant Report Information, please write in 'Report_Information.txt' or it will produce automatically.\n";
	my $POnum;
	if($localPATH=~/(C?(M|NG)SA?([0-9]{5,}))/){
		$POnum=$1;
		print "POnum of this case: $POnum\n";
	}else{
		die "POnumber not recognized from localPATH information!\n\n";
	}
	print "$Usage";
	return($POnum);
}

sub getCurrentTimeStr
{
    my ($loc)= @_;
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    my $nice_timestamp;
	if ($loc==0){
		$nice_timestamp = sprintf ( "%04d%02d%02d%02d%02d%02d",
                                   $year+1900,$mon+1,$mday,$hour,$min,$sec);
	}else{
		$nice_timestamp = sprintf ( "%04d-%02d-%02d %02d:%02d:%02d",
                                   $year+1900,$mon+1,$mday,$hour,$min,$sec);
	}
	return ($nice_timestamp);
}


return 1;
