#!/usr/bin/perl -w
use Cwd;
print "perl rand_trim.pl down_limit up_limit ex perl rand_trim.pl 8 10\n";
sub rand_n {
	my ($down_linm,$up_lime) = @_;
	my $random_number=0;
	while ($random_number <$down_linm ){
    $random_number=int (rand $up_lime);
	}
	my $random_number1 =rand();
	$random_number1=($random_number+$random_number1)*(10**9);
    return $random_number1;
};
my $localPATH=getcwd;
my $localPATH_ori=$localPATH."/rawFASTQ_ori";
chdir("$localPATH_ori");
my @raw=<*_R1_*>;
chdir("$localPATH");
foreach (@raw){
my $a1=int(&rand_n($ARGV[1],$ARGV[2]));
my $aa="extractPEfastqBybaseNum.sh rawFASTQ_ori rawFASTQ $_ $a1";
print $aa."\n";
system("$aa");
}
#print "asa$ARGV[0]\n";
if($ARGV[0]==1){
	#print "ass22\n";
	system("useAgentTrimHS2_all.sh");
}


my $localPATH_raw=$localPATH."/trimUMIFASTQ";
if(!-d $localPATH_raw){
	$localPATH_raw=$localPATH."/rawFASTQ";
}

chdir("$localPATH_raw");
my $errMsg;

system("Trimmomatic_PE.pl RNAseq");
my $errMsg;
$errMsg=`mkdir -p $localPATH/trimmedFASTQ`;
$errMsg=`mv *_R1.fastq.gz $localPATH/trimmedFASTQ`;
$errMsg=`mv *_R2.fastq.gz $localPATH/trimmedFASTQ`;
$errMsg=`md5sum *.gz > $localPATH_raw/md5.txt&`;

chdir("$localPATH/trimmedFASTQ");
$errMsg=`md5sum *.gz > $localPATH/trimmedFASTQ/md5.txt&`;
