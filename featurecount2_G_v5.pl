#!/usr/bin/perl -w
use Cwd;
use strict;
no strict 'refs';
use Encode qw/encode decode/;
require "/home/shihjielai/RNAseq_PIP/excelO.pl";
require "/home/shihjielai/RNAseq_PIP/StringTie_R.pl";
require "/home/shihjielai/RNAseq_PIP/StringTie_FPKM_TPM.pl";
require "/home/shihjielai/RNAseq_PIP/setting_RNA.pl";
require "/home/shihjielai/RNAseq_PIP/fastqc_RNA.pl";
require "/home/shihjielai/RNAseq_PIP/hisat2_RNA.pl";
require "/home/shihjielai/RNAseq_PIP/stringtie_pip.pl";
require "/home/shihjielai/RNAseq_PIP/information.pl";
require "/home/shihjielai/RNAseq_PIP/AS_GF.pl";
require "/home/shihjielai/RNAseq_PIP/enrichment_analysis.pl";
use Getopt::Long qw(GetOptions);

my $localPATH=`pwd`;
my $POnum=&infor_RNA($0,$localPATH);
my ($POinfoCheck,$Species)=();
my ($FC_Cut,$PValue,$expr,$Strandnesss,$re_report)=(2,0.05,0,1,0);
my ($Novel_option,$Threads)=("N",24);
my ($geneFusionTag,$alternativeSplicingTag,$seperateReport,$arribageneFusionTag,$GSTKsnpTag,$RunDE,$circRNA_CIRCexplorer2,$geneFusionTag_FFPE,$circ_RNA_CIRIquant,$rMATs)=("N","N","N","N","N","Y","N","N","N","N");

GetOptions(
		's=s'  => \$Species,
		'l=s'  => \$Strandnesss,
		'n=s'  => \$Novel_option,
		'f=s'  => \$FC_Cut,
		'p=s'  => \$PValue,
		'e=s'  => \$expr, ####1:FPKM or 0:TPM
		'r=s'  => \$re_report, #### 1: Not re cal p-value & 0: cal p-value 
		'gf=s'  => \$geneFusionTag,
		'as=s'  => \$alternativeSplicingTag,
		'sP=s'  => \$seperateReport,
		'af=s'  => \$arribageneFusionTag,
		'sn=s'  => \$GSTKsnpTag,
		'pc=s'  => \$POinfoCheck,
		'cr=s'	=> \$circRNA_CIRCexplorer2,
		'gfF=s' => \$geneFusionTag_FFPE,
		'ciri=s'	=> \$circ_RNA_CIRIquant,
		'rMATs=s' => \$rMATs,
		'Run=s'  => \$RunDE
);

if ( !defined($Species)){exit;}
open ( RUN, ">", "run".&getCurrentTimeStr(0).".log");
print RUN "featurecount2_G_v5.pl -s $Species -l $Strandnesss -f $FC_Cut -p $PValue -e $expr -r $re_report -gf $geneFusionTag -as $alternativeSplicingTag -sP $seperateReport -af $arribageneFusionTag -sn $GSTKsnpTag -cr $circRNA_CIRCexplorer2 -gfF $geneFusionTag_FFPE -ciri $circ_RNA_CIRIquant -Run $RunDE\n";
my ($hash_set_strand,$hash_set_path,$hash_set_genome)=&setting("$Strandnesss","$Species");

print RUN &getCurrentTimeStr(1)."write reportInfo\n";
if (!-e "${POnum}_reportInfo.txt" || defined($POinfoCheck)){
	&getLIMSinfo($POnum, ${$hash_set_genome}{ensemblDBname} ,$POinfoCheck); # write $reportInfo="${POnum}_reportInfo.txt";
}

###############################################
print RUN &getCurrentTimeStr(1)."cal base count\n";
my ($Paired_Check,$sample2fastq,$Sample_Name)=&baseCount(${$hash_set_path}{reformatPATH});
print RUN &getCurrentTimeStr(1)."run fastqc";
&fastq_qual(${$hash_set_path}{fastqcPATH},${$hash_set_path}{fastqcrevPATH},"$Paired_Check",@{$Sample_Name});
print RUN &getCurrentTimeStr(1)."run HiSat2 & Bam\n";
foreach (@$Sample_Name){
	if( $Paired_Check==0){
		&runHiSat2( "${$sample2fastq}{$_}{R1}","SE","$_","${$hash_set_path}{hisat2PATH}","${$hash_set_genome}{HISAT2_index}","${$hash_set_path}{samtoolsPATH}","${$hash_set_strand}{hisat2StrandPE}","${$hash_set_strand}{hisat2StrandSE}",*RUN )
	}else{
		&runHiSat2( "${$sample2fastq}{$_}{R1}","${$sample2fastq}{$_}{R2}","$_","${$hash_set_path}{hisat2PATH}","${$hash_set_genome}{HISAT2_index}","${$hash_set_path}{samtoolsPATH}","${$hash_set_strand}{hisat2StrandPE}","${$hash_set_strand}{hisat2StrandSE}",*RUN )
	}
}

print RUN &getCurrentTimeStr(1)."run StringTie & upload data\n";
&gdupload("$localPATH","$POnum");
&stringtie_RNA(${$hash_set_genome}{GFFGTF},${$hash_set_strand}{stringtieStrand},${$hash_set_path}{stringtiePATH},$0);

###############################################

my %hash_a=&Stringtie_gene("$expr");
my %hash_t=&Stringtie_Transcript2();
if($RunDE=~/^(Y)$/i){
&tabtoExcel("gene_count_matrix.csv","./Expression_Table/Gene_Count.xlsx");
&tabtoExcel("transcript_count_matrix.csv","./Expression_Table/TS_Count.xlsx");
&tabtoExcel("rawdata.txt","./Expression_Table/Gene_TPM.xlsx");
&tabtoExcel("T_rawdata.txt","./Expression_Table/TS_TPM.xlsx");
}
##############################################
my %hash_s=();
if (-e "sample2GroupInfo.txt")
{
	open(IN4, "<", "sample2GroupInfo.txt");
	my @IN4=<IN4>;close (IN4);
	foreach (@IN4){
		chomp $_;
		my @Temp=split /\t/, $_;
		push @{$hash_s{$Temp[1]}}, $Temp[0];
	}
}
########################################################
if($RunDE=~/^(Y)$/i){
	print RUN &getCurrentTimeStr(1)."run DE analysis\n";
	open(IN2, "<", "groupCompareInfo.txt") or die "Could not open \n\n";
	my @IN2=<IN2>;close (IN2);
	foreach my $ab(@IN2){
		chomp $ab;
		my @Temp2=split /\t/, $ab;
		my $nameR;
		my $nameR_T;
		if ( ($hash_s{$Temp2[0]} || $hash_s{$Temp2[1]}) && $re_report==0 ){
			print RUN "DESeq2 --$ab--$FC_Cut--$PValue\n";
			$nameR=&StringTie_R_DESeq2("$ab","$FC_Cut","gene_count_matrix.csv","rawdata.txt",0,"$PValue",%hash_s);# ab:compare; FC_cut:FC ,0 :gene, hash_s: group
			$nameR_T=&StringTie_R_DESeq2("$ab","$FC_Cut","transcript_count_matrix.csv","T_rawdata.txt",1,"$PValue",%hash_s); # ab:compare; FC_cut:FC ,1 :transcript, hash_s: group
		}elsif(!($hash_s{$Temp2[0]} || $hash_s{$Temp2[1]}) && $re_report==0) {
			print RUN "DESeq --$ab--$FC_Cut--$PValue\n";
			$nameR=&StringTie_R_DESeq("$ab","$FC_Cut","gene_count_matrix.csv","rawdata.txt",0,"$PValue");# ab:compare; FC_cut:FC ,0 :gene
			$nameR_T=&StringTie_R_DESeq("$ab","$FC_Cut","transcript_count_matrix.csv","T_rawdata.txt",1,"$PValue");# ab:compare; FC_cut:FC ,1 :transcript
		}else{
			$nameR=&StringTie_R_RE("$ab","$FC_Cut",0,"$PValue",%hash_s);
			$nameR_T=&StringTie_R_RE("$ab","$FC_Cut",1,"$PValue",%hash_s);
		}
		my $inlist=&outexcel(${$hash_set_genome}{Biomart_annotation},"${nameR}_DE1.tab","${nameR}_FC${FC_Cut}X1.tab","$nameR",0,"$expr","$FC_Cut",%hash_a);
		#biomart: annotation;DE1:all DE; FCX1: DE FCX;nameR:excel name; 0:gene;expr:exp value;FC:FC;hash_a:annotation locus
		&outexcel(${$hash_set_genome}{Biomart_annotation},"${nameR_T}_DE1.tab","${nameR_T}_FC${FC_Cut}X1.tab","$nameR_T",1,0,"$FC_Cut",%hash_t);# 1: transcript; 1:FPKM
		&clusterProfiler_GO_KEGG_ENSG_ID("$inlist", "${$hash_set_genome}{OrgDB}","$nameR","${$hash_set_genome}{KEGG_symbol}","$FC_Cut");#inlist: FC gene lsit
		&clusterProfiler_GO_KEGG_ENSG_ID2("${$hash_set_genome}{OrgDB}","$nameR","${$hash_set_genome}{KEGG_symbol}","$FC_Cut");#inlist: FC gene lsit
	}
#}	
	if($Species!='9'){
		print RUN &getCurrentTimeStr(1)."run DE GO & KEGG\n";
		&enrich_ANA();}
}
	&StringTie_R_plot("$expr");


###############################################

print RUN &getCurrentTimeStr(1)."run AS or GF\n";
if ($alternativeSplicingTag=~/^(Y)$/i){
	print RUN &getCurrentTimeStr(1)."alternative Splicing\n";
	&miso_AS(${$hash_set_path}{misoPATH},${$hash_set_genome}{misoGFFindex});
}
if ($geneFusionTag=~/^(Y)$/i){
	print RUN &getCurrentTimeStr(1)."Gene Fusion\n";
	&STAR_GF(${$hash_set_genome}{starFusionDBindex},${$hash_set_path}{starFusionPATH},$Paired_Check,%{$sample2fastq});
}
if ($arribageneFusionTag=~/^(Y)$/i){
	print RUN &getCurrentTimeStr(1)."arriba Gene Fusion\n";
	&arriba_GF(${$hash_set_genome}{starFusionDBindex},${$hash_set_path}{arriba},%{$sample2fastq});
}
if ($GSTKsnpTag=~/^(Y)$/i){
	print RUN &getCurrentTimeStr(1)."GATK SNP\n";
	&GATK_snp("${$hash_set_path}{samtoolsPATH}","${$hash_set_genome}{starFusionDBindex}");
}
if ($circRNA_CIRCexplorer2=~/^(Y)$/i){
	print RUN &getCurrentTimeStr(1)."CIRCexplorer2\n";
	&circ_RNA_CIRCexplorer2("${$hash_set_genome}{starFusionDBindex}");
}
if ($geneFusionTag_FFPE=~/^(Y)$/i){
	print RUN &getCurrentTimeStr(1)."Gene Fusion\n";
	&STAR_GF_FFPE(${$hash_set_genome}{starFusionDBindex},${$hash_set_path}{starFusionPATH},${$hash_set_path}{STAR},${$hash_set_path}{arriba},%{$sample2fastq});
}
if ($circ_RNA_CIRIquant=~/^(Y)$/i){
	print RUN &getCurrentTimeStr(1)."CIRCexplorer2\n";
	&circ_RNA_CIRCexplorer2("${$hash_set_genome}{starFusionDBindex}",%{$sample2fastq});
}
if ($rMATs=~/^(Y)$/i){
	print RUN &getCurrentTimeStr(1)."rMATs\n";
	&rmats_AS("${$hash_set_genome}{rmats}","${$hash_set_path}{STAR}","${$hash_set_genome}{starFusionDBindex}",%hash_s);
}


###############################################

print RUN &getCurrentTimeStr(1)."run E-report\n";
system("RNAseq_E-report_202008.pl $POnum _FC${FC_Cut}X $PValue $FC_Cut");
close RUN;
