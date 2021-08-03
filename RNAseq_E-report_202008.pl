#!/usr/bin/perl -w
use strict;
use DBI;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure(qw(posix_default no_ignore_case));
use FileHandle;
use utf8;
binmode(STDIN, ':encoding(utf8)');
binmode(STDOUT, ':encoding(utf8)');
binmode(STDERR, ':encoding(utf8)');


STDOUT->autoflush(1);

my $scriptTitle="Script to produce RNAseq E-report";
my $verinfo="Written by SeanChang\nModified by Keh-Ming Wu, Ver 1.0.8, 2019.08.08";
my $Usage =	
	"$scriptTitle\n".
	"=======================================\n".
	"$verinfo\n".
	"---------------------------------------\n".
	"* Something Wrong with the parameters *\n".
	"---------------------------------------\n\n".
	"Usage:\n".
	"$0 POnum\n".
	"i.e. $0 NGS1070818 _FC2X 1\n\n";
	
	
my $POnum=$ARGV[0] || "";
my $outStr=$ARGV[1] || ""; if ($outStr!~/^_/){$outStr="_".$outStr;}   # _FC2X
my $Pvalue=$ARGV[2];
my $FC=$ARGV[3];



my $libKit; 
# libKit: (1) Agilent's SureSelect Strand-Specific RNA Library Preparation Kit (2) Illumina's TruSeq Stranded Total RNA Library Prep Gold Kit

my $localPATH=`pwd`;
my $seperateReportTag;
if (defined($ARGV[2])){if ($ARGV[2]=~/^seperateReport$/i){$seperateReportTag="Y";}}

if (!defined($POnum) || $POnum eq "")
{
	if($localPATH=~/(C?(M|NG)SA?([0-9]{5,}))/)
	{
		$POnum=$1;
		print "POnum of this case: $POnum\nCorrect? (Y/N) ";
		chomp( my $POck = <STDIN> );
		if ($POck=~/N/i)
		{
			print "Please enter the POnum: ";
			chomp( $POnum = <STDIN> );
		}
	}
	else
	{
		print "POnumber not recognized from localPATH information!\n\n";
	}
}

my @FASTQ_All=<*.fastq*>;
my @bam_All=<*.bam>;
##### bam only #####
my @Sample_Name;
my %sample2GroupInfo=();

	my $S="";
	foreach my $aa (@bam_All){
		print "$aa\n";	
		my @Temp=split /\./,$aa;
		$S=&trim($Temp[0]);
		push (@Sample_Name, $S);
		print "@Sample_Name\n\n";
		my $mapping_log="$Temp[0]_mapping.log";
		print "$mapping_log";
		if (! -e $mapping_log){	
			my $flagstatStr=`samtools flagstat $aa`;
			my ($totalReadNum,$pairMateMappedNum,$singletonNum,$mappingRate)=();
			if ($flagstatStr=~/([0-9]+) \+ [0-9]* paired in sequencing/){$totalReadNum=$1;}
			if ($flagstatStr=~/([0-9]+) \+ [0-9]* with itself and mate mapped/){$pairMateMappedNum=$1;}
			if ($flagstatStr=~/([0-9]+) \+ [0-9]* singletons/){$singletonNum=$1;}
			$mappingRate=sprintf "%.2f%% overall alignment rate",($pairMateMappedNum+$singletonNum)/$totalReadNum *100;
			$mappingRate="$flagstatStr\n"."$mappingRate";
			open (ans2,">$mapping_log");
			print ans2 "$mappingRate";
			close ans2;
#48505970 + 0 paired in sequencing
#24252985 + 0 read1
#24252985 + 0 read2
#39422 + 0 properly paired (0.08%:-nan%)
#45558228 + 0 with itself and mate mapped
#1444121 + 0 singletons (2.98%:-nan%)
		}
	}

##### bam only #####

@FASTQ_All=sort naturalSort @FASTQ_All;
my $fastq_Num=scalar(@FASTQ_All);

my $Paired_Check=0;
my @Sample_Name_fastq;
foreach (@FASTQ_All){
	my $S="";
	if ($_=~/_R2.fastq/) {
		$Paired_Check=1;
		next;
	}elsif ($_=~/^(.+)_(R[12])\.fastq/){
		$S=&trim($1);
		push (@Sample_Name_fastq, $S) if ($2 eq "R1");
	}elsif ($_=~/^(.+)\.fastq/){
		$S=&trim($1);
		push (@Sample_Name_fastq, $S)
	}
}




my $sample2GroupInfoFile="sample2GroupInfo.txt";
my $replicateGroupNum=0;
if (-e $sample2GroupInfoFile)
{
	print getCurrentDateTime()." Found $sample2GroupInfoFile.\n";
	open S2Ginfo,"<","$sample2GroupInfoFile" || die "Could not open $sample2GroupInfoFile.\n\n";
	my @tmp=<S2Ginfo>;close S2Ginfo;
	foreach (@tmp)
	{
		$replicateGroupNum++;
		chomp $_;
		if ($_=~/^Sample.*\sGroup/i && $replicateGroupNum==1){print "====================\n$_\n"; $replicateGroupNum--; next;}
		elsif ($_=~/^(.+)\s+(.+)$/)
		{
			my $S=trim($1);
			my $G=trim($2);
			print "$replicateGroupNum. $S => $G\n";
			print "@Sample_Name\n\n";
			if ( !grep( /^$S$/, @Sample_Name ) ) {die "Typo?? Sample $S is not in the sampleList!!\n\n";}
			#if ( !(trim($1) ~~ @Sample_Name) ) {die "Typo?? Sample $1 is not in the sampleList!!\n\n";}
			
			push @{$sample2GroupInfo{$G}},$S;
			
			#if ( !grep( /^$G$/, @Group_Name ) ) {push @Group_Name,$G;}
			
			
		}
	}

}

###Get Information###
my $demoReportInfoStr=
"Organism	Homo sapiens
PONumber	NGS1070890
Institute	高雄醫學大學
Customer	許雅玲
Tel	07-3121101 ext. 2136
Sepc	150PE
Amount	9G
DBused	http://www.ensembl.org/Homo_sapiens/Info/Index";


my $reportInfo="${POnum}_reportInfo.txt";
if (-e "Report_Information.txt"){$reportInfo="Report_Information.txt";}

print "Reading $reportInfo\n";
my %hashinfor;
open(ReportInfor, "$reportInfo") or die "Please fill up $reportInfo as below (demo): \n$demoReportInfoStr\n\n";
my @ReportInfor=<ReportInfor>;close ReportInfor;
foreach (@ReportInfor){
	chomp $_;$_=~s/\r//g;
	next if ($_ eq "");
	my @Temp=split /\t/,$_;
	$hashinfor{$Temp[0]}=$Temp[1];
	
	print "  $Temp[0] => $hashinfor{$Temp[0]}\n";
}
print "===============================\n";
my $groupCompareInfoFile="groupCompareInfo.txt";
my %groupCompareInfo;
if (-e $groupCompareInfoFile)
{
	print "Found $groupCompareInfoFile.\n";
	open gCompare,"<","$groupCompareInfoFile" || die "Could not open $groupCompareInfoFile.\n\n";
	my @tmp=<gCompare>;close gCompare;
	my $i=0;
	print "   Exp(Treat) vs Ctrl\n";
	foreach (@tmp)
	{
		$i++;
		chomp $_;
		if ($_=~/^-?Treat-?\s-?Ctrl-?/i && $i==1){print "====================\n$_\n";next;}
		elsif ($_=~/^(.+)\s+(.+)$/)
		{
			my $T=trim($1);
			my $C=trim($2);
			print "$i. $T vs $C\n";
			
			$groupCompareInfo{$T."\t".$C}=$T."_vs_".$C;
		}
	}
}

else
{
	print("!!No sample\/group comparison design file \'groupCompareInfo.txt\' found!!\n\n");
}

###
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
$mon+=1;
$year+= 1900;
my $Date="$year-$mon-$mday";
### Get File ###
my $localEreportFolder="${POnum}_E-Report".$outStr;
my $reportTemplateFolder;

my $hostname=`hostname`; 
chomp $hostname;

if ($hostname eq "SuperMicro")
{
	$reportTemplateFolder="/Reports/E-Report-RNAseq-2019"; #####
}
elsif ($hostname eq "R930")
{
	$reportTemplateFolder="/home/shihjielai/RNAseq_PIP/E-Report-RNAseq-2020";
}

if (-d "./$localEreportFolder"){system("rm -rf ./$localEreportFolder");}
system("cp -rf $reportTemplateFolder ./$localEreportFolder");
system("mkdir -p ./$localEreportFolder/Content/fastqc");
system("cp -f ./fastqc/rev_*.html ./$localEreportFolder/Content/fastqc");
system("cp -r ./cummer ./$localEreportFolder/Content/.");
if (-e "groupCompareInfo.txt"){
	system("cp -rf ./R_dataFiles/* ./$localEreportFolder/Content/.");
	system("rm -rf ./$localEreportFolder/Content/*.r ./$localEreportFolder/Content/*.pathview ./$localEreportFolder/Content/*.txt");
	system("cp -r ./DE_Gene ./$localEreportFolder/Content/.");
	system("cp -r ./DE_TS ./$localEreportFolder/Content/.");
}
system("cp -r ./Expression_Table ./$localEreportFolder/Content/.");
system("cp -r ../6_Report/Gene_Fusion ./$localEreportFolder/Content/.") if (-d "../6_Report/Gene_Fusion");
system("cp -r ../6_Report/Alternative_Splicing ./$localEreportFolder/Content/.") if (-d "../6_Report/Alternative_Splicing");



#### Parser Index.html ####

print "Reading $localEreportFolder/index.html\n";

open(Index, "$localEreportFolder/index.html") or die;
my @Index_HTML=<Index>;close Index;system ("rm $localEreportFolder/index.html");
my $Index_HTML="";
open(IndexWrite, ">$localEreportFolder/index.html") or die;
foreach (@Index_HTML){
	if ($_=~/FOR_ORG1/) {
		$_=~s/FOR_ORG1/$hashinfor{"Organism"}/g;
	}elsif($_=~/FOR_PO1/){
		$_=~s/FOR_PO1/$hashinfor{"PONumber"}/g;
	}elsif($_=~/FOR_CUST_INS/){
		$_=~s/FOR_CUST_INS/$hashinfor{"Institute"}/g;
	}elsif($_=~/FOR_CUST_NAME/){
		$_=~s/FOR_CUST_NAME/$hashinfor{"Customer"}/g;
	}elsif($_=~/FOR_CUST_Tel/){
		$_=~s/FOR_CUST_Tel/$hashinfor{"Tel"}/g;		
	}elsif($_=~/FOR_DATE/){
		$_=~s/FOR_DATE/$Date/g;  
	}
	print IndexWrite $_;
	$Index_HTML.=$_;
}
close IndexWrite;

my $indexDEstr="
	<div class='row'>
		<a class='iframeCB marginDown2' href='Content/DE.html' title='Differential Expression'>
			<div class='flowTitleBox text-center marginDown1' data-scrollreveal='enter lefts'>
				<h3 class='flowTitleBoxFont'>Differential Expression</h3>
			</div>
		</a>
  		<div class='flow-line-bottom'>
  		<div class='flow-circle-bottom'>
		</div>
		</div>
			<div class='col-md-5 flow-box marginDown1'>
				<p>Differential gene/transcript expression analysis results.</p>
			</div>
	</div>
";

if ((!-d "./$localEreportFolder/Content/Gene_Fusion") && (!-d "./$localEreportFolder/Content/Alternative_Splicing"))
{
	$Index_HTML=~ s/<!----------Set_Write_Here_DE----------->.+<!----------Set_Write_Here_End----------->/$indexDEstr/s;
	open(IndexWrite, ">$localEreportFolder/index.html") or die;
	print IndexWrite $Index_HTML;
	close IndexWrite;
	
	print "NOT found Gene_Fusion & Alternative_Splicing results!!\n";
	print "Discard 'Gene_Fusion & Alternative_Splicing results' section in index.html\n\n";
	`rm -rf  ./$localEreportFolder/Content/FusionAS.html`;
}
else
{
	print "Found Gene_Fusion & Alternative_Splicing results!!\n";
	print "Keep Gene_Fusion & Alternative_Splicing results section in index.html\n\n";
}



#### Parsing Sequencing.html ####

# libkit==1   =>  "https://www.agilent.com/en/sureselect-strand-specific-rna-library-prep"
# <a data-fancybox data-type="iframe" data-src="https://www.agilent.com/en/sureselect-strand-specific-rna-library-prep" href="javascript:void(0);"><strong>Agilent's SureSelect Strand-Specific RNA Library Preparation Kit</strong></a>




if ($hashinfor{'kit'}=~/Illumina/) # For 1 as default kit, no need to change. For 2, change to "Illumina's TruSeq Stranded Total RNA Library Prep Gold Kit"
{
	print "Reading $localEreportFolder/Content/Sequencing.html\n";

	open(SeqInfor, "$localEreportFolder/Content/Sequencing.html") or die;
	my @SeqInfor=<SeqInfor>;close SeqInfor;system ("rm -rf $localEreportFolder/Content/Sequencing.html");
	open(SeqInforWrite, ">>$localEreportFolder/Content/Sequencing.html") or die;
	foreach my $line (@SeqInfor){
		if ($line=~/Agilent's SureSelect Strand-Specific RNA Library Preparation Kit/)
		{
			print SeqInforWrite '<a href="https://www.illumina.com/products/by-type/sequencing-kits/library-prep-kits/truseq-stranded-total-rna.html" target="_blank"><strong>Illumina\'s TruSeq Stranded Total RNA Library Prep Gold Kit (Cat. 20020598)</strong></a>';
			print SeqInforWrite "\n";
		}
		else
		{
			print SeqInforWrite $line;
		}
	}

}


if ($hashinfor{'kit'}=~/HS2/) # For 1 as default kit, no need to change. For 2, change to "Illumina's TruSeq Stranded Total RNA Library Prep Gold Kit"
{
	print "Reading $localEreportFolder/Content/Sequencing.html\n";

	open(SeqInfor, "$localEreportFolder/Content/Sequencing.html") or die;
	my @SeqInfor=<SeqInfor>;close SeqInfor;system ("rm -rf $localEreportFolder/Content/Sequencing.html");
	open(SeqInforWrite, ">>$localEreportFolder/Content/Sequencing.html") or die;
	foreach my $line (@SeqInfor){
		if ($line=~/Agilent's SureSelect Strand-Specific RNA Library Preparation Kit/)
		{
			print SeqInforWrite '<a href= "https://www.agilent.com/en/product/next-generation-sequencing/hybridization-based-next-generation-sequencing-ngs/rna-seq-library-preparation-kits/sureselect-xt-hs2-mrna-library-preparation-kit-1049125" target="_blank"><strong>SureSelect XT HS2 mRNA Library Preparation kit (Agilent, USA)</strong></a>';
			print SeqInforWrite "\n";
		}
		else
		{
			print SeqInforWrite $line;
		}
	}

}

close(SeqInforWrite);

#### Parsing ExpInfor.html ####

print "Reading $localEreportFolder/Content/ExpInfor.html\n";

open(ExpInfor, "$localEreportFolder/Content/ExpInfor.html") or die;
my @ExpInfor=<ExpInfor>;close ExpInfor;system ("rm -rf $localEreportFolder/Content/ExpInfor.html");
open(ExpInforWrite, ">>$localEreportFolder/Content/ExpInfor.html") or die;
foreach my $line (@ExpInfor){
	if($line=~/<!--Start_here_with_Sample-->/){
		foreach my $SS(@Sample_Name){
			my $SP=$hashinfor{"Organism"};
			my $ST=$hashinfor{"Spec"};
			my $ED=$hashinfor{"Amount"};
			print ExpInforWrite "
			<tr align=\"center\"><td width=\"28%\">$SS</td>
				<td width=\"18%\">$SP</td>
				<td width=\"18%\">$ST</td>
				<td width=\"27%\">$ED</td>
			</tr>
			";
			
			print "
			<tr align=\"center\"><td width=\"28%\">$SS</td>
				<td width=\"18%\">$SP</td>
				<td width=\"18%\">$ST</td>
				<td width=\"27%\">$ED</td>
			</tr>\n
			";
		}
	}
	elsif($line=~/<!--For_Group_Information-->/){
		if ($replicateGroupNum>0)
		{
			print ExpInforWrite "
			<table width='650' border='3'>
			<tbody>
			<tr align='center'>
			  <td colspan='2'><div align='center'><strong><font color='red'>Group Information</font></strong></div></td>
			</tr>
			<tr align='center'>
			  <td width='50%'><b><font color='red'>Group</font></b></td>
			  <td width='50%'><b><font color='red'>Sample(s)</font></b></td>
			</tr>
			";
			
			foreach my $G (sort keys %sample2GroupInfo)
			{
				my $Snum=scalar(@{$sample2GroupInfo{$G}});
				my $Sstr=join("<BR>",@{$sample2GroupInfo{$G}});
				if ($Snum>1 || $sample2GroupInfo{$G}[0] ne $G )
				{
					print ExpInforWrite
		"<tr>
			<td width='50%' align='center'>$G</td>
			<td width='50%' align='center'>$Sstr</td>
		</tr>
		";
				}
			}
		
			print ExpInforWrite "
			</tbody>
			</table>
			";
		
		}
	}
	elsif($line=~/<!--For_Exp_Design-->/){
		foreach my $Set (sort keys %groupCompareInfo)
		{
			my @Temp=split /\t/, $Set;
			print ExpInforWrite "
			<tr>
				<td width=\"50%\" align=\"center\">$Temp[0]</td>
				<td width=\"50%\" align=\"center\">$Temp[1]</td>
			</tr>";
			print "
			<tr>
				<td width=\"50%\" align=\"center\">$Temp[0]</td>
				<td width=\"50%\" align=\"center\">$Temp[1]</td>
			</tr>\n";
		}
	}
	else{
		print ExpInforWrite $line;
	}
}
close(ExpInforWrite);

###########
print "Reading $localEreportFolder/Content/Data_QC.html\n";

open(QCFile, "$localEreportFolder/Content/Data_QC.html") or die;
my @QCHTML=<QCFile>;close QCFile;
open(QCFile_Write, ">$localEreportFolder/Content/Data_QC.html") or die;

#my @Stats=<*.stats>;
#system ("cp -rf fastqc/* $localEreportFolder/Content/Data_QC.html");


foreach (@QCHTML)
{
    if ($_=~/<!----------Sample_Reads_information----------->/) 
	{
        for(my $K=0;$K<=$#Sample_Name_fastq;$K++){
            open(STATS, "$Sample_Name_fastq[$K].rbStat.txt") || next "$Sample_Name_fastq[$K].rbStat.txt not found!!\n\n";
            my @STATS=<STATS>;
            chomp $STATS[0];
            my @Temp=split /\t/, $STATS[0];
			my $serial=$K+1;
			@Temp[1..4]=map { &commify($_) } @Temp[1..4];
			
            print QCFile_Write "
                <tr>
				<td align='right'>$serial. </td>
                <td align='center'>$Temp[0]</td>
                <td align='right'>$Temp[1]</td>
                <td align='right'>$Temp[2]</td>
                <td align='right'>$Temp[3]</td>
                <td align='right'>$Temp[4]</td>
                </tr>                                                
            ";
			
			print "
                <tr>
                <td align='right'>$serial. </td>
		<td align='center'>$Temp[0]</td>
                <td align='right'>$Temp[1]</td>
                <td align='right'>$Temp[2]</td>
                <td align='right'>$Temp[3]</td>
                <td align='right'>$Temp[4]</td>
                </tr> \n                                              
            ";
        }
    }
	elsif ($_=~/<!--------------------QC------------------->/) 
	{
		
		my $serial=0;
		foreach my $SS(@Sample_Name_fastq)
		{
			my $R1=$SS."_R1_fastqc";
			my $R2=$SS."_R2_fastqc";
			
			$serial++;
			my $revHtmlR1="$localEreportFolder/Content/fastqc/rev_$R1.html";
			if (!-e $revHtmlR1)
			{
				if ($Paired_Check==0)
				{
					$R1=$SS."_fastqc";
					$revHtmlR1="$localEreportFolder/Content/fastqc/rev_$R1.html";
				}
			}
			
			
			print QCFile_Write "
			<tr>
			<td align='right'>$serial. </td>
			<td align='center'>$SS</td>
			<td align='center'>";
			
			print QCFile_Write "
			<a data-fancybox data-type=\"iframe\" data-src=\"fastqc/rev_$R1.html\" href=\"javascript:;\"><img  src=\"../assets/img/icon/html.png\" height=\"64\" width=\"64\">R1</a>
			<a data-fancybox data-type=\"iframe\" data-src=\"fastqc/rev_$R2.html\" href=\"javascript:;\"><img  src=\"../assets/img/icon/html.png\" height=\"64\" width=\"64\">R2</a>
			</td>
			</tr>" if ($Paired_Check==1);
			
			print QCFile_Write "
			<a data-fancybox data-type=\"iframe\" data-src=\"fastqc/rev_$R1.html\" href=\"javascript:;\"><img src=\"../assets/img/icon/html.png\" height=\"64\" width=\"64\"></a>
			</td>
			</tr>" if ($Paired_Check==0);
			#if (!-e $revHtmlR1){ "!! Paired_Check=$Paired_Check|$revHtmlR1 not found !!!\n\n";}
		}
		
	}else{
		print QCFile_Write $_;
	}
	
}
close(QCFile_Write);
#####################

print "Reading $localEreportFolder/Content/Align.html\n";

open(Align, "$localEreportFolder/Content/Align.html") or die;
my @Align=<Align>;close Align;

open(Align2, ">$localEreportFolder/Content/Align.html");
foreach my $line(@Align){
	if($line=~/<!--Here_is_mapping_information_insert-->/)
	{
		my @mapping=<*_mapping.log>;
		for my $log(@mapping)
		{
				open(Mlog, "$log") or die;
				my @Mlog=<Mlog>;close Mlog;
				my $lastline=pop @Mlog;
				my $sample=$log;
				$sample=~s/(1st)?_mapping.log//;
				my @Temp=split / /,$lastline;
				my $mapping_rate=$Temp[0];
				print Align2  "
					<tr>
					<td width=\"50\%\" align='center'><b>$sample</b></td>
					<td width=\"50\%\" align='center'><b>$mapping_rate</b></td>
					</tr>
				";
				print "
					<tr>
					<td width=\"50\%\" align='center'><b>$sample</b></td>
					<td width=\"50\%\" align='center'><b>$mapping_rate</b></td>
					</tr>\n
				";
			
		}
	}
	elsif ($line=~/<!--FOR_Database_Web-->/)
	{
		if (!defined($hashinfor{'DBused'})){die "DBused info. not provided in Report_Information.txt\n\n";}
		$line=~s/<!--FOR_Database_Web-->/$hashinfor{'DBused'}/;
		print $line;
		print Align2 $line;
	}

	else{
		print Align2 $line;
	}
}
close(Align2);
#########################

print "Reading $localEreportFolder/Content/DE.html\n";


open(DE_Gene, "<","$localEreportFolder/Content/DE.html") or die;
my @DE=<DE_Gene>;close DE_Gene;
open(DE_Gene2, ">","$localEreportFolder/Content/DE.html");
foreach my $html_de(@DE){
	if ($html_de=~/<!----------Set_Write_Here----------->/) {
		
		#open(Design, "Exp_Design.txt") or die;
		#my @Design=<Design>;close Design;
		#$groupCompareInfo{$T."\t".$C}=$T."_vs_".$C;
		foreach my $Set (sort values %groupCompareInfo){
		#foreach my $Set(@Design){
			#$Set=~s/\t/_/g ;
			my $DEstr=$Set.$outStr;
			if (-e "$localEreportFolder/Content/DE_Gene/${DEstr}_detail.xlsx" && -e "$localEreportFolder/Content/DE_TS/TS_${DEstr}_detail.xlsx")
			{$DEstr.="_detail";}
			
			$html_de.="
					<tr>
					<td align='center'>$Set</td>						
					<td align='center'><a href=\"DE_Gene/$DEstr.xlsx\" download><img src=\"../assets/img/icon/XLSX.png\" alt=\"\" width=\"60\" height=\"60\" border=\"0\"></a></td>
					<td align='center'><a href=\"DE_TS/TS_$DEstr.xlsx\" download><img src=\"../assets/img/icon/XLSX.png\" alt=\"\" width=\"60\" height=\"60\" border=\"0\"></a></td>
					<td align='center'>
					<a href=\"Enrich_$DEstr.xlsx\" download><img src=\"../assets/img/icon/XLSX.png\" alt=\"\" width=\"60\" height=\"60\" border=\"0\"></a>&nbsp;&nbsp;
					<a href=\"Enrich_${DEstr}_GSEA.xlsx\" target=\"_blank\"><img src=\"../assets/img/icon/XLSX.png\" alt=\"\" width=\"60\" height=\"60\" border=\"0\"></a>
					</td>
					<td align='center'><a href=\"EnrichPlots_$DEstr.html\" target=\"_blank\"><img src=\"../assets/img/icon/html.png\" alt=\"\" width=\"60\" height=\"60\" border=\"0\"></a>&nbsp;&nbsp;
<a href=\"Enrich_$DEstr.html\" target=\"_blank\"><img src=\"../assets/img/icon/html.png\" alt=\"\" width=\"60\" height=\"60\" border=\"0\"></a>
</td>
			</tr>
				";
				
		}
		$html_de=~s/<!----------Set_Write_Here----------->//;
		print DE_Gene2 $html_de;
		print "$html_de\n";	
		
	}elsif ($html_de=~/<td colspan=\"7\">By default/){
		#print "<td colspan=\"7\">By default, we set p-value cut off to $Pvalue and abs(log${FC} FC) cut off to ".log2(${FC})." for discovering differential genes.<br>";
		print DE_Gene2 "<td colspan=\"7\">By default, we set p-value cut off to $Pvalue and abs(log${FC}FC) cut off to".log2(${FC})."for discovering differential genes.<br>";
	}elsif ($html_de=~/<!----------Set_pFC_Here----------->/){
		#print "<li>Genes with <em>p<\/em> value ≤ $Pvalue and ≥ ${FC}-fold changes were  considered significantly differentially expressed. <br>";
		print DE_Gene2 "<li>Genes with <em>p</em> value < $Pvalue and > ${FC}-fold changes were  considered significantly differentially expressed.<br>";
	}else{
		print DE_Gene2 $html_de;
	}
}
close(DE_Gene2);

print "Done!!\n\n";


exit;
#================

sub log2
{
	my $n= shift;
	return sprintf("%.2f",log($n)/log(2));
}



sub getCurrentTimeStr
{
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    my $nice_timestamp = sprintf ( "%04d%02d%02d%02d%02d%02d",
                                   $year+1900,$mon+1,$mday,$hour,$min,$sec);
    return $nice_timestamp;
}
sub getCurrentDateTime 
{
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    my $nice_timestamp = sprintf ( "%04d-%02d-%02d %02d:%02d:%02d",
                                   $year+1900,$mon+1,$mday,$hour,$min,$sec);
    return "[$nice_timestamp]";
}
sub trim
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

sub naturalSortInner {
	no strict;
	no warnings;
     $x = uc( shift );
    $y = uc( shift );
    if( !($x =~ /\d+(\.\d+)?/) ) {
        return $x cmp $y;
    }
    $xBefore = $`;
    $xMatch = $&;
    $xAfter = $';
    if( !($y =~ /\d+(\.\d+)?/) ) {
        return $x cmp $y;
    }
    if( $xBefore eq $` ) {
        if( $xMatch == $& ) {
            return naturalSortInner( $xAfter, $' );
        } else {
            return $xMatch <=> $&;
        }
    } else {
        return $x cmp $y;
    }
    print "\n<before: '$xBefore', match: '$xMatch', after: '$xAfter'>\n";
}
sub naturalSort {
	no strict;
	no warnings;
    naturalSortInner( $a, $b );
}
sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text;
}
