sub setting {	
	
	($Strandnesss,$Species)= @_;
	my %hash_set_strand;
	

	######## Strand ########
	
	#my $Strand_opt; #for cuffdiff
	#my $hisat2StrandPE;
	#my $hisat2StrandSE;
	if ($Strandnesss==1) {
		$hash_set_strand{"Strand_opt"}="fr-firststrand"; #Agilent SureSelect Stranded RNA library kit 
		$hash_set_strand{"hisat2StrandPE"}="RF";
		$hash_set_strand{"hisat2StrandSE"}="R";
		$hash_set_strand{"stringtieStrand"}="--rf";
		
	}elsif($Strandnesss==2){
		$hash_set_strand{"Strand_opt"}="fr-secondstrand";
		$hash_set_strand{"hisat2StrandPE"}="FR";
		$hash_set_strand{"hisat2StrandSE"}="F";
	}elsif($Strandnesss==3){
		$hash_set_strand{"hisat2StrandSE"}="fr-unstranded"; #for NexteraXT library kit
		$hash_set_strand{"hisat2StrandSE"}="";
		$hash_set_strand{"hisat2StrandSE"}="";
	}
	
	####### Path #######
	my $hostname=`hostname`; 
	chomp $hostname;
	my %hash_set_path;
	print "$hostname\n";
	
	if ($hostname eq "SuperMicro")
	{	
		$hash_set_path{"HISAT2indexPATH"}="/NGS/Database/HISAT2/Ensembl_r89";
		$hash_set_path{"STARfusionIndexPATH"}="/NGS/Database/STAR-Fusion";
		$hash_set_path{"GTFPATH"}="/NGS/Database/Annotation/GTF/Release_89";
		$hash_set_path{"BiomartAnnPATH"}="/home/shihjielai/RNAseq_PIP/Biomart_ann"; #/mouse_89_biomart.txt";
		$hash_set_path{"hisat2PATH"}="/software/hisat2-2.1.0";
		$hash_set_path{"fastqcPATH"}="/software/FastQC";
		$hash_set_path{"fastqcrevPATH"}="/software/Tools";
		$hash_set_path{"samtoolsPATH"}="/software/samtools-1.4";
		$hash_set_path{"stringtiePATH"}="/software/stringtie-1.3.3b.Linux_x86_64"; #1.3.4d
		$hash_set_path{"misoPATH"}="/miniconda3/bin";
		$hash_set_path{"starFusionPATH"}="/opt/miniconda3/envs/RNAseq/bin/STAR-Fusion";
		$hash_set_path{"arriba"}="";		
		$hash_set_path{"STAR"}="";
		#$hash_set_path{"gffcomparePATH"}="/usr/bin";
		$hash_set_path{"cuffdiffPATH"}="/usr/local/bin";
		$hash_set_path{"reformatPATH"}="/usr/local/genome/bbmap";
		
	}
	elsif ($hostname eq "R930")
	{
		$hash_set_path{"HISAT2indexPATH"}="/export/md1/databases/HISAT2_index/Ensembl_r101";
		$hash_set_path{"STARfusionIndexPATH"}="/export/md1/databases/_RNAseq_related";
		$hash_set_path{"GTFPATH"}="/export/md1/databases/GTF_Ann/Release_101";
		#$hash_set_path{"BiomartAnnPATH"}="/home/shihjielai/RNAseq_PIP/Biomart_ann/101"; #/mouse_89_biomart.txt";
		$hash_set_path{"BiomartAnnPATH"}="/home/shihjielai/RNAseq_PIP/Biomart_ann/101";
		$hash_set_path{"hisat2PATH"}="/opt/miniconda3/envs/RNAseq2/bin";
		$hash_set_path{"fastqcPATH"}="/usr/local/genome/FastQC";
		$hash_set_path{"fastqcrevPATH"}="/usr/local/genome/bin";
		$hash_set_path{"samtoolsPATH"}="/opt/miniconda3/envs/RNAseq2/bin";
		$hash_set_path{"stringtiePATH"}="/usr/local/genome/_RNAseq/stringtie-2.1.4.Linux_x86_64"; #1.3.4d
		$hash_set_path{"misoPATH"}="/usr/local/bin";
		$hash_set_path{"starFusionPATH"}="/opt/miniconda3/envs/RNAseq2/lib/STAR-Fusion";
		$hash_set_path{"arriba"}="/home/shihjielai/RNAseq_PIP/arriba_v2.0.0";		
		$hash_set_path{"STAR"}="/opt/miniconda3/envs/RNAseq2/bin";
		#$hash_set_path{"gffcomparePATH"}="/usr/local/genome/_RNAseq/gffcompare";
		$hash_set_path{"cuffdiffPATH"}="/usr/local/genome/_RNAseq/cufflinks-2.2.1.Linux_x86_64";
		$hash_set_path{"reformatPATH"}="/usr/local/genome/bbmap";
		$hash_set_path{"rmats"}="/export/md2/_RNAcases/A20190910/RNGS1091190_fusion/rmats_turbo_v4_1_1";
	}
	
	######## index & fa ########
	my %hash_set_genome;
	
	my %species2index=(
	1 => 'Homo_sapiens.GRCh38.dna.primary_assembly',
	2 => 'Mus_musculus.GRCm38.dna.primary_assembly',
	3 => 'Rattus_norvegicus.Rnor_6.0.dna.toplevel',
	4 => 'Arabidopsis_thaliana.TAIR10.dna.toplevel',
	5 => 'Caenorhabditis_elegans.WBcel235.dna.toplevel',
	6 => 'Danio_rerio.GRCz11.dna.primary_assembly',
	7 => 'Drosophila_melanogaster.BDGP6.28.dna.toplevel',
	8 => 'Oryza_sativa.IRGSP-1.0.dna.toplevel',
	9 => 'Cust'
	);
	
	$hash_set_genome{"HISAT2_index"}="$hash_set_path{\"HISAT2indexPATH\"}/$species2index{$Species}";
	$hash_set_genome{"FAREF"}="$hash_set_path{\"HISAT2indexPATH\"}/$species2index{$Species}.fa";
	
	######## gtf ########
	my %species2gtf=(
	1 => 'Homo_sapiens.GRCh38.101',
	2 => 'Mus_musculus.GRCm38.101',
	3 => 'Rattus_norvegicus.Rnor_6.0.101',
	4 => 'Arabidopsis_thaliana.TAIR10.48',
	5 => 'Caenorhabditis_elegans.WBcel235.101',
	6 => 'Danio_rerio.GRCz11.101',
	7 => 'Drosophila_melanogaster.BDGP6.28.101',
	8 => 'Oryza_sativa.IRGSP-1.0.48',
	9 => 'Cust'
	);
	$hash_set_genome{"GFFGTF"}="$hash_set_path{\"GTFPATH\"}/$species2gtf{$Species}.gtf";
	
	######## biomart ########
	my %species2biomart=(
	1 => 'Homo_sapiens.GRCh38.p13_biomart.txt',
	2 => 'Mus_musculus.GRCm38.p6_biomart.txt',
	#2 => 'mouse_89_biomart.txt',
	3 => 'Rattus_norvegicus.Rnor6.0_biomart.txt',
	4 => 'Arabidopsis_thaliana.TAIR10_biomart.txt',
	5 => 'Caenorhabditis_elegans.WBcel235_biomart.txt',
	6 => 'Danio_rerio.GRCz11_biomart.txt',
	#6 => 'dre_89_biomart.txt',
	7 => 'Drosophila_melanogaster.BDGP6.28_biomart.txt',
	8 => 'Oryza_sativa.Japonica_Group.IRGSP-1.0_biomart.txt',
	9 => 'Cust'
	);
	
	$hash_set_genome{"Biomart_annotation"}="$hash_set_path{\"BiomartAnnPATH\"}/$species2biomart{$Species}";
	
	######## Name ########
	# my %SpeciesName=(
	# 1 => 'Homo sapiens',
	# 2 => 'Mus musculus',
	# 3 => 'Rattus norvegicus',
	# 4 => 'Arabidopsis thaliana',
	# 5 => 'Caenorhabditis elegans',
	# 6 => 'Danio rerio',
	# 7 => 'Drosophila melanogaster',
	# 8 => 'Oryza sativa',
	# 9 => 'Custom'
	# );
	
	######## DBweb ########
	my %ensemblDBname=(
	1 => 'http://www.ensembl.org/Homo_sapiens/Info/Index',
	2 => 'http://nov2020.archive.ensembl.org/Mus_musculus/Info/Index',
	3 => 'http://www.ensembl.org/Rattus_norvegicus/Info/Index',
	4 => 'http://plants.ensembl.org/Arabidopsis_thaliana/Info/Index',
	5 => 'http://www.ensembl.org/Caenorhabditis_elegans/Info/Index',
	6 => 'http://www.ensembl.org/Danio_rerio/Info/Index',
	7 => 'http://www.ensembl.org/Drosophila_melanogaster/Info/Index',
	8 => 'http://plants.ensembl.org/Oryza_sativa/Info/Index',
	9 => 'Cust'
	);
	$hash_set_genome{"ensemblDBname"}="$ensemblDBname{$Species}";
	######## GF index ########
	my %starFusionDBindex=(
	1 => 'GRCh38_v27_CTAT_lib_Feb092018/ctat_genome_lib_build_dir.star2.7.5c',
	2 => '',
	3 => '',
	4 => '',
	5 => '',
	6 => '',
	7 => '',
	8 => '',
	9 => 'Cust'
	);
	$hash_set_genome{"starFusionDBindex"}="$hash_set_path{\"STARfusionIndexPATH\"}/$starFusionDBindex{$Species}";
	
	##human##
	if($Species==1){
		$hash_set_genome{"KEGG_symbol"}="hsa";
		$hash_set_genome{"OrgDB"}="org.Hs.eg.db";
		#$hash_set_genome{"misoGFFindex"}="$hash_set_path{\"GTFPATH\"}/Homo_sapiens.GRCh38.89";
		$hash_set_genome{"misoGFFindex"}="/export/md2/_RNAcases/CRO_TU/Homo_sapiens.GRCh38.89";
	##Mouse##	
	}elsif($Species==2){
		$hash_set_genome{"KEGG_symbol"}="mmu";
		$hash_set_genome{"OrgDB"}="org.Mm.eg.db";
	##Rat##
	}elsif($Species==3){
		$hash_set_genome{"KEGG_symbol"}="rno";
		$hash_set_genome{"OrgDB"}="org.Rn.eg.db";
	##Arabidopsis##
	}elsif($Species==4){
		$hash_set_genome{"KEGG_symbol"}="ath";
		$hash_set_genome{"OrgDB"}="org.At.tair.db";
	##C.elgean##
	}elsif($Species==5){
		$hash_set_genome{"KEGG_symbol"}="cel";
		$hash_set_genome{"OrgDB"}="org.Ce.eg.db";
	##ZebraFish##
	}elsif($Species==6){
		$hash_set_genome{"KEGG_symbol"}="dre";
		$hash_set_genome{"OrgDB"}="org.Dr.eg.db";
	##FruitFly##
	}elsif($Species==7){
		$hash_set_genome{"KEGG_symbol"}="dme";
		$hash_set_genome{"OrgDB"}="org.Dm.eg.db";
	##OSA##
	}elsif($Species==8){
		$hash_set_genome{"KEGG_symbol"}="osa";
		$hash_set_genome{"OrgDB"}="MeSH.Osa.eg.db";
	}elsif($Species==9){
		($hash_set_genome{"FAREF"},
		$hash_set_genome{"HISAT2_index"},
		$hash_set_genome{"GFFGTF"},
		$hash_set_genome{"Biomart_annotation"},
		$hash_set_genome{"KEGG_symbol"},
		$hash_set_genome{"OrgDB"})=cust("$hash_set_path{hisat2PATH}");
		#($FAREF,$HISAT2_index,$GFFGTF,$Biomart_annotation,$KEGG_symbol,$OrgDB);
	}
	
	return(\%hash_set_strand,\%hash_set_path,\%hash_set_genome);
}

sub cust{
	
	my ($hisat2PATH)= @_;
	my @fa=<*.fasta *.fa *fna>;  # custom reference seq
	my $FAREF;
	if (-e $fa[0]){$FAREF=$fa[0];}
	else
	{
		print "No reference fasta/fa file found!!\n\n"; 
		die "No reference fasta/fa file found!!\n\n";
	}
	my $FAREFprefix=$FAREF;  $FAREFprefix=~s/(\.(fa|fasta|fas|fna))$//;
	if (!-e "$FAREFprefix.1.ht2")
	{
		print "CUST: $FAREF\n";
		my $cmdStr="$hisat2PATH/hisat2-build $FAREF $FAREFprefix";
		print "$cmdStr\n";
		system("$cmdStr");
	}
	my $HISAT2_index="./$FAREFprefix";
	my $KEGG_symbol="Cust";
	my $OrgDB="org.MeSH.Cust.db";
	my @ann=<*_ann.txt>;
	my $Biomart_annotation;
	if (scalar(@ann))
	{
		$Biomart_annotation=$ann[0];
	}
	my @gtf=<*.gtf *.gff>;
	my $GFFGTF;
	if (scalar(@gtf))
	{
		$GFFGTF=$gtf[0];
		print " Annotation GTF/GFF file found: $GFFGTF\n";  
	}
	else
	{ 
		die "No Annotation GTF/GFF file found!!\n\n";
	}

	my @faa=<*.faa1>;
	if (scalar(@fan))
	{
		system("bash /usr/local/genome/_RNAseq/interproscan-5.36-75.0/interproscan.sh -goterms -i $faa[0] -f tsv -appl Pfam,PANTHER,Hamap");
		my @tsv=<*.tsv>;
		open(CP, "$tsv[0]") or die "NO tsv";
		open(ANS, ">$tsv[0].txt");
		my @CP=<CP>; close CP;
		my %pname;
		foreach (@CP){
			chomp $_;
			my @Temp=split /\t/, $_;
			if (exists $pname{$Temp[0]}){
			}else{$pname{$Temp[0]}="";}
			foreach my $a(@Temp){
				if($Temp[3]=~/^(Pfam|PANTHER|Hamap)/){
					if ($a=~/GO:/){
						my @Temp1=split /\|/, $a;
						foreach $b (@Temp1){
							if($pname{$Temp[0]}=~/$b/){
							}else{
								$pname{$Temp[0]}=$pname{$Temp[0]}.",$b";
							}
						}				
					}					
				}
		
			}
		}
		foreach my $key (keys %pname){
			$pname{$key}=~s/,//;
			print ANS "$key\t$pname{$key}\n";
		}
		close ANS;
	}

	return ($FAREF,$HISAT2_index,$GFFGTF,$Biomart_annotation,$KEGG_symbol,$OrgDB);
}

sub getLIMSinfo{
	use DBI;
	binmode(STDIN, ':encoding(utf8)');
	binmode(STDOUT, ':encoding(utf8)');
	binmode(STDERR, ':encoding(utf8)'); 
	my ($localPATH,$DBused,$POinfoCheck)= @_;
	$POnum=$localPATH;
	my %POrecord;
	my %sequencing;
	my %examination;
	#Connect to LIMS database
	my $dbh= &connectLIMS();
	my $commonSampleName;
	#############################################################################
	my $SQLstr="select * from `POrecord` where `poNUM`='$POnum'";
	my $sth = $dbh->prepare($SQLstr);
	$sth->execute or die "Unable to execute query: $SQLstr\n *$dbh->errstr\n"; 
	my $rows = $sth->rows;
	while(my $ref = $sth->fetchrow_hashref()) 
	{
		$POrecord{$POnum}{'POid'}=$ref->{'id'};
		$POrecord{$POnum}{'PJid'}=$ref->{'project_id'};
		$POrecord{$POnum}{'POstatus'}=$ref->{'status'};
	}
	print "LIMS info===============================\n";
	if ($rows==0) # this POnum not in LIMS database => suggest it is a fake number
	{
		$POrecord{$POnum}{'POid'}="NA";
		$POrecord{$POnum}{'PJid'}="NA";
		$POrecord{$POnum}{'POstatus'}="NA";
		$POrecord{$POnum}{'Customer_id'}="NA";
		$POrecord{$POnum}{'projectNum'}="NA";
		$POrecord{$POnum}{'PJtype'}="NA";
		$POrecord{$POnum}{'PJclass'}="NA";
		$POrecord{$POnum}{'PJservice'}="NA";
		$POrecord{$POnum}{'PJspecies'}="NA";
		$POrecord{$POnum}{'PJnote'}="internalTEST?";
		$POrecord{$POnum}{'Customer_Name'}="NA";
		$POrecord{$POnum}{'Customer_Institute'}="NA";
		$POrecord{$POnum}{'sampleNumSeq'}=0;
		$POrecord{$POnum}{'sampleNumExm'}=0;
		$POrecord{$POnum}{'NewStatus'}="";
		print "!! No data found for POnum='$POnum' !!\n\n";
	}
	else{
		$SQLstr="select * from `projects` where `id`=$POrecord{$POnum}{'PJid'}";
		#print "$SQLstr\n";
		$sth = $dbh->prepare($SQLstr);
		$sth->execute or die "Unable to execute query: $SQLstr\n *$dbh->errstr\n"; 
		$rows = $sth->rows;
		while(my $ref = $sth->fetchrow_hashref()){
			$POrecord{$POnum}{'Customer_id'}=$ref->{'Customer_id'};
			$POrecord{$POnum}{'projectNum'}=$ref->{'projectNum'};
			$POrecord{$POnum}{'PJtype'}=$ref->{'PJtype'};
			$POrecord{$POnum}{'PJclass'}=$ref->{'PJclass'};
			$POrecord{$POnum}{'PJservice'}=$ref->{'PJservice'} || "";
			$POrecord{$POnum}{'PJspecies'}=$ref->{'PJspecies'} || "";
			$POrecord{$POnum}{'PJnote'}=$ref->{'note'} || "";
			$SQLstr="SELECT * FROM `customers` WHERE `id` = $POrecord{$POnum}{'Customer_id'}";
			#print "$SQLstr\n";
			my $tth = $dbh->prepare($SQLstr);
			$tth->execute or die "Unable to execute query: $SQLstr\n *$dbh->errstr\n"; 
			$rows = $tth->rows;
			while(my $ref = $tth->fetchrow_hashref()) 
			{
				$POrecord{$POnum}{'Customer_Name'}=$ref->{'c_name'};
				$POrecord{$POnum}{'Customer_Institute'}=$ref->{'institute'};
			}
		}
		print "$POrecord{$POnum}{'projectNum'}|$POnum|$POrecord{$POnum}{'Customer_Name'}|$POrecord{$POnum}{'Customer_Institute'}\n";
		print "$POrecord{$POnum}{'PJtype'}|$POrecord{$POnum}{'PJclass'}|$POrecord{$POnum}{'PJservice'}|$POrecord{$POnum}{'PJspecies'}\n";
		
		
		
		$SQLstr="select * from `sequencing` where `POrecord_id`='$POrecord{$POnum}{'POid'}'";
		#print "$SQLstr\n";
		$sth = $dbh->prepare($SQLstr);
		$sth->execute or die "Unable to execute query: $SQLstr\n *$dbh->errstr\n"; 
		#$rows = $sth->rows;
		$POrecord{$POnum}{'sampleNumSeq'}=$sth->rows;
		while(my $ref = $sth->fetchrow_hashref()) {
			$commonSampleName=$ref->{'sampleName'}; 
			$sequencing{$POnum}{$commonSampleName}{'id'}=$ref->{'id'};
			$sequencing{$POnum}{$commonSampleName}{'length_bp'}=$ref->{'length_bp'};
			$sequencing{$POnum}{$commonSampleName}{'type'}=$ref->{'type'} ;  #PE or SE
			$sequencing{$POnum}{$commonSampleName}{'targetQuantity'}=$ref->{'targetQuantity'};
			$sequencing{$POnum}{$commonSampleName}{'quantityUnit'}=$ref->{'quantityUnit'};
			$sequencing{$POnum}{$commonSampleName}{'targetQuantityBp'}=$ref->{'targetQuantityBp'};
			
			print "[上機定序]\t$commonSampleName\t$sequencing{$POnum}{$commonSampleName}{'targetQuantity'}$sequencing{$POnum}{$commonSampleName}{'quantityUnit'}\t$sequencing{$POnum}{$commonSampleName}{'targetQuantityBp'} bases\t$sequencing{$POnum}{$commonSampleName}{'length_bp'}$sequencing{$POnum}{$commonSampleName}{'type'}\n";
		
		}
		$SQLstr="select * from `library` where `POrecord_id`='$POrecord{$POnum}{'POid'}'";
		#print "$SQLstr\n";
		$sth = $dbh->prepare($SQLstr);
		$sth->execute or die "Unable to execute query: $SQLstr\n *$dbh->errstr\n"; 
		#$rows = $sth->rows;
		$POrecord{$POnum}{'sampleNumSeq'}=$sth->rows;
		while(my $ref = $sth->fetchrow_hashref()) {
	
			$sequencing{$POnum}{$commonSampleName}{'libraryKit'}=$ref->{'libraryKit'};

			print "[library]\t$sequencing{$POnum}{$commonSampleName}{'libraryKit'}\n";
		}
	}		
	my $reportInfoStr=
		"Organism	$POrecord{$POnum}{'PJspecies'}\n".
		"PONumber	$POnum\n".
		"Institute	$POrecord{$POnum}{'Customer_Institute'}\n".
		"Customer	$POrecord{$POnum}{'Customer_Name'}\n".
		"Tel	-\n".
		"Spec	$sequencing{$POnum}{$commonSampleName}{'length_bp'}$sequencing{$POnum}{$commonSampleName}{'type'}\n".
		"Amount	$sequencing{$POnum}{$commonSampleName}{'targetQuantity'}$sequencing{$POnum}{$commonSampleName}{'quantityUnit'}\n".
		"kit	$sequencing{$POnum}{$commonSampleName}{'libraryKit'}\n";
			
	
	if (defined($DBused)){
		$reportInfoStr.="DBused	$DBused\n";
			
		my $reportInfo="${POnum}_reportInfo.txt";
		if (!-e $reportInfo){
			print "Writing '$reportInfo'...\n";
			open ReportInfo,">","$reportInfo" || die "Cannot open $reportInfo!!\n\n";
			binmode(ReportInfo, ':encoding(utf8)');
			print ReportInfo $reportInfoStr;
			close (ReportInfo);
		}
		else{
			print "Found '$reportInfo'...Skip!\n";
		}
	}
		print "========================================\n".$reportInfoStr."========================================\n\n";	
	if (defined($POinfoCheck)){exit;}
}

sub connectLIMS
{
	#Connect to LIMS database
	my $DBhost="192.168.5.114"; #`DBhost`;
	my ($username,$password);
	chomp $DBhost;

		$username="ngs";
		$password='bio@NGS80158777';
	
	my $dbName="LIMS";
	my $dbh;
	$dbh = DBI->connect("DBI:mysql:$dbName:$DBhost",$username, $password, { mysql_enable_utf8 => 1,}) or die "Unable to connect to contacts Database: $dbh->errstr\n"; 
	$dbh->{RaiseError} = 1; 
	$dbh->{mysql_auto_reconnect} = 1; 

	my $SQLstr="SET NAMES 'utf8'";
	my $sth = $dbh->prepare($SQLstr);
	$sth->execute or die "Unable to execute query: $SQLstr\n *$dbh->errstr\n";
	#############################################################################
	return $dbh;

}

sub gdupload
{
	my ($localPATH,$POID)= @_;
	open(RSCRIPT,">googleupload.sh") or die;
	$localPATH=~s/\/Run*+\n$//g;
	print RSCRIPT "
		#!/bin/bash
		if [ ! -d $localPATH/data/ ];then
			mkdir $localPATH/Bam
			mv *.bam *.bai $localPATH/Bam
			mkdir $localPATH/data
			mv $localPATH/*rawFASTQ*/ $localPATH/*trimUMIFASTQ*/ $localPATH/*trimmedFASTQ*/ $localPATH/*Bam*/ $localPATH/data
			ln -s $localPATH/data/Bam/*.bam ./
			rm -rf *.gz
			ln -s $localPATH/data/*trimmedFASTQ*/*.gz ./
		fi
		if [ -d $localPATH/data/rawFASTQ/ ];then
			mkid=\$(echo \$(gdrive mkdir -p 0B-I5dOIaBUFJfnB2b0lqUXQxYXFBVlhTRFliVkYtYS1HVl84aEZnM043MHVHOUdGYlR1ZUE $POID) | awk '{print \$2}')
			echo \$mkid
			gdrive upload --recursive -p \$mkid \"$localPATH/data\"&
		fi
		";
	close RSCRIPT;
	if(! -d "$localPATH/data"){
		system("bash googleupload.sh > gd.log");
	}
}




return 1;
