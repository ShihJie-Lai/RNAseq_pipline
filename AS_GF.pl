require "/home/shihjielai/RNAseq_PIP/excelO.pl";
use Cwd;

sub STAR_GF{
	my ($starFusionRef,$starFusionPATH,$Paired_Check,%sample2fastq)=@_;#,@Sample_Name
	if (!-e "$starFusionRef"){die "STAR-Fusion Not Supported for Species other then Right Now!!\n\n";}
	system("mkdir -p ../4_fusion");
	system("mkdir -p ../6_Report/Gene_Fusion");
	my $cmdStr;
	my @Sample_Name= keys %sample2fastq;
	for ( my $i=0;$i<=$#Sample_Name;$i++)
	{
		my $sampleName=$Sample_Name[$i];
		
		if (!-e "../4_fusion/${sampleName}_geneFusion_outdir/star-fusion.fusion_predictions.tsv" || !-e "../4_fusion/${sampleName}_geneFusion_outdir/star-fusion.fusion_predictions.abridged.tsv")
		{
			if ($Paired_Check==0) #SE
			{
				$cmdStr="$starFusionPATH/STAR-Fusion --genome_lib_dir $starFusionRef --left_fq $sample2fastq{$sampleName}{R1} --CPU 12 --output_dir ../4_fusion/${sampleName}_geneFusion_outdir &> ../4_fusion/${sampleName}_starFusion.log";
			}
			else
			{
				$cmdStr="$starFusionPATH/STAR-Fusion --genome_lib_dir $starFusionRef --left_fq $sample2fastq{$sampleName}{R1} --right_fq $sample2fastq{$sampleName}{R2} --CPU 12 --output_dir ../4_fusion/${sampleName}_geneFusion_outdir &> ../4_fusion/${sampleName}_starFusion.log";
			}
			print " $cmdStr\n";
			system("$cmdStr");
		}
		
		if (!-e "../6_Report/Gene_Fusion/${sampleName}.fusion_predictions.tsv")
		{
			$cmdStr="cp ../4_fusion/${sampleName}_geneFusion_outdir/star-fusion.fusion_predictions.tsv ../6_Report/Gene_Fusion/${sampleName}.fusion_predictions.tsv";
			print "$cmdStr\n";
			system("$cmdStr");
			if (!-e "../6_Report/Gene_Fusion/${sampleName}.fusion_predictions.xlsx")
			{
				&tabtoExcel("../6_Report/Gene_Fusion/${sampleName}.fusion_predictions.tsv","../6_Report/Gene_Fusion/${sampleName}.fusion_predictions.xlsx");	
			}
		}
		if (!-e "../6_Report/Gene_Fusion/${sampleName}.fusion_predictions.abridged.tsv")
		{
			$cmdStr="cp ../4_fusion/${sampleName}_geneFusion_outdir/star-fusion.fusion_predictions.abridged.tsv ../6_Report/Gene_Fusion/${sampleName}.fusion_predictions.abridged.tsv";
			print "$cmdStr\n";
			print RUN "$cmdStr\n";
			system("$cmdStr");
			if (!-e "../6_Report/Gene_Fusion/${sampleName}.fusion_predictions.abridged.xlsx")
			{
				&tabtoExcel("../6_Report/Gene_Fusion/${sampleName}.fusion_predictions.abridged.tsv","../6_Report/Gene_Fusion/${sampleName}.fusion_predictions.abridged.xlsx");	
			}
		}
	}
}

sub miso_AS{

	my ($misoPATH,$misoGFFindex)=@_;
	if (!-d $misoGFFindex){die "misoGFFindex: $misoGFFindex not found!!\n\n";}
	system("mkdir -p ../5_misoAS");
	system("mkdir -p ../6_Report/Alternative_Splicing");
	my $cmdStr;
	my @Sample_Name=<*.sorted.bam>;
	#my @Sample_Name= keys %sample2fastq;
	for ( my $i=0;$i<=$#Sample_Name;$i++)
	{
		my $sampleName=$Sample_Name[$i];
		my @Temp2=split /\./, $sampleName;
		$sampleName=$Temp2[0];
		my $BAM="../5_misoAS/${sampleName}.sorted.filtered.bam";
		if (-e "${sampleName}.sorted.bam" && !-e $BAM)
		{
			$cmdStr="samtools view -h ${sampleName}.sorted.bam | awk 'length(\$10) > 149 || \$1 ~ /^\@/' | samtools view -bS - > $BAM";
			print " $cmdStr\n";
			system("$cmdStr");
		}
		if (-e $BAM && !-e "$BAM.bai")
		{
			$cmdStr="samtools index $BAM";
			print " $cmdStr\n";
			system("$cmdStr");
		}
		
		if (!-e "../5_misoAS/${sampleName}_MISO_summary/summary/${sampleName}_MISO.miso_summary")
		{
			$cmdStr="$misoPATH/miso --run $misoGFFindex $BAM --output-dir ../5_misoAS/${sampleName}_MISO --read-len 150 --paired-end 225 150 -p 10  > ../5_misoAS/${sampleName}_MISO.log";			
			print " $cmdStr\n";
			system("$cmdStr");
			
			$cmdStr="$misoPATH/summarize_miso --summarize-samples ../5_misoAS/${sampleName}_MISO/ ../5_misoAS/${sampleName}_MISO_summary/  > ../5_misoAS/${sampleName}_MISO_summary.log";
			print " $cmdStr\n";
			system("$cmdStr");

		}
		
		if (!-e "../6_Report/Alternative_Splicing/${sampleName}.misoAS.tsv")
		{
			$cmdStr="cp -a ../5_misoAS/${sampleName}_MISO_summary/summary/${sampleName}_MISO.miso_summary ../6_Report/Alternative_Splicing/${sampleName}.misoAS.tsv";
			print "$cmdStr\n";
			system("$cmdStr");
		}
		if (-e "../6_Report/Alternative_Splicing/${sampleName}.misoAS.tsv" && !-e "../6_Report/Alternative_Splicing/${sampleName}.misoAS.xlsx")
		{
			&tabtoExcel("../6_Report/Alternative_Splicing/${sampleName}.misoAS.tsv","../6_Report/Alternative_Splicing/${sampleName}.misoAS.xlsx");
		}
	}
	
	
}

sub sep_Report{

	my $rawFolder="../1_rawFASTQ";
	my $trimmedFolder="../2_trimmedFASTQ";
	my $targetFolder="../6_Report";

	my @SampleName=();
	my $SampleName1;
	my %sample2fastq;
	my %fastq2sample;

	my($day, $month, $year) = (localtime)[3,4,5];
	$month = sprintf '%02d', $month+1;
	$day   = sprintf '%02d', $day;
	my $todayStr=($year-100).$month.$day;   #190314

	my %FileSize;
	my @FASTQ_R1=<$rawFolder/*_R1*.fastq*>;
	my @FASTQ_R2=<$rawFolder/*_R2*.fastq*>;
	my @FASTQ_R1t=();
	my @FASTQ_R2t=();
	
	for(my $i=0;$i<=$#FASTQ_R1;$i++)
	{	

		if (!-e $FASTQ_R2[$i]){die "$FASTQ_R2[$i] not found!!\n\n";}
		#A010002RSS00_S127_R1_001.fastq.gz
		if ($FASTQ_R1[$i]=~/$rawFolder\/(.+)_S[0-9]+(_L00[1-9]+)?_R1(_001)?\.fastq*/)
		{
			$SampleName1=$1;
		}
		else
		{
			die "$FASTQ_R1[$i]: unknown format!!\n\n";
		}
		$fastq2sample{$FASTQ_R1[$i]}=$SampleName1;
		$fastq2sample{$FASTQ_R2[$i]}=$SampleName1;
		
		push @SampleName,$SampleName1;
		
		$FileSize{$FASTQ_R1[$i]} = (stat $FASTQ_R1[$i])[7];
		$FileSize{$FASTQ_R2[$i]} = (stat $FASTQ_R2[$i])[7];
		
		print "$FASTQ_R1[$i] size = $FileSize{$FASTQ_R1[$i]}\n";
		print "$FASTQ_R2[$i] size = $FileSize{$FASTQ_R2[$i]}\n";

	}
	
	open(ckList, ">checkList_RSS.txt") or die "Cannot open checkList_RSS.txt!!\n\n";
	
	print "\n\nEfileCode\tType\tRawFileNum\tRawFileSize\tReportFileNum\tReportFileSize\tCheck\n";
	print ckList "\n\nEfileCode\tType\tRawFileNum\tRawFileSize\tReportFileNum\tReportFileSize\tCheck\n";
	for(my $i=0;$i<=$#FASTQ_R1;$i++)
	{		
		my $uploadFolderName=$SampleName[$i]."$todayStr";
		$FASTQ_R1t[$i]=$FASTQ_R1[$i]; $FASTQ_R1t[$i]=~s/$rawFolder/$targetFolder\/$uploadFolderName/;
		$FASTQ_R2t[$i]=$FASTQ_R2[$i]; $FASTQ_R2t[$i]=~s/$rawFolder/$targetFolder\/$uploadFolderName/;
		
		`mkdir -p $targetFolder/$uploadFolderName`;
		if (!-e $FASTQ_R1t[$i]) {print "copy $FASTQ_R1t[$i]...\n";`ln -s $FASTQ_R1[$i] $targetFolder/$uploadFolderName/.`;}
		if (!-e $FASTQ_R2t[$i]) {print "copy $FASTQ_R2t[$i]...\n";`ln -s $FASTQ_R2[$i] $targetFolder/$uploadFolderName/.`;}
		
		my $reportFileName=$SampleName[$i]."_Report.rar";
		if (-e "$targetFolder/$rarReportFile")
		{	
			print "copy $targetFolder/$rarReportFile ...\n";
			`cp -ar $targetFolder/$rarReportFile $targetFolder/$uploadFolderName/$reportFileName`;
		}
		elsif (-e "$rarReportFile")
		{
			print "copy $rarReportFile ...\n";
			`cp -ar $rarReportFile $targetFolder/$uploadFolderName/$reportFileName`;
		}
		
		$FileSize{$reportFileName} = (stat "$targetFolder/$uploadFolderName/$reportFileName")[7];
		my $FASTQsize=$FileSize{$FASTQ_R1[$i]} + $FileSize{$FASTQ_R2[$i]};
		my $checkRow="$uploadFolderName\tRSS\t2\t$FASTQsize\t1\t$FileSize{$reportFileName}";
		
		print "$checkRow\n";
		print ckList "$checkRow\n";
	}
	print "\n\n";
}


sub arriba_GF{
	my ($starFusionRef,$arribaPATH,%sample2fastq)=@_;
	if (!-e "$starFusionRef"){die "arriba Not Supported for Species other then Right Now!!\n\n";}
	system("mkdir -p ../7_arriba");
	my $cmdStr;
	my @Sample_Name= keys %sample2fastq;
	for ( my $i=0;$i<=$#Sample_Name;$i++)
	{
		my $sampleName=$Sample_Name[$i];
		
		if (!-e "../7_arriba/${sampleName}.fusions.tsv" )
		{
			$cmdStr="$arribaPATH/run_arriba.sh $starFusionRef/ref_genome.fa.star.idx/ $starFusionRef/ref_annot.gtf $starFusionRef/ref_genome.fa $arribaPATH/database/blacklist_hg38_GRCh38_v2.0.0.tsv $arribaPATH/database/known_fusions_hg38_GRCh38_v2.0.0.tsv $arribaPATH/database/protein_domains_hg38_GRCh38_v2.0.0.gff3 16 $sampleName $sample2fastq{$sampleName}{R1} $sample2fastq{$sampleName}{R2} >${sampleName}_arriba.log";
			print "$cmdStr\n";
			system("$cmdStr");
			system("mv $sampleName.*.tsv $sampleName.Aligned.sortedByCoord.out.* ../7_arriba");
			$cmdStr="$arribaPATH/draw_fusions.R --fusions=../7_arriba/${sampleName}.fusions.tsv --alignments=../7_arriba/${sampleName}.Aligned.sortedByCoord.out.bam --output=${sampleName}.fusions.pdf --annotation=$starFusionRef/ref_annot.gtf --cytobands=$arribaPATH/database/cytobands_hg38_GRCh38_v2.0.0.tsv --proteinDomains=$arribaPATH/database/protein_domains_hg38_GRCh38_v2.0.0.gff3";
			print "$cmdStr\n";
			system("$cmdStr");	
		}
		if (! -e "../7_arriba/${sampleName}_mapping.log"){
			&tabtoExcel("../7_arriba/${sampleName}.fusions.tsv","../7_arriba/${sampleName}.fusions.xlsx");	
			my $flagstatStr=`samtools flagstat ../7_arriba/${sampleName}.Aligned.sortedByCoord.out.bam`;
			my ($totalReadNum,$pairMateMappedNum,$singletonNum,$mappingRate)=();
			if ($flagstatStr=~/([0-9]+) \+ [0-9]* paired in sequencing/){$totalReadNum=$1;}
			if ($flagstatStr=~/([0-9]+) \+ [0-9]* with itself and mate mapped/){$pairMateMappedNum=$1;}
			if ($flagstatStr=~/([0-9]+) \+ [0-9]* singletons/){$singletonNum=$1;}
			$mappingRate=sprintf "%.2f%% overall alignment rate",($pairMateMappedNum+$singletonNum)/$totalReadNum *100;
			$mappingRate="$flagstatStr\n"."$mappingRate";
			open (ans2,">../7_arriba/${sampleName}_mapping.log");
			print ans2 "$mappingRate";
			close ans2;
		}

		if (-e "../7_arriba/${sampleName}.fusions.tsv"){
			open(IN2, "<", "../7_arriba/${sampleName}.fusions.tsv");
			open ( bamlist3, ">", "../7_arriba/${sampleName}.fusions_high.tsv");
			my @IN3=<IN2>;close (IN2);
			foreach my $ab(@IN3){
				chomp $ab;
				my @Temp2=split /\t/, $ab;
				if( $Temp2[0]=~/gene1/ ){
					print bamlist3 "$ab\n";
				}
				if( $Temp2[14]=~/high/ ){
					print bamlist3 "$ab\n";
				}
			}
			close (bamlist3);
			open(IN2, "<", "../7_arriba/${sampleName}.fusions_high.tsv");
			@IN3=<IN2>;close (IN2);
			if ($#IN3<5){
				system("cat ../7_arriba/${sampleName}.fusions.tsv | head -n 20 > ../7_arriba/${sampleName}.fusions_high.tsv");
			}
			my $cmdStr="$arribaPATH/draw_fusions.R --fusions=../7_arriba/${sampleName}.fusions_high.tsv --alignments=../7_arriba/${sampleName}.Aligned.sortedByCoord.out.bam --output../7_arriba/=${sampleName}.fusions.pdf --annotation=$starFusionRef/ref_annot.gtf --cytobands=$arribaPATH/database/cytobands_hg38_GRCh38_v2.0.0.tsv --proteinDomains=$arribaPATH/database/protein_domains_hg38_GRCh38_v2.0.0.gff3";
			print "$cmdStr\n";
			#system("$cmdStr");
		}
	}
}


sub GATK_snp{
	my ($sam,$starFusionRef)=@_;
	#system("ln -s $starFusionRef/ref_genome.fa ./");
	#my $bas=basename($starFusionRef/ref_genome.fa,'.fa');
	#system("$sam/samtools faidx ./$bas");
	#system("mkdir -p ../8_GATK");
	
	open(ANS,">GATK_snp.sh") or die;
	print ANS "
#!/bin/bash

ln -s $starFusionRef/ref_genome.fa ./
$sam/samtools faidx ./ref_genome.fa
mkdir -p ../8_GATK
java -jar /usr/local/genome/_toolSource/picard-tools-2.8.3/build/libs/picard.jar CreateSequenceDictionary REFERENCE=ref_genome.fa OUTPUT=ref_genome.dict

ls | grep _R1.fastq.gz | while read ids ;
do 
echo \$ids
aa=\$(echo \$ids |sed  \"s/_R1.fastq.gz//g\")
echo \$aa

STAR --genomeDir $starFusionRef/ref_genome.fa.star.idx/ --readFilesIn \${aa}_R1.fastq.gz \${aa}_R2.fastq.gz --runThreadN 20 --readFilesCommand zcat ######

genomeDir_2=../8_GATK/human38_\$aa
mkdir \$genomeDir_2

STAR --runMode genomeGenerate --genomeDir \$genomeDir_2 --genomeFastaFiles $starFusionRef/ref_genome.fa --sjdbFileChrStartEnd ./SJ.out.tab --sjdbOverhang 75 --runThreadN 20

STAR --genomeDir \$genomeDir_2 --readFilesIn \${aa}_R1.fastq.gz \${aa}_R2.fastq.gz --runThreadN 20 --readFilesCommand zcat

java -jar /usr/local/genome/_toolSource/picard-tools-2.8.3/build/libs/picard.jar AddOrReplaceReadGroups I=Aligned.out.sam O=../8_GATK/rg_\${aa}.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample

mv Aligned.out.sam ../8_GATK/\${aa}_Aligned.out.sam

java -jar /usr/local/genome/_toolSource/picard-tools-2.8.3/build/libs/picard.jar MarkDuplicates I=../8_GATK/rg_\${aa}.bam O=../8_GATK/dedupped_\${aa}\.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=../8_GATK/\${aa}.metrics

/usr/local/genome/gatk-4.1.1.0/gatk SplitNCigarReads --reference ref_genome.fa -I ../8_GATK/dedupped_\${aa}.bam -O ../8_GATK/split_\${aa}.bam

/usr/local/genome/gatk-4.1.1.0/gatk HaplotypeCaller -R ref_genome.fa -I ../8_GATK/split_\${aa}.bam --dont-use-soft-clipped-bases true -stand-call-conf 20.0 -O ../8_GATK/\${aa}.vcf

/usr/local/genome/gatk-4.1.1.0/gatk VariantFiltration -R ref_genome.fa -V ../8_GATK/\${aa}.vcf -window 35 -cluster 3 --filter-name FS -filter \"FS > 30.0\" --filter-name QD -filter \"QD < 2.0\" -O ../8_GATK/\${aa}_filt.vcf

done
";
close ANS;
system("bash GATK_snp.sh");

}

sub circ_RNA_CIRCexplorer2{
	my ($starFusionRef)=@_;	
	open(ANS,">circ_RNA.sh") or die;
	print ANS "
#!/bin/bash
if [ ! -d ../9_circRNA];then
	mkdir -p ../9_circRNA
fi
ls | grep _R1_fastq.gz | while read ids ;
do
	echo \$ids
	aa=\$(echo \$ids |sed \"s/_R1.fastq.gz//g\")
	if [ ! -d ../9_circRNA/\$aa];then
		echo \$aa
		STAR --chimSegmentMin 10 --runThreadN 10 --genomeDir $starFusionRef/ref_genome.fa.star.idx/ --readFilesIn \${aa}_R1.fastq.gz \${aa}_R2.fastq.gz --readFilesCommand zcat
		CIRCexplorer2 parse -t STAR Chimeric.out.junction -b ../9_circRNA/\$aa/\${aa}_junction.bed -f > CIRCexplorer2_parse.log
		mv Chimeric.out.junction ../9_circRNA/\$aa
		CIRCexplorer2 annotate -r GRCh38_ref.txt -g $starFusionRef/ref_genome.fa -b ../9_circRNA/\${aa}_junction.bed -o ../9_circRNA/\$aa/\${aa}_circularRNA_known.txt > CIRCexplorer2_annotate.log
		sed -i '1i\\chrom\\tstart\\tend\\tname\\tscore\\tstrand\\tthickStart\\tthickEnd\\titemRgb\\texonCount\\texonSizes\\texonOffsets\\treadNumber\\tcircType\\tgeneName\\tisoformName\\tindex\\tflankIntron' ../9_circRNA/\$aa/\${aa}_circularRNA_known.txt
	fi
done
";
	close ANS;
	system("bash circ_RNA.sh");
	my $dir = "../8_GATK";
	opendir DIR, $dir;
	my @etc_all_file = readdir DIR;
	closedir DIR;
	@etc_all_file = grep {!/\./} @etc_all_file;
	foreach my $bb (@etc_all_file){
		$aa="$dir/${bb}/${bb}_circularRNA_known.txt";
		&tabtoExcel("$aa","$dir/${bb}/${bb}_circularRNA_known.xlsx");
	}
}



sub STAR_GF_FFPE{
	my ($starFusionRef,$starFusionPATH,$STARPATH,$arribaPATH,%sample2fastq)=@_;#,@Sample_Name
	if (!-e "$starFusionRef"){die "STAR-Fusion Not Supported for Species other then Right Now!!\n\n";}
	system("mkdir -p ../10_FFPE_fusion");
	system("mkdir -p ../6_Report/Gene_Fusion_FFPE");
	my $cmdStr;
	my @Sample_Name= keys %sample2fastq;
	for ( my $i=0;$i<=$#Sample_Name;$i++)
	{
		my $sampleName=$Sample_Name[$i];
		
		if (!-e "../10_FFPE_fusion/${sampleName}_geneFusion_outdir/star-fusion.fusion_predictions.tsv" || !-e "../10_FFPE_fusion/${sampleName}_geneFusion_outdir/star-fusion.fusion_predictions.abridged.tsv")
		{	
			$cmdStr="$STARPATH/STAR --runThreadN 16 --genomeDir $starFusionRef/ref_genome.fa.star.idx/ --genomeLoad NoSharedMemory --readFilesIn ${sampleName}_R1.fastq.gz ${sampleName}_R2.fastq.gz --readFilesCommand zcat --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 --outFilterMultimapNmax 100 --outFilterMismatchNmax 33 --seedSearchStartLmax 12 --alignSJoverhangMin 15 --outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0.3 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 1 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --outReadsUnmapped None --chimOutJunctionFormat 1 --alignSJDBoverhangMin 10 --alignMatesGapMax 200000 --alignIntronMax 200000 --outSAMstrandField intronMotif --outSAMattrRGline ID:GRPundef --chimMultimapScoreRange 10 --chimNonchimScoreDropMin 10 --chimMultimapNmax 50 --outFilterType BySJout --alignSplicedMateMapLminOverLmate 0.5 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --twopassMode Basic  > ${sampleName}_fusion.bam";
			print " $cmdStr\n";
			system("$cmdStr");
			$cmdStr="samtools sort -@ 24 ${sampleName}_fusion.bam -o ${sampleName}_fusion_sort.bam";
			print " $cmdStr\n";
			system("$cmdStr");
			print " rm -rf ${sampleName}_fusion.bam\n";
			system ("rm -rf ${sampleName}_fusion.bam");
			$cmdStr="samtools index ${sampleName}_fusion_sort.bam";
			print " $cmdStr\n";
			system("$cmdStr");			
			system("mkdir -p ../10_FFPE_fusion/${sampleName}_STAR_FUSION");
			system("sed -i '1d' Chimeric.out.junction");
			system("mv ${sampleName}_fusion_sort.bam.bai ${sampleName}_fusion_sort.bam Chimeric.out.junction ../10_FFPE_fusion/${sampleName}_STAR_FUSION");
			#chdir("${sampleName}_STAR_FUSION");
			$cmdStr="/opt/miniconda3/envs/RNAseq/bin/STAR-Fusion --genome_lib_dir /export/md1/databases/_RNAseq_related/GRCh38_v27_CTAT_lib_Feb092018/ctat_genome_lib_build_dir/ -J ${sampleName}_STAR_FUSION/Chimeric.out.junction --output_dir ../10_FFPE_fusion/${sampleName}_geneFusion_outdir > ${sampleName}_SF.log";######
			print " $cmdStr\n";
			#system("$cmdStr");
			#system("cat ../10_FFPE_fusion/${sampleName}_geneFusion_outdir/star-fusion.fusion_predictions.tsv | head -n 20 > ../10_FFPE_fusion/${sampleName}_geneFusion_outdir/star-fusion.fusion_predictions_top20.tsv");
			$cmdStr="$arribaPATH/draw_fusions.R --fusions=../10_FFPE_fusion/${sampleName}_geneFusion_outdir/star-fusion.fusion_predictions.tsv --alignments=../10_FFPE_fusion/${sampleName}_geneFusion_outdir/${sampleName}_fusion_sort.bam --output=${sampleName}.fusions.pdf --annotation=$starFusionRef/ref_annot.gtf --cytobands=$arribaPATH/database/cytobands_hg38_GRCh38_v2.0.0.tsv --proteinDomains=$arribaPATH/database/protein_domains_hg38_GRCh38_v2.0.0.gff3";
			print "$cmdStr\n";
			#system("$cmdStr");
		}
		
		#if (!-e "../6_Report/Gene_Fusion_FFPE/${sampleName}.fusion_predictions.tsv")
		#{
		#	$cmdStr="cp ../10_FFPE_fusion/${sampleName}_geneFusion_outdir/star-fusion.fusion_predictions.tsv ../6_Report/Gene_Fusion_FFPE/${sampleName}.fusion_predictions.tsv";
		#	print "$cmdStr\n";
		#	system("$cmdStr");
		#	$cmdStr="cp ../10_FFPE_fusion/${sampleName}_geneFusion_outdir/${sampleName}.fusions.pdf ../6_Report/Gene_Fusion_FFPE/";
		#	print "$cmdStr\n";
		#	system("$cmdStr");			
		#	if (!-e "../6_Report/Gene_Fusion_FFPE/${sampleName}.fusion_predictions.xlsx")
		#	{
		#		&tabtoExcel("../6_Report/Gene_Fusion_FFPE/${sampleName}.fusion_predictions.tsv","../6_Report/Gene_Fusion_FFPE/${sampleName}.fusion_predictions.xlsx");	
		#	}
		#}
		#if (!-e "../6_Report/Gene_Fusion_FFPE/${sampleName}.fusion_predictions.abridged.tsv")
		#{
		#	$cmdStr="cp ../10_FFPE_fusion/${sampleName}_geneFusion_outdir/star-fusion.fusion_predictions.abridged.tsv ../6_Report/Gene_Fusion_FFPE/${sampleName}.fusion_predictions.abridged.tsv";
		#	print "$cmdStr\n";
		#	system("$cmdStr");
		#	if (!-e "../6_Report/Gene_Fusion_FFPE/${sampleName}.fusion_predictions.abridged.xlsx")
		#	{
		#		&tabtoExcel("../6_Report/Gene_Fusion_FFPE/${sampleName}.fusion_predictions.abridged.tsv","../6_Report/Gene_Fusion_FFPE/${sampleName}.fusion_predictions.abridged.xlsx");	
		#	}
		#}
	}
}

sub circ_RNA_CIRIquant{
	my ($ciriqPATH,%sample2fastq)=@_;#,@Sample_Name
	system("mkdir -p ../11_CIRC_ciri");
	system("mkdir -p ../6_Report/ciri");
	my $cmdStr;
	my @Sample_Name= keys %sample2fastq;
	for ( my $i=0;$i<=$#Sample_Name;$i++)
	{
		my $sampleName=$Sample_Name[$i];
		$cmdStr="$ciriqPATH/CIRIquant -t 24 -1 ${sampleName}_R1.fastq.gz -2 ${sampleName}_R2.fastq.gz --config /home/shihjielai/RNAseq_PIP/GRCh38.yml -o ${sampleName}_ciri_bam -p ${sampleName} -e ${sampleName}_ciri.log --bam ${sampleName}.sorted.bam";
		print "$cmdStr\n";
		system("$cmdStr");
	}
	if( !-e "sample2GroupInfo.txt"){
		open(IN2, "<", "groupCompareInfo.txt") or die "Could not open \n\n";
		my @IN2=<IN2>;close (IN2);
		foreach my $ab(@IN2){
			chomp $ab;
			my @Temp2=split /\t/, $ab;
			$cmdStr="$ciriqPATH/CIRI_DE -n $Temp2[1]/$Temp2[1].gtf -c $Temp2[0]/$Temp2[0].gtf -o ../11_CIRC_ciri/$Temp2[0]_vs_$Temp2[1].tsv";
			print "$cmdStr\n";
			system("$cmdStr");
			&tabtoExcel("../11_CIRC_ciri/$Temp2[0]_vs_$Temp2[1].tsv","../6_Report/ciri/$Temp2[0]_vs_$Temp2[1].xlsx");	
		}
	}else{
		open(IN2, "<", "groupCompareInfo.txt") or die "Could not open \n\n";
		my @IN2=<IN2>;close (IN2);
		open(IN4, "<", "sample2GroupInfo.txt");
		my @IN4=<IN4>;close (IN4);
		my %hash_s=();
		foreach (@IN4){
			chomp $_;
			my @Temp=split /\t/, $_;
			push @{$hash_s{$Temp[1]}}, $Temp[0];
		}
		foreach my $ab(@IN2){
			my @Temp2=split /\t/, $ab;
			open ( Genelist, ">", "sample.lst");
			open ( Genelist2, ">", "sample_gene.lst");
			my $num=1;
			foreach my $ab1(@{$hash_s{$Temp2[0]}}){
				print Genelist "$ab1 ./$ab1/${ab1}.gtf T $num\n";
				print Genelist2 "$ab1 ./$ab1/gene/${ab1}_out.gtf\n";
				$num=$num+1;
			}
			$num=1;
			foreach my $ab1(@{$hash_s{$Temp2[1]}}){
				print Genelist "$ab1 ./$ab1/${ab1}.gtf C $num\n";
				print Genelist2 "$ab1 ./$ab1/gene/${ab1}_out.gtf\n";
				$num=$num+1;
			}
			close Genelist;close Genelist2;
			$cmdStr="$ciriqPATH/prep_CIRIquant -i sample.lst --lib ../11_CIRC_ciri/$Temp2[0]_$Temp2[1]_library_info.csv --circ ../11_CIRC_ciri/$Temp2[0]_$Temp2[1]_circRNA_info.csv --bsj ../11_CIRC_ciri/$Temp2[0]_$Temp2[1]_circRNA_bsj.csv --ratio ../11_CIRC_ciri/$Temp2[0]_$Temp2[1]_circRNA_ratio.csv";
			print "$cmdStr\n";
			system("$cmdStr");
			system("prepDE.py -i sample_gene.lst");	
			$cmdStr="$ciriqPATH/CIRI_DE_replicate --lib ../11_CIRC_ciri/$Temp2[0]_$Temp2[1]_library_info.csv --bsj ../11_CIRC_ciri/$Temp2[0]_$Temp2[1]_circRNA_bsj.csv --gene gene_count_matrix.csv --our ../11_CIRC_ciri/$Temp2[0]_$Temp2[1]_circRNA_de.tsv";
			print "$cmdStr\n";
			system("$cmdStr");
			&tabtoExcel("../11_CIRC_ciri/$Temp2[0]_$Temp2[1]_circRNA_de.tsv","../6_Report/ciri/$Temp2[0]_$Temp2[1]_circRNA_de.xlsx");		
		}
	}
}


sub rmats_AS{
	my ($rmats_pip,$STARPATH,$starFusionRef,%hash_s)=@_;
	my $localPATH=`pwd`;
	chdir("..");
	my $localPATH2=`pwd`;
	chdir($localPATH);
	system("mkdir -p ../12_rmats_AS");
	system("mkdir -p ../6_Report/rmats_AS");	

	if( -e "sample2GroupInfo.txt" & ! -e "b1.txt"){
		open(IN2, "<", "groupCompareInfo.txt") or die "Could not open \n\n";
		#open(IN3, "<", "sample2GroupInfo.txt") or die "Could not open \n\n";
		my @IN2=<IN2>;close (IN2);
		#my @IN3=<IN3>;close (IN3);
		#my %has_group;
		my $cmdStr;
		foreach my $ab(@IN2){
			chomp $ab;
			my @Temp2=split /\t/, $ab;
			open ( bamlist1, ">", "../12_rmats_AS/$Temp2[0].txt");
			open ( bamlist2, ">", "../12_rmats_AS/$Temp2[1].txt");
			my $ii=0;
			my $sample1;
			foreach my $sampleName (@{$hash_s{$Temp2[0]}}){
				print "AAA$sampleName\t$Temp2[0]\n";
				if( !-e "../12_rmats_AS/${sampleName}_fusion.bam" && "../7_arriba/${sampleName}.Aligned.sortedByCoord.out.bam"){
					$cmdStr="$STARPATH/STAR --runThreadN 16 --genomeDir $starFusionRef/ref_genome.fa.star.idx/ --genomeLoad NoSharedMemory --readFilesIn ${sampleName}_R1.fastq.gz ${sampleName}_R2.fastq.gz --readFilesCommand zcat --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 --outFilterMultimapNmax 100 --outFilterMismatchNmax 33 --seedSearchStartLmax 12 --alignSJoverhangMin 15 --outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0.3 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 1 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --outReadsUnmapped None --chimOutJunctionFormat 1 --alignSJDBoverhangMin 10 --alignMatesGapMax 200000 --alignIntronMax 200000 --outSAMstrandField intronMotif --outSAMattrRGline ID:GRPundef --chimMultimapScoreRange 10 --chimNonchimScoreDropMin 10 --chimMultimapNmax 50 --outFilterType BySJout --alignSplicedMateMapLminOverLmate 0.5 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --twopassMode Basic  > ../12_rmats_AS/${sampleName}_fusion.bam";
					print "$cmdStr\n";
					#system("$cmdStr");
					$cmdStr="samtools sort -\@ 24 ../12_rmats_AS/${sampleName}_fusion.bam -o ../12_rmats_AS/${sampleName}_sort_fusion.bam.bai";
					print "$cmdStr\n";
					#system("$cmdStr");
					$cmdStr="samtools index ../12_rmats_AS/${sampleName}_sort_fusion.bam ../12_rmats_AS/${sampleName}_sort_fusion.bam.bai";
					print "$cmdStr\n";
					#system("$cmdStr");
				}

				if(!-e "../7_arriba/${sampleName}.Aligned.sortedByCoord.out.bam"){
					if($ii==0){
						$sample1="$localPATH2/12_rmats_AS/${sampleName}_sort_fusion.bam";
						print bamlist1 "$localPATH2/12_rmats_AS/${sampleName}_sort_fusion.bam";
						$ii=$ii+1;
					}else{
						$sample1=$sample1.",$localPATH2/12_rmats_AS/${sampleName}_sort_fusion.bam";
						print bamlist1 ",$localPATH2/12_rmats_AS/${sampleName}_sort_fusion.bam";	
					}		
				}else{
					if($ii==0){
						$sample1="$localPATH2/7_arriba/${sampleName}.Aligned.sortedByCoord.out.bam";
						print bamlist1 "$localPATH2/7_arriba/${sampleName}.Aligned.sortedByCoord.out.bam";
						$ii=$ii+1;
					}else{
						$sample1=$sample1.",$localPATH2/7_arriba/${sampleName}.Aligned.sortedByCoord.out.bam";
						print bamlist1 ",$localPATH2/7_arriba/${sampleName}.Aligned.sortedByCoord.out.bam";	
					}
				}
			}
			$ii=0;
			my $sample2;
			foreach my $sampleName (@{$hash_s{$Temp2[1]}}){
				print "AAA$sampleName\t$Temp2[1]\n";
				if( !-e "../12_rmats_AS/${sampleName}_fusion.bam" && "../7_arriba/${sampleName}.Aligned.sortedByCoord.out.bam"){
					$cmdStr="$STARPATH/STAR --runThreadN 16 --genomeDir $starFusionRef/ref_genome.fa.star.idx/ --genomeLoad NoSharedMemory --readFilesIn ${sampleName}_R1.fastq.gz ${sampleName}_R2.fastq.gz --readFilesCommand zcat --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 --outFilterMultimapNmax 100 --outFilterMismatchNmax 33 --seedSearchStartLmax 12 --alignSJoverhangMin 15 --outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0.3 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 1 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --outReadsUnmapped None --chimOutJunctionFormat 1 --alignSJDBoverhangMin 10 --alignMatesGapMax 200000 --alignIntronMax 200000 --outSAMstrandField intronMotif --outSAMattrRGline ID:GRPundef --chimMultimapScoreRange 10 --chimNonchimScoreDropMin 10 --chimMultimapNmax 50 --outFilterType BySJout --alignSplicedMateMapLminOverLmate 0.5 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --twopassMode Basic  > ../12_rmats_AS/${sampleName}_fusion.bam";
					print "$cmdStr\n";
					#system("$cmdStr");
					$cmdStr="samtools sort -\@ 24 ../12_rmats_AS/${sampleName}_fusion.bam -o ../12_rmats_AS/${sampleName}_sort_fusion.bam.bai";
					print "$cmdStr\n";
					#system("$cmdStr");
					$cmdStr="samtools index ../12_rmats_AS/${sampleName}_sort_fusion.bam ../12_rmats_AS/${sampleName}_sort_fusion.bam.bai";
					print "$cmdStr\n";
					#system("$cmdStr");
					
				}
				if(!-e "../7_arriba/${sampleName}.Aligned.sortedByCoord.out.bam"){
					if($ii==0){
						$sample1="$localPATH2/12_rmats_AS/${sampleName}_sort_fusion.bam";
						print bamlist1 "$localPATH2/12_rmats_AS/${sampleName}_sort_fusion.bam";
						$ii=$ii+1;
					}else{
						$sample1=$sample1.",$localPATH2/12_rmats_AS/${sampleName}_sort_fusion.bam";
						print bamlist1 ",$localPATH2/12_rmats_AS/${sampleName}_sort_fusion.bam";	
					}		
				}else{
					if($ii==0){
						$sample1="$localPATH2/7_arriba/${sampleName}.Aligned.sortedByCoord.out.bam";
						print bamlist1 "$localPATH2/7_arriba/${sampleName}.Aligned.sortedByCoord.out.bam";
						$ii=$ii+1;
					}else{
						$sample1=$sample1.",$localPATH2/7_arriba/${sampleName}.Aligned.sortedByCoord.out.bam";
						print bamlist1 ",$localPATH2/7_arriba/${sampleName}.Aligned.sortedByCoord.out.bam";	
					}
				}
			}
			close bamlist1;
			close bamlist2;
			my $cmsStr="python $rmats_pip/rmats.py --b1 ../12_rmats_AS/$Temp2[0].txt --b2 ../12_rmats_AS/$Temp2[1].txt --gtf $starFusionRef/ref_annot.gtf -t paired --nthread 14 --od ../12_rmats_AS/output_$Temp2[0]_$Temp2[1] --tmp ../12_rmats_AS/tmp_output_$Temp2[0]_$Temp2[1] --libType fr-firststrand --novelSS --allow-clipping --paired-stats --readLength 100 --variable-read-length";
			print "$cmsStr\n";
			#system("$cmdStr");
			if (-e "../12_rmats_AS/output_$Temp2[0]_$Temp2[1]/SE.MATS.JCEC.txt"){
				open(IN2, "<", "../12_rmats_AS/output_$Temp2[0]_$Temp2[1]/SE.MATS.JCEC.txt");
				open ( bamlist3, ">", "../12_rmats_AS/output_$Temp2[0]_$Temp2[1]/SE.MATS.JCEC.FDR0.05.txt");
				my @IN3=<IN2>;close (IN2);
				foreach my $ab(@IN3){
					chomp $ab;
					my @Temp2=split /\t/, $ab;
					if( $Temp2[0]=~/ID/ ){
						print bamlist3 "$ab\n";
					}
					if( ($Temp2[19] <0.05 && $Temp2[19]!=/NA/) &&($Temp2[22] > 0.3 || $Temp2[22] < -0.3)){
						print bamlist3 "$ab\n";
					}
				}
				close (bamlist3);
				my $cmsStr="rmats2sashimiplot --b1 $sample1 --b2 $sample2 -t SE -e ../12_rmats_AS/output_$Temp2[0]_$Temp2[1]/SE.MATS.JCEC.FDR0.05.txt --l1 $Temp2[0] --l2 $Temp2[1] --exon_s 1 --intron_s 5 -o ../12_rmats_AS/SE_$Temp2[0]_$Temp2[1]";
				print "$cmsStr";
				#system("$cmdStr");
			}
			

		}
	}
}

sub DEXseq_exon{
	my ($gtf2gff,$exoncount,$gtf,%hash_s)=@_;
	system("mkdir -p 13_DEXseq_exon");
	my @bmafile=<*.sorted.bam>;
	my $cmsStr="python $gtf2gff/dexseq_prepare_annotation.py $gtf2gff dexseq.gff";
	system("$cmdStr");
	foreach (@bmafile){
		my @Temp2=split /\t/, $_;
		my $cmsStr="python $gtf2gff/dexseq_count.py -s reverse -p yes -a 5 -f bam -r pos human_Star.gff ${_} ${Temp2[0]}_DEXseq.txt";
		system("$cmdStr");
	}
	open(IN2, "<", "groupCompareInfo.txt") or die "Could not open \n\n";
	my @IN2=<IN2>;close (IN2);
	my %has_group;
	my $localPATH=`pwd`;
	chdir("..");
	my $localPATH2=`pwd`;
	chdir($localPATH);
	foreach (@IN2){
		my @Temp2=split /\t/, $_;
		system("mkdir -p ../13_DEXseq_exon/$Temp2[0]_$Temp2[1]");
		open ( bamlist1, ">", "$Temp2[0]_$Temp2[1]_bamlist.txt");
		foreach my $nam (@{$hash_s{$Temp2[0]}}){
			system("mv ${nam}_DEXseq.txt ../13_DEXseq_exon/$Temp2[0]_$Temp2[1]");
			print bamlist1 "$localPATH2/13_DEXseq_exon/$Temp2[0]_$Temp2[1]/${nam}_DEXseq.txt\n";
		}

		foreach my $nam (@{$hash_s{$Temp2[1]}}){
			system("mv ${nam}_DEXseq.txt ../13_DEXseq_exon/$Temp2[0]_$Temp2[1]");
			print bamlist1 "$localPATH2/13_DEXseq_exon/$Temp2[0]_$Temp2[1]/${nam}_DEXseq.txt\n";
		}
		close bamlist1;
		if(!$hash_s{$Temp2[0]}){push @{$hash_s{$Temp2[0]}}, $Temp2[0];};
		if(!$hash_s{$Temp2[1]}){push @{$hash_s{$Temp2[1]}}, $Temp2[1];};
		my $treat = join ',' , @{$hash_s{$Temp2[0]}};$treat=~s/,/\",\"/g;$treat=~s/-/\./g;
		my $cont = join ',' , @{$hash_s{$Temp2[1]}};$cont=~s/,/\",\"/g;$cont=~s/-/\./g;
		open ( bamlist1, ">", "$Temp2[0]_$Temp2[1]_bamlist.r");
		print bamlist1 "

		countFiles <- read.delim(\"$Temp2[0]_$Temp2[1]_bamlist.txt\",header=F)
		group = c( rep(\"treated\",$#{$hash_s{$Temp2[0]}}+1),rep(\"control\",$#{$hash_s{$Temp2[1]}}+1))
		
		sampleTable = data.frame(
  		row.names = c(\"$treat\",\"$cont\"),
  		condition = group)
		sampleTable

		library( \"DEXSeq\" )

		dxd = DEXSeqDataSetFromHTSeq(
  		countFiles,
  		sampleData=sampleTable,
  		design= ~ sample + exon + condition:exon,
  		flattenedfile=flattenedFile )

		colData(dxd)

		head( counts(dxd), 5 )

		split( seq_len(ncol(dxd)), colData(dxd)\$exon )
		head( featureCounts(dxd), 5 )
		head( rowRanges(dxd), 3 )
		sampleAnnotation( dxd )
		dxd = estimateSizeFactors( dxd )
		dxd = estimateDispersions( dxd )
		dxd = testForDEU( dxd )


		dxd = estimateExonFoldChanges( dxd, fitExpToVar=\"condition\",denominator=\"control\")
		dxr1 = DEXSeqResults( dxd )

		aa<-aaaM[which(aaaM\$padj<0.05),]
		aa1<-aa\$groupID[order(aa\$padj)]
		data5 = unique(aa1)
		png(filename =paste0(data5[1],\"_DEXSeq.png\"),width = 9000, height = 5000,res=500)
		plotDEXSeq( dxr1, data5[1], legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, splicing=T,norCounts=T, FDR=0.05, displayTranscripts=TRUE,names=T)
		dev.off()
		png(filename =paste0(data5[2],\"_DEXSeq.png\"),width = 9000, height = 5000,res=500)
		plotDEXSeq( dxr1, data5[2], legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, splicing=T,norCounts=T, FDR=0.05, displayTranscripts=TRUE,names=T)
		dev.off()
		png(filename =paste0(data5[3],\"_DEXSeq.png\"),width = 9000, height = 5000,res=500)
		plotDEXSeq( dxr1, data5[3], legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, splicing=T,norCounts=T, FDR=0.05, displayTranscripts=TRUE,names=T)
		dev.off()
		png(filename =paste0(data5[4],\"_DEXSeq.png\"),width = 9000, height = 5000,res=500)
		plotDEXSeq( dxr1, data5[4], legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, splicing=T,norCounts=T, FDR=0.05, displayTranscripts=TRUE,names=T)
		dev.off()
		png(filename =paste0(data5[5],\"_DEXSeq.png\"),width = 9000, height = 5000,res=500)
		plotDEXSeq( dxr1, data5[5], legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, splicing=T,norCounts=T, FDR=0.05, displayTranscripts=TRUE,names=T)
		dev.off()";
	close bamlist1;
	#system("Rscript $Temp2[0]_$Temp2[1]_bamlist.r");
	}
}





return 1;
