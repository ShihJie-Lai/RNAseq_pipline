
sub outexcel {
	use strict;
	use Spreadsheet::ParseExcel;
	use Term::ANSIColor;
	use Excel::Writer::XLSX;
	my ($ensembl_file,$geneNFfile,$geneNFfile1,$XLSX2,$gene,$expr,$FC,%hash_a)= @_; 
	#ensembl_file: annotation data；geneNFfile:input all data；geneNFfile1: input FC"$FC"X data；XLSX2:EXCEL name；expr: FPKM or TPM；$FC:FC；hash_a: stringtie output locus；
	
	my @annInfoContent=();
	my %hash_ensembl;
	my %hash_enrich;
	my %hash_enrich2;
	open(ensembl_file,"<","$ensembl_file") or die "Could not open $ensembl_file!\n";
	@annInfoContent=<ensembl_file>; close ensembl_file;
	
	my $exp_n;
	if($expr==1){
		$exp_n="FPKM";
	}else{
		$exp_n="TPM";
	}
	my $loc;
	my $loc_R="R_dataFiles";
	if(! -d "R_dataFiles"){
		mkdir ("R_dataFiles");
	}
	
	###### annotation #########
	foreach (@annInfoContent)
	{
		$_=~s/\r\n?/\n/g;chomp $_;
		my @annItems=split /\t/, $_;
		if ($gene==0){
			$hash_ensembl{$annItems[0]}=join("\t", @annItems[2..6]);
			$hash_enrich{$annItems[0]}=$annItems[5];
			$hash_enrich2{$annItems[0]}=$annItems[5];
			$loc="DE_Gene";
			if(! -d "DE_Gene"){
				mkdir ("DE_Gene");
			}
		}else{
			$hash_ensembl{$annItems[1]}=join("\t","$annItems[0]",@annItems[2..6]);
			$hash_enrich{$annItems[1]}=$annItems[5];
			$hash_enrich2{$annItems[1]}=$annItems[5];
			$loc="DE_TS";
			if(! -d "DE_TS"){
				mkdir ("DE_TS");	
			}
		}	
	}
	
	########## combine annotation & value ###############
	my %hash_data=();
	my @solar=(); 
	my $sample_n=();
	my $sample_n1=();
	if (!-e $geneNFfile){ print "$geneNFfile not found!! Seems no gene passing the filters for this set.\n"; print "$geneNFfile not found!! Seems no gene passing the filters for this set.\n"; next;}
	open(WriteData,"<",$geneNFfile) or die "!!Cannot open file: $geneNFfile!!\n";
	my @Writedata=<WriteData>;
	close WriteData;
	for (my $i=0;$i<=$#Writedata; $i++){#$#Writedata
		chomp $Writedata[$i];
		my @row_data=split /\t/, $Writedata[$i];
		if($row_data[0]=~/^gene_id/){
			if ($gene==0){
				$hash_data{$row_data[0]}=join("\t","Protein stable ID","Gene description","Gene name","NCBI gene ID","Gene type","locus",@row_data[1..($#row_data)]);
			}else{
				$row_data[0]="Transcript_id";
				$hash_data{$row_data[0]}=join("\t","Gene stable ID","Protein stable ID","Gene description","Gene name","NCBI gene ID","Gene type","locus",@row_data[1..($#row_data)]);
			}
			$sample_n=$#row_data;
			my @row_n=split /\t/, $hash_data{$row_data[0]};
			$sample_n1=$#row_n-$sample_n+2;
		}else{
			my $row_data2;
			
			if ($gene==0){
				if (!$hash_ensembl{$row_data[0]}){$hash_ensembl{$row_data[0]}="-\t-\t-\t-\t-";$hash_enrich2{$row_data[0]}="-";}
				$row_data2=$hash_a{$row_data[0]}[1];
			}else{
				if (!$hash_ensembl{$row_data[0]}){$hash_ensembl{$row_data[0]}="-\t-\t-\t-\t-\t-";$hash_enrich2{$row_data[0]}="-";}			
				$row_data2=$hash_a{$row_data[0]}[1];
			}
			$hash_data{$row_data[0]}=join("\t",$hash_ensembl{$row_data[0]},$row_data2,@row_data[1..($#row_data)]);
			$hash_enrich2{$row_data[0]}="$hash_enrich2{$row_data[0]}\t$row_data[($#row_data-3)]";
		}
		push @solar, $row_data[0];
	}
	my %hash_data_FC2X=();

	my @solar2;
	if (!-e $geneNFfile1){ print "$geneNFfile1 not found!! Seems no gene passing the filters for this set.\n"; print "$geneNFfile1 not found!! Seems no gene passing the filters for this set.\n"; next;}
	open(WriteData1,"<",$geneNFfile1) or die "!!Cannot open file: $geneNFfile1!!\n";
	my @Writedata1=<WriteData1>;
	close WriteData1;
	for (my $i=0;$i<=$#Writedata1; $i++){##$#Writedata
		chomp $Writedata1[$i];
		my @row_data=split /\t/, $Writedata1[$i];
		if($row_data[0]=~/^gene_id/){
			if ($gene==0){
				$hash_data_FC2X{$row_data[0]}=join("\t","Protein stable ID","Gene description","Gene name","NCBI gene ID","Gene type","locus",@row_data[1..($#row_data)]);
			}else{
				$row_data[0]="Transcript_id";
				$hash_data_FC2X{$row_data[0]}=join("\t","Gene stable ID","Protein stable ID","Gene description","Gene name","NCBI gene ID","Gene type","locus",@row_data[1..($#row_data)]);
			}
		}else{
			my $row_data2;
			if ($gene==0){
				if (!$hash_enrich{$row_data[0]}){$hash_enrich{$row_data[0]}="-";}
				$row_data2=$hash_a{$row_data[0]}[1];
			}else{
				if (!$hash_enrich{$row_data[0]}){$hash_enrich{$row_data[0]}="-";}
				$row_data2=$hash_a{$row_data[0]}[1];
			}
			$hash_data_FC2X{$row_data[0]}=join("\t",$hash_ensembl{$row_data[0]},$row_data2,@row_data[1..($#row_data)]);
			$hash_enrich{$row_data[0]}="$hash_enrich{$row_data[0]}\t$row_data[($#row_data-3)]";
		}
		push @solar2, $row_data[0];
	}
	
	################# to excel - set  ################################
	my $XLSX1="${XLSX2}_FC${FC}X.xlsx";####
	my $workbook=Excel::Writer::XLSX->new("$XLSX1");
	my $header_format_1_merge = $workbook->add_format(align   => 'center',font  => 'Arial');
	$header_format_1_merge->set_right();
	$header_format_1_merge->set_bottom();
	$header_format_1_merge->set_bold();
	$header_format_1_merge->set_color( 'red' );
	$header_format_1_merge->set_size(10);
 
	my $header_format_1 = $workbook->add_format(	
		font  => 'Arial'
	);
	$header_format_1->set_bold();
	$header_format_1->set_right();
	$header_format_1->set_bottom();
	$header_format_1->set_top();
	$header_format_1->set_color( 'red' );
	$header_format_1->set_size(10);
    
	my $header_format_2 = $workbook->add_format(	
		font  => 'Arial'
	);
	$header_format_2->set_border();
	$header_format_2->set_size(10);
	
	################### input data ###############################
	my $worksheet = $workbook->add_worksheet("Gene_Expression");
	$worksheet->set_column(0, 0, 12); # Column  A
	$worksheet->set_column(1, 1, 12); # Column  B
	$worksheet->set_column(2, 3, 15); # Column  C
	$worksheet->set_column(4, 9, 10); # Column  D-J
	$worksheet->set_column(10, 10, 12); # Column  D-J
	
	$worksheet->merge_range(0,0,0,($sample_n1-1), '',$header_format_1_merge);
	$worksheet->merge_range(0,$sample_n1,0,$sample_n1+($sample_n-6), $exp_n,$header_format_1_merge);
	$worksheet->merge_range(0,$sample_n1+($sample_n-6)+1,0,$sample_n1+($sample_n-4), 'Ratio & Log2 Ratio',$header_format_1_merge);
	$worksheet->merge_range(0,$sample_n1+($sample_n-4)+1,0,$sample_n1+($sample_n-1), 'Statistics',$header_format_1_merge);
	
	for (my $i=0;$i<=$#solar; $i++){
		my @row_data=split /\t/, $hash_data{$solar[$i]};
		if($i==0){
			$worksheet->write($i+1, 0, $solar[$i],$header_format_1);
			for (my $j=0;$j<=$#row_data;$j++){
				$worksheet->write($i+1, $j+1, $row_data[$j],$header_format_1);
			}
		}else{
			$worksheet->write($i+1, 0, $solar[$i],$header_format_2);
			for (my $j=0;$j<=$#row_data;$j++){
				$worksheet->write($i+1, $j+1, $row_data[$j],$header_format_2);
			}
		}
	}
	
	my $worksheetDE = $workbook->add_worksheet("DE_Genes_FC${FC}X");###${FC_Cut}
	
	$worksheetDE->set_column(0, 0, 15); # Column  A
	$worksheetDE->set_column(1, 1, 15); # Column  B
	$worksheetDE->set_column(2, 3, 15); # Column  C
	$worksheetDE->set_column(4, 9, 10); # Column  D-J
	$worksheetDE->set_column(10, 10, 12); # Column  D-J

	$worksheetDE->merge_range(0,0,0,($sample_n1-1), '',$header_format_1_merge);
	$worksheetDE->merge_range(0,$sample_n1,0,$sample_n1+($sample_n-6), $exp_n,$header_format_1_merge);#6=statistic(3) + ratio(2) +1
	$worksheetDE->merge_range(0,$sample_n1+($sample_n-6)+1,0,$sample_n1+($sample_n-4), 'Ratio & Log2 Ratio',$header_format_1_merge); #4=statistic+1
	$worksheetDE->merge_range(0,$sample_n1+($sample_n-4)+1,0,$sample_n1+($sample_n-1), 'Statistics',$header_format_1_merge);
	
	for (my $i=0;$i<=$#solar2; $i++){
		my @row_data=split /\t/, $hash_data_FC2X{$solar2[$i]};
		if($i==0){
			$worksheetDE->write($i+1, 0, $solar2[$i],$header_format_1);
			for (my $j=0;$j<=$#row_data;$j++){
				$worksheetDE->write($i+1, $j+1, $row_data[$j],$header_format_1);
			}
		}else{
			$worksheetDE->write($i+1, 0, $solar2[$i],$header_format_2);
			for (my $j=0;$j<=$#row_data;$j++){
				$worksheetDE->write($i+1, $j+1, $row_data[$j],$header_format_2);
			}
		}
	}
	
	###### input plot #######
	my $worksheetplot = $workbook->add_worksheet("plot");
	my @png=<${XLSX2}*.png>;
	my $l=0;
	foreach my $plot(@png){
		$worksheetplot->insert_image(0,$l,$plot,10,10);
		$l=$l+10;
    }
	$workbook->close();
	############################
	system ("mv $XLSX1 $loc/$XLSX1");
	
	open(Temp1, ">$loc_R/KK_${XLSX2}_FC${FC}X.pathview") or die;
	my $input_ID=",\"-\"";
	for (my $i=1;$i<=$#solar2; $i++){
		my @Temp=split /\t/,$hash_enrich{$solar2[$i]};
		if ($input_ID!~/$Temp[0]/){
			print Temp1 $hash_enrich{$solar2[$i]}."\n";
			$input_ID=$input_ID.",\"".$Temp[0]."\"";
		}
	}
	close Temp1;

	open(Temp1, ">$loc_R/All_${XLSX2}_FC${FC}X.pathview") or die;
	my $input_ID1=",\"-\"";
	for (my $i=1;$i<=$#solar; $i++){
		my @Temp=split /\t/,$hash_enrich2{$solar[$i]};
		if ($input_ID1!~/$Temp[0]/){
			print Temp1 $hash_enrich2{$solar[$i]}."\n";
			$input_ID1=$input_ID1.",\"".$Temp[0]."\"";
		}
	}
	close Temp1;

	$input_ID=~s/,//;
	return $input_ID;
};

sub clusterProfiler_GO_KEGG_ENSG_ID{ # csv & txt to excel
	my ($input_ID, $OrgDb,$EGO,$Pathway_SP,$FC) = @_; 
	#input_ID: gene list ； OrgDb:R genome database；EGO: R name； Pathway_SP:KEGG geneome id；FC:FC
	my $CP= $EGO."_FC${FC}X";
	my $CP1=$OrgDb;
	my $loc_R="R_dataFiles";
	my @Temp=split /\./,$CP1;
	$CP1="$Temp[0].$Temp[1].$Temp[2]";
	open(R_enrich, ">${EGO}_FC${FC}X.r") or die;
        print R_enrich "
		library(clusterProfiler)
		library(enrichplot)
		library(DOSE)
		library($OrgDb)
		library(ggplot2)
		gene<-c($input_ID)\n
		
		egoMF <- enrichGO(gene= gene,OrgDb=$OrgDb,keyType=\"ENTREZID\",ont= \"MF\",pAdjustMethod = \"BH\",pvalueCutoff  = 1,qvalueCutoff  = 1,readable=TRUE)
		egoBP <- enrichGO(gene= gene,OrgDb=$OrgDb,keyType=\"ENTREZID\",ont= \"BP\",pAdjustMethod = \"BH\",pvalueCutoff  = 1,qvalueCutoff  = 1,readable=TRUE)
		egoCC <- enrichGO(gene= gene,OrgDb=$OrgDb,keyType=\"ENTREZID\",ont= \"CC\",pAdjustMethod = \"BH\",pvalueCutoff  = 1,qvalueCutoff  = 1,readable=TRUE)
		kk <- enrichKEGG(gene= gene,organism= \"$Pathway_SP\",pvalueCutoff = 1,qvalueCutoff=1,keyType = \"ncbi-geneid\")
		
		FC<-read.delim(\"KK_$CP.pathview\",header=F,row.names=1)
		genelist1<-FC\$V2
		names(genelist1)<-row.names(FC)
		
		################################################
		b<-data.frame()
		for (i in c(1:10)){
		  if(length(egoMF\@result\$ID)<i){break;}
		  a<-data.frame(egoMF\@geneSets[egoMF\@result\$ID[i]])
	  	row.names(a)<-a[,1]
		  aa<-merge(a,FC,by=\"row.names\")
		  b<-rbind(b,matrix(c(egoMF\@result\$Description[i],nrow(aa[aa\$V2>0,]),\"UP\"),1,3))
		  b<-rbind(b,matrix(c(egoMF\@result\$Description[i],nrow(aa[aa\$V2<0,]),\"DOWN\"),1,3))
		}
		b\$V2<-as.numeric(as.character(b\$V2))
		cc<-max(b\$V2)
		if(cc<5){
			cc2<-seq(-cc,cc,1)
		}else{
			cc1<-(cc*2)/10
			cc2<-round(seq(-cc,cc,cc1),0)
		}
		g<-ggplot(data = b, aes(x = V1, fill = V3)) +
		  geom_bar(stat = \"identity\", data = subset(b, V3 == \"UP\"), aes(y = V2)) +
		  geom_bar(stat = \"identity\", data = subset(b, V3 == \"DOWN\"), aes(y=V2 * (-1)) ) +
		  scale_y_continuous(limits = c(-cc, cc), breaks = cc2, labels = abs(cc2)) + 
		  theme(axis.text = element_text(colour = \"black\")) + 
		  coord_flip()+ylab(\"\") + xlab(\"\")+ guides(fill = guide_legend(title = \"\"))
		png(filename = \"MF_".$CP."_py.png\",width = 2000, height = 1600,res=100)
		print(g)
		dev.off()
		#ggsave(file=\"MF_".$CP."_py.svg\",plot=g,width=10,height=8)
		##################################################
		b<-data.frame()
		for (i in c(1:10)){
		  if(length(egoBP\@result\$ID)<i){break;}
		  a<-data.frame(egoBP\@geneSets[egoBP\@result\$ID[i]])
	  	row.names(a)<-a[,1]
		  aa<-merge(a,FC,by=\"row.names\")
		  b<-rbind(b,matrix(c(egoBP\@result\$Description[i],nrow(aa[aa\$V2>0,]),\"UP\"),1,3))
		  b<-rbind(b,matrix(c(egoBP\@result\$Description[i],nrow(aa[aa\$V2<0,]),\"DOWN\"),1,3))
		}
		b\$V2<-as.numeric(as.character(b\$V2))
		cc<-max(b\$V2)
		if(cc<5){
			cc2<-seq(-cc,cc,1)
		}else{
			cc1<-(cc*2)/10
			cc2<-round(seq(-cc,cc,cc1),0)
		}
		g<-ggplot(data = b, aes(x = V1, fill = V3)) +
		  geom_bar(stat = \"identity\", data = subset(b, V3 == \"UP\"), aes(y = V2)) +
		  geom_bar(stat = \"identity\", data = subset(b, V3 == \"DOWN\"), aes(y=V2 * (-1)) ) +
		  scale_y_continuous(limits = c(-cc, cc), breaks =cc2, labels = abs(cc2)) + 
		  theme(axis.text = element_text(colour = \"black\")) + 
		  coord_flip()+ylab(\"\") + xlab(\"\")+ guides(fill = guide_legend(title = \"\"))
		png(filename = \"BP_".$CP."_py.png\",width = 2000, height = 1600,res=100)
		print(g)
		dev.off()
		#ggsave(file=\"BP_".$CP."_py.svg\",plot=g,width=10,height=8)
		##################################################
		b<-data.frame()
		for (i in c(1:10)){
		  if(length(egoCC\@result\$ID)<i){break;}
		  a<-data.frame(egoCC\@geneSets[egoCC\@result\$ID[i]])
	  	row.names(a)<-a[,1]
		  aa<-merge(a,FC,by=\"row.names\")
		  b<-rbind(b,matrix(c(egoCC\@result\$Description[i],nrow(aa[aa\$V2>0,]),\"UP\"),1,3))
		  b<-rbind(b,matrix(c(egoCC\@result\$Description[i],nrow(aa[aa\$V2<0,]),\"DOWN\"),1,3))
		}
		b\$V2<-as.numeric(as.character(b\$V2))
		cc<-max(b\$V2)
		if(cc<5){
			cc2<-seq(-cc,cc,1)
		}else{
			cc1<-(cc*2)/10
			cc2<-round(seq(-cc,cc,cc1),0)
		}
		g<-ggplot(data = b, aes(x = V1, fill = V3)) +
		  geom_bar(stat = \"identity\", data = subset(b, V3 == \"UP\"), aes(y = V2)) +
		  geom_bar(stat = \"identity\", data = subset(b, V3 == \"DOWN\"), aes(y=V2 * (-1)) ) +
		  scale_y_continuous(limits = c(-cc, cc), breaks = cc2, labels = abs(cc2)) + 
		  theme(axis.text = element_text(colour = \"black\")) + 
		  coord_flip()+ylab(\"\") + xlab(\"\")+ guides(fill = guide_legend(title = \"\"))
		png(filename = \"CC_".$CP."_py.png\",width = 2000, height = 1600,res=100)
		print(g)
		dev.off()
		#ggsave(file=\"CC_".$CP."_py.svg\",plot=g,width=10,height=8)
		##################################################
		b<-data.frame()
		for (i in c(1:10)){
		  if(length(kk\@result\$ID)<i){break;}
		  a<-data.frame(kk\@geneSets[kk\@result\$ID[i]])
		  row.names(a)<-a[,1]
		  aa<-merge(a,FC,by=\"row.names\")
		  b<-rbind(b,matrix(c(kk\@result\$Description[i],nrow(aa[aa\$V2>0,]),\"UP\"),1,3))
		  b<-rbind(b,matrix(c(kk\@result\$Description[i],nrow(aa[aa\$V2<0,]),\"DOWN\"),1,3))
		}
		b\$V2<-as.numeric(as.character(b\$V2))
		cc<-max(b\$V2)
		if(cc<5){
			cc2<-seq(-cc,cc,1)
		}else{
			cc1<-(cc*2)/10
			cc2<-round(seq(-cc,cc,cc1),0)
		}
		g<-ggplot(data = b, aes(x = V1, fill = V3)) +
		  geom_bar(stat = \"identity\", data = subset(b, V3 == \"UP\"), aes(y = V2)) +
		  geom_bar(stat = \"identity\", data = subset(b, V3 == \"DOWN\"), aes(y=V2 * (-1)) ) +
		  scale_y_continuous(limits = c(-cc, cc), breaks = cc2, labels = abs(cc2)) + 
		  theme(axis.text = element_text(colour = \"black\")) + 
		  coord_flip()+ylab(\"\") + xlab(\"\")+ guides(fill = guide_legend(title = \"\"))
		png(filename = \"KK_".$CP."_py.png\",width = 2000, height = 1600,res=100)
		print(g)
		dev.off()
		#ggsave(file=\"KK_".$CP."_py.svg\",plot=g,width=10,height=8)	
		##################################################		
		cc<-\"$Pathway_SP\"		
		if(cc==\"hsa\"){
			egoDO <- enrichDO(gene=gene,ont=\"DO\",pvalueCutoff  = 1,pAdjustMethod = \"BH\",minGSSize=5,maxGSSize=500,qvalueCutoff= 1,readable=TRUE)
			write.table(data.frame(egoDO), file = \"DO_$CP.tab\",sep=\"\\t\")
			png(filename = \"DO_".$CP."_bar.png\",width = 2000, height = 1600,res=100)
			g<-barplot(egoDO, showCategory=20)
			print(g)
			dev.off()
			#ggsave(file=\"DO_".$CP."_bar.svg\",plot=g,width=10,heigth=8)
			png(filename = \"DO_".$CP."_dot.png\",width = 2000, height = 1600,res=100)
			g<-dotplot(egoDO, showCategory=20)
			print(g)
			dev.off()
			#ggsave(file=\"DO_".$CP."_dot.svg\",plot=g,width=10,heigth=8)
			png(filename = \"DO_".$CP."_emapplot.png\",width = 3200, height = 3200,res=380)
			g<-emapplot(egoDO)
			print(g)
			dev.off()
			#ggsave(file=\"DO_".$CP."_emapplot.svg\",plot=g,width=10,heigth=8)
			png(filename = \"DO_".$CP."_cent.png\",width = 1600, height = 1600,res=100)
			g<-cnetplot(egoDO,foldChange=genelist1)
			print(g)
			dev.off()
			#ggsave(file=\"DO_".$CP."_cent.svg\",plot=g,width=10,heigth=8)
			b<-data.frame()
			for (i in c(1:10)){
			if(length(egoDO\@result\$ID)<i){break;}
			a<-data.frame(egoDO\@geneSets[egoDO\@result\$ID[i]])
			row.names(a)<-a[,1]
			aa<-merge(a,FC,by=\"row.names\")
			b<-rbind(b,matrix(c(egoDO\@result\$Description[i],nrow(aa[aa\$V2>0,]),\"UP\"),1,3))
			b<-rbind(b,matrix(c(egoDO\@result\$Description[i],nrow(aa[aa\$V2<0,]),\"DOWN\"),1,3))
			}
			b\$V2<-as.numeric(as.character(b\$V2))
			cc<-max(b\$V2)
			if(cc<5){
				cc2<-seq(-cc,cc,1)
			}else{
				cc1<-(cc*2)/10
				cc2<-round(seq(-cc,cc,cc1),0)
			}
			g<-ggplot(data = b, aes(x = V1, fill = V3)) +
			geom_bar(stat = \"identity\", data = subset(b, V3 == \"UP\"), aes(y = V2)) +
			geom_bar(stat = \"identity\", data = subset(b, V3 == \"DOWN\"), aes(y=V2 * (-1)) ) +
			scale_y_continuous(limits = c(-cc, cc), breaks = cc2, labels = abs(cc2)) + 
			theme(axis.text = element_text(colour = \"black\")) + 
			coord_flip()+ylab(\"\") + xlab(\"\")+ guides(fill = guide_legend(title = \"\"))
			png(filename = \"DO_".$CP."_py.png\",width = 2000, height = 1600,res=100)
			print(g)
			dev.off()
			#ggsave(file=\"DO_".$CP."_py.svg\",plot=g,width=10,heigth=8)
		}
		
		kk<- setReadable(kk, '$OrgDb' ,'ENTREZID' )		

		entrezIDs <- mget(gene, ".$CP1."SYMBOL , ifnotfound=NA)       #change database
		aa<-t(data.frame(entrezIDs))
		bb<-data.frame(aa,gene)
		write.table(bb,\"b.tab\",sep=\"\\t\", row.names = F,col.names = F,na = \"-\", quote = F)
		
		write.table(data.frame(egoMF), file = \"MF_$CP.tab\",sep=\"\\t\")
		write.table(data.frame(egoBP), file = \"BP_$CP.tab\",sep=\"\\t\")
		write.table(data.frame(egoCC), file = \"CC_$CP.tab\",sep=\"\\t\")
		write.table(data.frame(kk), file = \"KK_$CP.tab\",sep=\"\\t\")


		png(filename = \"MF_".$CP."_bar.png\",width = 2000, height = 1600,res=100)
		g<-barplot(egoMF, showCategory=20)
		print(g)
		dev.off()
		#ggsave(file=\"MF_".$CP."_bar.svg\",plot=g,width=10,height=8)			
		png(filename = \"BP_".$CP."_bar.png\",width = 2000, height = 1600,res=100)
		g<-barplot(egoBP, showCategory=20)
		print(g)
		dev.off()
		#ggsave(file=\"BP_".$CP."_bar.svg\",plot=g,width=10,height=8)			
		png(filename = \"CC_".$CP."_bar.png\",width = 2000, height = 1600,res=100)
		g<-barplot(egoCC, showCategory=20)
		print(g)
		dev.off()
		#ggsave(file=\"CC_".$CP."_bar.svg\",plot=g,width=10,height=8)
		png(filename = \"KK_".$CP."_bar.png\",width = 2000, height = 1600,res=100)
		g<-barplot(kk, showCategory=20)
		print(g)
		dev.off()
		#ggsave(file=\"KK_".$CP."_bar.svg\",plot=g,width=10,heigth=8)

		png(filename = \"MF_".$CP."_dot.png\",width = 2000, height = 1600,res=100)
		g<-dotplot(egoMF, showCategory=20)
		print(g)
		dev.off()
		#ggsave(file=\"MF_".$CP."_dot.svg\",plot=g,width=10,height=8)				
		png(filename = \"BP_".$CP."_dot.png\",width = 2000, height = 1600,res=100)
		g<-dotplot(egoBP, showCategory=20)
		print(g)
		dev.off()
		#ggsave(file=\"BP_".$CP."_dot.svg\",plot=g,width=10,height=8)				
		png(filename = \"CC_".$CP."_dot.png\",width = 2000, height = 1600,res=100)
		g<-dotplot(egoCC, showCategory=20)
		print(g)
		dev.off()
		#ggsave(file=\"CC_".$CP."_dot.svg\",plot=g,width=10,height=8)
		png(filename = \"KK_".$CP."_dot.png\",width = 2000, height = 1600,res=100)
		g<-dotplot(kk, showCategory=20)
		print(g)
		dev.off()
		#ggsave(file=\"KK_".$CP."_dot.svg\",plot=g,width=10,height=8)
		

		png(filename = \"MF_".$CP."_emapplot.png\",width = 3200, height = 3200,res=380)
		g<-emapplot(egoMF)
		print(g)
		dev.off()
		#ggsave(file=\"MF_".$CP."_emapplot.svg\",plot=g,width=10,height=8)			
		png(filename = \"BP_".$CP."_emapplot.png\",width = 3200, height = 3200,res=380)
		g<-emapplot(egoBP)
		print(g)
		dev.off()
		#ggsave(file=\"BP_".$CP."_emapplot.svg\",plot=g,width=10,height=8)				
		png(filename = \"CC_".$CP."_emapplot.png\",width = 3200, height = 3200,res=380)
		g<-emapplot(egoCC)
		print(g)
		dev.off()
		#ggsave(file=\"CC_".$CP."_emapplot.svg\",plot=g,width=10,height=8)
		png(filename = \"KK_".$CP."_emapplot.png\",width = 3200, height = 3200,res=380)
		g<-emapplot(kk)
		print(g)
		dev.off()
		#ggsave(file=\"KK_".$CP."_emapplot.svg\",plot=g,width=10,height=8)

		png(filename = \"MF_".$CP."_cent.png\",width = 3000, height = 3000,res=300)
		g<-cnetplot(egoMF,foldChange=genelist1)
		print(g)
		dev.off()
		#ggsave(file=\"MF_".$CP."_cent.svg\",plot=g,width=10,height=8)			
		png(filename = \"BP_".$CP."_cent.png\",width = 3000, height = 3000,res=300)
		g<-cnetplot(egoBP,foldChange=genelist1)
		print(g)
		dev.off()
		#ggsave(file=\"BP_".$CP."_cent.svg\",plot=g,width=10,height=8)			
		png(filename = \"CC_".$CP."_cent.png\",width = 3000, height = 3000,res=300)
		cnetplot(egoCC,foldChange=genelist1)
		print(g)
		dev.off()
		#ggsave(file=\"CC_".$CP."_cent.svg\",plot=g,width=10,height=8)
		png(filename = \"KK_".$CP."_cent.png\",width = 3000, height = 3000,res=300)
		g<-cnetplot(kk,foldChange=genelist1)
		print(g)
		dev.off()
		#ggsave(file=\"KK_".$CP."_cent.svg\",plot=g,width=10,height=8)

		png(filename = \"MF_".$CP."_heat.png\",width = 3000, height = 3000,res=300)
		g<-heatplot(egoMF,foldChange=genelist1)
		print(g)
		dev.off()
		#ggsave(file=\"MF_".$CP."_heat.svg\",plot=g,width=10,height=8)				
		png(filename = \"BP_".$CP."_heat.png\",width = 3000, height = 3000,res=300)
		g<-heatplot(egoBP,foldChange=genelist1)
		print(g)
		dev.off()
		#ggsave(file=\"BP_".$CP."_heat.svg\",plot=g,width=10,height=8)			
		png(filename = \"CC_".$CP."_heat.png\",width = 3000, height = 3000,res=300)
		g<-heatplot(egoCC,foldChange=genelist1)
		print(g)
		dev.off()
		#ggsave(file=\"CC_".$CP."_heat.svg\",plot=g,width=10,height=8)
		png(filename = \"KK_".$CP."_heat.png\",width = 3000, height = 3000,res=300)
		g<-heatplot(kk,foldChange=genelist1)
		print(g)
		dev.off()
		#ggsave(file=\"KK_".$CP."_heat.svg\",plot=g,width=10,height=8)

		#png(filename = \"MF_".$CP."_upset.png\",width = 3000, height = 3000,res=300)
		#upsetplot(egoMF,10)
		#dev.off()				
		#png(filename = \"BP_".$CP."_upset.png\",width = 3000, height = 3000,res=300)
		#upsetplot(egoBP,10)
		#dev.off()				
		#png(filename = \"CC_".$CP."_upset.png\",width = 3000, height = 3000,res=300)
		#upsetplot(egoCC,10)
		#dev.off()
		#png(filename = \"KK_".$CP."_upset.png\",width = 3000, height = 3000,res=300)
		#upsetplot(kk,10)
		#dev.off()						

		png(filename = \"BP_".$CP."_GOP.png\",width = 3600, height = 3600,res=400)
		plotGOgraph(egoBP)
		dev.off()			
		png(filename = \"MF_".$CP."_GOP.png\",width = 3600, height = 3600,res=400)
		plotGOgraph(egoMF)
		dev.off()				
		png(filename = \"CC_".$CP."_GOP.png\",width = 3600, height = 3600,res=400)
		plotGOgraph(egoCC)
		dev.off()

		png(filename = \"BP_".$CP."_GO.png\",width = 3600, height = 3600,res=400)
		g<-goplot(egoBP)
		print(g)
		dev.off()
		#ggsave(file=\"BP_".$CP."_GO.svg\",plot=g,width=10,height=8)				
		png(filename = \"MF_".$CP."_GO.png\",width = 3600, height = 3600,res=400)
		g<-goplot(egoMF)
		print(g)
		dev.off()
		#ggsave(file=\"MF_".$CP."_GO.svg\",plot=g,width=10,height=8)				
		png(filename = \"CC_".$CP."_GO.png\",width = 3600, height = 3600,res=400)
		g<-goplot(egoCC)
		print(g)
		dev.off()
		#ggsave(file=\"CC_".$CP."_GO.svg\",plot=g,width=10,height=8)



";
	close R_enrich;
	system ("mv ${EGO}_FC${FC}X.r $loc_R/${EGO}_FC${FC}X.r");
	
}


sub clusterProfiler_GO_KEGG_ENSG_ID2{ # csv & txt to excel
	my ($OrgDb,$EGO,$Pathway_SP,$FC) = @_; 
	#input_ID: gene list ； OrgDb:R genome database；EGO: R name； Pathway_SP:KEGG geneome id；FC:FC
	my $CP= $EGO."_FC${FC}X_GSEA";
	my $CP1=$OrgDb;
	my $loc_R="R_dataFiles";
	my @Temp=split /\./,$CP1;
	$CP1="$Temp[0].$Temp[1].$Temp[2]";
	open(R_enrich, ">${EGO}_FC${FC}X_GSEA.r") or die;
        print R_enrich "
		library(clusterProfiler)
		library(enrichplot)
		library($OrgDb)
		library(ggplot2)

		FC<-read.delim(\"All_${EGO}_FC${FC}X.pathview\",header=F,row.names=1)
		genelist1<-FC\$V2
		names(genelist1)<-row.names(FC)
		genelist1 <- sort(genelist1,decreasing=T)

		
		egoMF <- tryCatch({gseGO(geneList= genelist1,OrgDb=$OrgDb,keyType=\"ENTREZID\",ont= \"MF\",pAdjustMethod = \"BH\",pvalueCutoff  = 1,nPerm=1000,minGSSize=100,maxGSSize=500)}, error=function(e){gseGO(geneList= genelist1[which(abs(genelist1)>2)],OrgDb=$OrgDb,keyType=\"ENTREZID\",ont= \"MF\",pAdjustMethod = \"BH\",pvalueCutoff  = 1,nPerm=1000,minGSSize=100,maxGSSize=500)})
		egoBP <- tryCatch({gseGO(geneList= genelist1,OrgDb=$OrgDb,keyType=\"ENTREZID\",ont= \"BP\",pAdjustMethod = \"BH\",pvalueCutoff  = 1,nPerm=1000,minGSSize=100,maxGSSize=500)}, error=function(e){gseGO(geneList= genelist1[which(abs(genelist1)>2)],OrgDb=$OrgDb,keyType=\"ENTREZID\",ont= \"BP\",pAdjustMethod = \"BH\",pvalueCutoff  = 1,nPerm=1000,minGSSize=100,maxGSSize=500)})
		egoCC <- tryCatch({gseGO(geneList= genelist1,OrgDb=$OrgDb,keyType=\"ENTREZID\",ont= \"CC\",pAdjustMethod = \"BH\",pvalueCutoff  = 1,nPerm=1000,minGSSize=100,maxGSSize=500)}, error=function(e){gseGO(geneList= genelist1[which(abs(genelist1)>2)],OrgDb=$OrgDb,keyType=\"ENTREZID\",ont= \"CC\",pAdjustMethod = \"BH\",pvalueCutoff  = 1,nPerm=1000,minGSSize=100,maxGSSize=500)})
		kk <- tryCatch({gseKEGG(geneList= genelist1,organism= \"$Pathway_SP\",pvalueCutoff = 1,keyType = \"ncbi-geneid\",nPerm=1000,minGSSize=100,maxGSSize=500)}, error=function(e){gseKEGG(geneList= genelist1[which(abs(genelist1)>2)],organism= \"$Pathway_SP\",pvalueCutoff = 1,keyType = \"ncbi-geneid\",nPerm=1000,minGSSize=100,maxGSSize=500)})
		
		png(filename = \"MF_".$CP."_ridge.png\",width = 3000, height = 3000,res=300)
		g<-ridgeplot(egoMF,10)
		print(g)
		dev.off()
		#ggsave(file=\"MF_".$CP."_ridge.svg\",plot=g,width=10,height=8)			
		png(filename = \"BP_".$CP."_ridge.png\",width = 3000, height = 3000,res=300)
		g<-ridgeplot(egoBP,10)
		print(g)
		dev.off()				
		#ggsave(file=\"BP_".$CP."_ridge.svg\",plot=g,width=10,height=8)
		png(filename = \"CC_".$CP."_ridge.png\",width = 3000, height = 3000,res=300)
		g<-ridgeplot(egoCC,10)
		print(g)
		dev.off()
		#ggsave(file=\"CC_".$CP."_ridge.svg\",plot=g,width=10,height=8)
		png(filename = \"KK_".$CP."_ridge.png\",width = 3000, height = 3000,res=300)
		g<-ridgeplot(kk,10)
		print(g)
		dev.off()
		#ggsave(file=\"KK_".$CP."_ridge.svg\",plot=g,width=10,height=8)

		for(i in 1:5){
			png(filename =paste0(\"MF_".$CP."_gsea\",i,\".png\"),width = 3000, height = 3000,res=300)
			cc<-gseaplot2(egoMF,geneSetID=i,title=egoMF\$Description[i])
			print(cc)			
			dev.off()				
			png(filename =paste0(\"BP_".$CP."_gsea\",i,\".png\"),width = 3000, height = 3000,res=300)
			cc<-gseaplot2(egoBP,geneSetID=i,title=egoBP\$Description[i])
			print(cc)			
			dev.off()				
			png(filename =paste0(\"CC_".$CP."_gsea\",i,\".png\"),width = 3000, height = 3000,res=300)
			cc<-gseaplot2(egoCC,geneSetID=i,title=egoCC\$Description[i])
			print(cc)			
			dev.off()
			png(filename =paste0(\"KK_".$CP."_gsea\",i,\".png\"),width = 3000, height = 3000,res=300)
			cc<-gseaplot2(kk,geneSetID=i,title=kk\$Description[i])
			print(cc)			
			dev.off()
		}
		
		lapply(1:5 ,function(i) {
			anno <- egoMF[i , c(\"NES\",\"pvalue\",\"p.adjust\")]
			lab <- paste0(names(anno),\"=\",round(anno, 3), collapse=\"\n\")
			png(filename =paste0(\"MF_".$CP."_gsearank\",i,\".png\"),width = 3000, height = 3000,res=300)
			cc<-gsearank(egoMF,i,egoMF[i,2]) + annotate(\"text\",0,egoMF[i,\"enrichmentScore\"]*.9,label = lab , hjust=0,vjust=0)
			print(cc)			
			dev.off()

			anno <- egoBP[i , c(\"NES\",\"pvalue\",\"p.adjust\")]
			lab <- paste0(names(anno),\"=\",round(anno, 3), collapse=\"\n\")
			png(filename =paste0(\"BP_".$CP."_gsearank\",i,\".png\"),width = 3000, height = 3000,res=300)
			cc<-gsearank(egoBP,i,egoBP[i,2]) + annotate(\"text\",0,egoBP[i,\"enrichmentScore\"]*.9,label = lab , hjust=0,vjust=0)
			print(cc)			
			dev.off()

			anno <- egoCC[i , c(\"NES\",\"pvalue\",\"p.adjust\")]
			lab <- paste0(names(anno),\"=\",round(anno, 3), collapse=\"\n\")
			png(filename =paste0(\"CC_".$CP."_gsearank\",i,\".png\"),width = 3000, height = 3000,res=300)
			cc<-gsearank(egoCC,i,egoCC[i,2]) + annotate(\"text\",0,egoCC[i,\"enrichmentScore\"]*.9,label = lab , hjust=0,vjust=0)
			print(cc)			
			dev.off()

			anno <- kk[i , c(\"NES\",\"pvalue\",\"p.adjust\")]
			lab <- paste0(names(anno),\"=\",round(anno, 3), collapse=\"\n\")
			png(filename =paste0(\"KK_".$CP."_gsearank\",i,\".png\"),width = 3000, height = 3000,res=300)
			cc<-gsearank(kk,i,kk[i,2]) + annotate(\"text\",0,kk[i,\"enrichmentScore\"]*.9,label = lab , hjust=0,vjust=0)
			print(cc)			
			dev.off()
		})
		
		egoMF<- setReadable(egoMF, '$OrgDb' ,'ENTREZID' )
		egoBP<- setReadable(egoBP, '$OrgDb' ,'ENTREZID' )
		egoCC<- setReadable(egoCC, '$OrgDb' ,'ENTREZID' )
		kk<- setReadable(kk, '$OrgDb' ,'ENTREZID' )		

		entrezIDs <- mget(names(genelist1), ".$CP1."SYMBOL , ifnotfound=NA)       #change database
		aa<-t(data.frame(entrezIDs))
		bb<-data.frame(aa,names(genelist1))
		write.table(bb,\"b.tab\",sep=\"\\t\", row.names = F,col.names = F,na = \"-\", quote = F)
		
		write.table(data.frame(egoMF), file = \"MF_$CP.tab\",sep=\"\\t\")
		write.table(data.frame(egoBP), file = \"BP_$CP.tab\",sep=\"\\t\")
		write.table(data.frame(egoCC), file = \"CC_$CP.tab\",sep=\"\\t\")
		write.table(data.frame(kk), file = \"KK_$CP.tab\",sep=\"\\t\")

		png(filename = \"MF_".$CP."_dot.png\",width = 2000, height = 1600,res=100)
		g<-dotplot(egoMF, showCategory=20)
		print(g)
		dev.off()
		#ggsave(file=\"MF_".$CP."_dot.svg\",plot=g,width=10,height=8)			
		png(filename = \"BP_".$CP."_dot.png\",width = 2000, height = 1600,res=100)
		g<-dotplot(egoBP, showCategory=20)
		print(g)
		dev.off()
		#ggsave(file=\"BP_".$CP."_dot.svg\",plot=g,width=10,height=8)			
		png(filename = \"CC_".$CP."_dot.png\",width = 2000, height = 1600,res=100)
		g<-dotplot(egoCC, showCategory=20)
		print(g)
		dev.off()
		#ggsave(file=\"CC_".$CP."_dot.svg\",plot=g,width=10,height=8)
		png(filename = \"KK_".$CP."_dot.png\",width = 2000, height = 1600,res=100)
		g<-dotplot(kk, showCategory=20)
		print(g)
		dev.off()
		#ggsave(file=\"KK_".$CP."_dot.svg\",plot=g,width=10,height=8)		

		png(filename = \"MF_".$CP."_emapplot.png\",width = 3200, height = 3200,res=380)
		g<-emapplot(egoMF)
		print(g)
		dev.off()
		#ggsave(file=\"MF_".$CP."_emapplot.svg\",plot=g,width=10,height=8)			
		png(filename = \"BP_".$CP."_emapplot.png\",width = 3200, height = 3200,res=380)
		g<-emapplot(egoBP)
		print(g)
		dev.off()				
		#ggsave(file=\"BP_".$CP."_emapplot.svg\",plot=g,width=10,height=8)
		png(filename = \"CC_".$CP."_emapplot.png\",width = 3200, height = 3200,res=380)
		g<-emapplot(egoCC)
		print(g)
		dev.off()
		#ggsave(file=\"CC_".$CP."_emapplot.svg\",plot=g,width=10,height=8)
		png(filename = \"KK_".$CP."_emapplot.png\",width = 3200, height = 3200,res=380)
		g<-emapplot(kk)
		print(g)
		dev.off()
		#ggsave(file=\"KK_".$CP."_emapplot.svg\",plot=g,width=10,height=8)

		png(filename = \"MF_".$CP."_cent.png\",width = 3000, height = 3000,res=300)
		g<-cnetplot(egoMF,foldChange=genelist1)
		print(g)
		dev.off()
		#ggsave(file=\"MF_".$CP."_cent.svg\",plot=g,width=10,height=8)				
		png(filename = \"BP_".$CP."_cent.png\",width = 3000, height = 3000,res=300)
		g<-cnetplot(egoBP,foldChange=genelist1)
		print(g)
		dev.off()
		#ggsave(file=\"BP_".$CP."_cent.svg\",plot=g,width=10,height=8)				
		png(filename = \"CC_".$CP."_cent.png\",width = 3000, height = 3000,res=300)
		g<-cnetplot(egoCC,foldChange=genelist1)
		print(g)
		dev.off()
		#ggsave(file=\"CC_".$CP."_cent.svg\",plot=g,width=10,height=8)
		png(filename = \"KK_".$CP."_cent.png\",width = 3000, height = 3000,res=300)
		g<-cnetplot(kk,foldChange=genelist1)
		print(g)
		dev.off()
		#ggsave(file=\"KK_".$CP."_cent.svg\",plot=g,width=10,height=8)

		png(filename = \"MF_".$CP."_heat.png\",width = 3000, height = 3000,res=300)
		g<-heatplot(egoMF,foldChange=genelist1)
		print(g)
		dev.off()
		#ggsave(file=\"MF_".$CP."_heat.svg\",plot=g,width=10,height=8)			
		png(filename = \"BP_".$CP."_heat.png\",width = 3000, height = 3000,res=300)
		g<-heatplot(egoBP,foldChange=genelist1)
		print(g)
		dev.off()
		#ggsave(file=\"BP_".$CP."_heat.svg\",plot=g,width=10,height=8)				
		png(filename = \"CC_".$CP."_heat.png\",width = 3000, height = 3000,res=300)
		g<-heatplot(egoCC,foldChange=genelist1)
		print(g)
		dev.off()
		#ggsave(file=\"CC_".$CP."_heat.svg\",plot=g,width=10,height=8)	
		png(filename = \"KK_".$CP."_heat.png\",width = 3000, height = 3000,res=300)
		g<-heatplot(kk,foldChange=genelist1)
		print(g)
		dev.off()	
		#ggsave(file=\"KK_".$CP."_heat.svg\",plot=g,width=10,height=8)	

		#png(filename = \"MF_".$CP."_upset.png\",width = 3000, height = 3000,res=300)
		#upsetplot(egoMF,10)
		#dev.off()				
		#png(filename = \"BP_".$CP."_upset.png\",width = 3000, height = 3000,res=300)
		#upsetplot(egoBP,10)
		#dev.off()				
		#png(filename = \"CC_".$CP."_upset.png\",width = 3000, height = 3000,res=300)
		#upsetplot(egoCC,10)
		#dev.off()
		#png(filename = \"KK_".$CP."_upset.png\",width = 3000, height = 3000,res=300)
		#upsetplot(kk,10)
		#dev.off()						

		png(filename = \"BP_".$CP."_GOP.png\",width = 3600, height = 3600,res=400)
		plotGOgraph(egoBP)
		dev.off()				
		png(filename = \"MF_".$CP."_GOP.png\",width = 3600, height = 3600,res=400)
		plotGOgraph(egoMF)
		dev.off()				
		png(filename = \"CC_".$CP."_GOP.png\",width = 3600, height = 3600,res=400)
		plotGOgraph(egoCC)
		dev.off()
			
";
	close R_enrich;
	system ("mv ${EGO}_FC${FC}X_GSEA.r $loc_R/${EGO}_FC${FC}X_GSEA.r");
	
}

sub tabtoExcel{# csv & txt to excel
	use Spreadsheet::ParseExcel;
	use Term::ANSIColor;
	use Excel::Writer::XLSX;
	use File::Basename;

	my $loc_R="Expression_Table";
	if(! -d "Expression_Table"){
		mkdir ("Expression_Table");
	}
	my ($ensembl_file,$excelname) = @_;
	open(ensembl_file,"<",$ensembl_file) or die "Could not open $ensembl_file!\n";
	my @annInfoContent=<ensembl_file>; close ensembl_file;
	my $excelname2=basename($excelname,".xlsx");
	
	my $workbook=Excel::Writer::XLSX->new("$excelname");
	my $header_format_1 = $workbook->add_format(	
		font  => 'Arial'
	);
	$header_format_1->set_bold();
	$header_format_1->set_right();
	$header_format_1->set_bottom();
	$header_format_1->set_top();
	$header_format_1->set_color( 'red' );
	$header_format_1->set_size(10);
    
	my $header_format_2 = $workbook->add_format(	
		font  => 'Arial'
	);
	$header_format_2->set_border();
	$header_format_2->set_size(10);
	my $worksheetDE = $workbook->add_worksheet();#"$excelname2"
	for (my $i=0;$i<=$#annInfoContent; $i++){
		
		chomp $annInfoContent[$i];
		$annInfoContent[$i]=~s/\s+$//;
		my @row_data;
		if($ensembl_file =~/csv/){
			@row_data=split /,/, $annInfoContent[$i];
		}else{
			@row_data=split /\t/, $annInfoContent[$i];
		}
		if($i==0){
			for (my $j=0;$j<=$#row_data;$j++){
				$worksheetDE->write($i, $j, $row_data[$j],$header_format_1);
			}
		}else{
			for (my $j=0;$j<=$#row_data;$j++){
				$worksheetDE->write($i, $j, $row_data[$j],$header_format_2);
				
			}
		}
	}
	$workbook->close();
}
return 1;
