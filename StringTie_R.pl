sub StringTie_R_DESeq2 { # Run DEseq2
	($ab,$FC,$rawread,$rawdata,$gene,$PValue,%hash_s)= @_;
	#ab: groupCompareInfo；FC:FC :  rawread: raw read data； rawdata: FPKM or TPM；gene: TS or Gene ；hash_s :sample2GroupInfo
	chomp $ab;
	my @Temp2=split /\t/, $ab;
	if(!$hash_s{$Temp2[0]}){push @{$hash_s{$Temp2[0]}}, $Temp2[0];};
	if(!$hash_s{$Temp2[1]}){push @{$hash_s{$Temp2[1]}}, $Temp2[1];};
	my	$treat = join ',' , @{$hash_s{$Temp2[0]}};$treat=~s/,/\",\"/g;$treat=~s/-/\./g;
	my	$cont = join ',' , @{$hash_s{$Temp2[1]}};$cont=~s/,/\",\"/g;$cont=~s/-/\./g;
	my $nameR;
	if ($gene==0){
		$nameR="$Temp2[0]_vs_$Temp2[1]";
	}else{
		$nameR="TS_$Temp2[0]_vs_$Temp2[1]";
	}
	
	my $p=$PValue;

	
	open(RSCRIPT,">$nameR.R") or die;
	
		print RSCRIPT
		"readdata = read.csv(\"$rawread\", row.names=1)
		rawdata <- read.delim(\"$rawdata\", row.names=1)
		
		head( readdata )
		library(DESeq2)
		library(ggplot2)
		library(EnhancedVolcano)
		
		
		datause<-readdata[,c(\"$treat\",\"$cont\")]
		group = c( rep(\"treated\",$#{$hash_s{$Temp2[0]}}+1),rep(\"control\",$#{$hash_s{$Temp2[1]}}+1))
		pasillaDesign = data.frame(row.names = colnames( datause ),group)
		dds <- DESeqDataSetFromMatrix(countData = datause, colData = pasillaDesign, design = ~ group)
		dds <- DESeq(dds, minReplicatesForReplace=Inf)
		res <- results(dds, cooksCutoff=FALSE,independentFiltering=FALSE)
		datares<-as.data.frame(res)
		ansd<-datares[,c(5,6)]
		#ansd<-round(ansd,4)
		write.table(ansd,\"${nameR}.tab\", sep = \"\t\",na = \"1\",col.names = F, row.names = T, quote = F)
		
		datause2<-rawdata[,c(\"$treat\",\"$cont\")]
		mena_treat_nor<-rowMeans(as.matrix(datause2[,which(group==\"treated\")]))
		mena_con_nor<-rowMeans(as.matrix(datause2[,which(group==\"control\")]))
		
		sd_treat_nor<-rowSds(as.matrix(datause2[,which(group==\"treated\")]))
		sd_con_nor<-rowSds(as.matrix(datause2[,which(group==\"control\")]))
		
		FC<-mena_treat_nor/mena_con_nor
		log2FC<-log2(FC)
		datause3<-data.frame(datause2[,which(group==\"treated\")],mena_treat_nor,sd_treat_nor,datause2[,which(group==\"control\")],mena_con_nor,sd_con_nor,FC,log2FC)
		#datause3<-round(datause3,4)

		datause4<-merge(datause3, ansd, by = \"row.names\")
		datause4[is.na(datause4)]<-1
		datause4\$aa<-ifelse(datause4\$padj<0.05,\"yes\",\"no\")
		

		

		################
		datause6<-datause4[which((datause4\$mena_treat_nor>0.3 | datause4\$mena_con_nor>0.3) ),]#| datause4\$padj<0.05
		
		keyvals <- ifelse(
		  datause6\$log2FC < -log2($FC) & datause6\$pvalue<=$p, 'royalblue',
		  ifelse(datause6\$log2FC > log2($FC) & datause6\$pvalue<=$p, 'red',
		         'black'))
		keyvals[is.na(keyvals)] <- 'black'
		names(keyvals)[keyvals == 'red'] <- 'Up signification'
		names(keyvals)[keyvals == 'black'] <- 'No signification'
		names(keyvals)[keyvals == 'royalblue'] <- 'Down signification'
		family<-c(length(which(keyvals == 'red')),length(which(keyvals == 'royalblue')))
		datause6\$pvalue[is.na(datause6\$pvalue)]<-1
		datause6\$pvalue[is.na(datause6\$pvalue)]<-1
		
		df <- data.frame(x = c(max(datause6\$log2FC)*(2/3),min(datause6\$log2FC)*(2/3)), y = quantile(seq(0,-log10(min(datause6\$pvalue[datause6\$pvalue>0])),1),0.7), family = family,col=c(\"red\",\"royalblue\"))

		mim_p<- -log10(min(datause6\$pvalue[datause6\$pvalue>0]))
		if(mim_p\%\/\%5<1){
			min_p2<-1
			mim_p<-4
		}else{
			min_p2<-mim_p\%\/\%5
		}
		png(filename = \"${nameR}_Volcano.png\",width = 3000, height = 2000,res=300)
		g<-EnhancedVolcano(datause6,lab = \"\" ,selectLab = rownames(datause6)[which(names(keyvals) \%in\% c('high', 'low'))],FCcutoff=log2($FC), colCustom = keyvals,x = \"log2FC\",y = \"pvalue\",pCutoff =$p)+ scale_y_continuous(trans = \"log1p\",breaks = round(c(0, -log10($p),seq(2,mim_p,min_p2)),2))+geom_text(aes(x, y, label = family),data =df,size = 8,color=c(\"red\",\"royalblue\"))
		print(g)
		dev.off()	
		
		################

		colnames(datause4)<-c(\"gene_id\",\"$treat\",\"Mean_$Temp2[0]\",\"SD_$Temp2[0]\",\"$cont\",\"Mean_$Temp2[1]\",\"SD_$Temp2[1]\",\"ratio($Temp2[0]\/$Temp2[1])\",\"log2ratio($Temp2[0]\/$Temp2[1])\",\"p_value\",\"q_value\",\"significant\")

		datause5<-datause4[which(datause4\$p_value<=$p & ( datause4\$\"log2ratio($Temp2[0]\/$Temp2[1])\" > log2($FC) | datause4\$\"log2ratio($Temp2[0]\/$Temp2[1])\" < -log2($FC) ) & (datause4\$Mean_$Temp2[0]>0.3 | datause4\$Mean_$Temp2[1]>0.3)),]

		write.table(datause4,\"${nameR}_DE1.tab\", sep = \"\t\",na = \"1\",col.names = T, row.names = F, quote = F)
		write.table(datause5,\"${nameR}_FC${FC}X1.tab\", sep = \"\t\",na = \"1\",col.names = T, row.names = F, quote = F)
		
		
		png(filename = \"${nameR}_MA.png\",width = 6000, height = 6000,res=1000)
		g<-ggplot(datause4,aes(x=(Mean_$Temp2[0]+Mean_$Temp2[1])/2, y=log2(Mean_$Temp2[0]/Mean_$Temp2[1]), col=ifelse(p_value<=$p,\"yes\",\"no\"))) + geom_point() + scale_x_log10() + ggtitle(\"MA plot\")+ labs(col = \"p_value(<=$p)\")
		print(g)
		dev.off()
		
		png(filename = \"${nameR}_scatter.png\",width = 6000, height = 6000,res=1000)
		g<-ggplot(datause4,aes(x=Mean_$Temp2[0], y=Mean_$Temp2[1], col=ifelse(p_value<=$p,\"yes\",\"no\"))) + geom_point() + scale_x_log10(limits = c(0.0001,100000)) + scale_y_log10(limits = c(0.0001,100000))+ ggtitle(\"$Temp2[0]_vs_$Temp2[1]\")+ labs(col = \"p_value(<=$p)\")
		print(g)
		dev.off()
		
		
		";
	close RSCRIPT;
	system("Rscript $nameR.R");
	return $nameR;
}

sub StringTie_R_DESeq { # Run DEseq
	($ab,$FC,$rawread,$rawdata,$gene,$PValue)= @_;
	#ab: groupCompareInfo；FC:FC :  rawread: raw read data； rawdata: FPKM or TPM；gene: TS or Gene
	chomp $ab;
	my @Temp2=split /\t/, $ab;
	my	$treat =$Temp2[0];
	my	$cont = $Temp2[1];
	my $nameR;
	if ($gene==0){
		$nameR="$Temp2[0]_vs_$Temp2[1]";
	}else{
		$nameR="TS_$Temp2[0]_vs_$Temp2[1]";
	}
	
	my $p=$PValue;

	open(RSCRIPT,">$nameR.R") or die;
	
		print RSCRIPT
		"
		pasillaCountTable = read.csv(\"$rawread\", row.names=1)
		rawdata <- read.delim(\"$rawdata\", row.names=1)
		head( pasillaCountTable )	
	
		library( \"DESeq\" )
		library(ggplot2)
		library(EnhancedVolcano)
		
		pasillaCountTable1=pasillaCountTable
		pasillaDesign = data.frame(
		row.names = colnames(pasillaCountTable1),
		condition = colnames(pasillaCountTable1))
		pasillaDesign
		countTable = pasillaCountTable1
		condition = pasillaDesign\$condition
		cds = newCountDataSet( countTable, condition )
		cds = estimateSizeFactors( cds )
		
		cds3 =tryCatch({estimateDispersions( cds, method=\"blind\", sharingMode=\"fit-only\" )},error=function(e){estimateDispersions( cds, method=\"blind\", sharingMode=\"fit-only\",fitType=\"local\" )})
		#cds3 = estimateDispersions( cds, method=\"blind\", sharingMode=\"fit-only\" )
		res2 = nbinomTest( cds3, \"$Temp2[1]\", \"$Temp2[0]\" )
		ansd<-data.frame(res2\$id,res2\$pval,res2\$padj)
		ansd<-ansd[order(ansd[,1]),]
		row.names(ansd)<-ansd[,1]
		ansd<-ansd[,-1]
		write.table(ansd,\"$nameR.tab\", sep = \"\t\",na = \"1\",col.names = F, row.names = F, quote = F)
		datause2<-rawdata[,c(\"$treat\",\"$cont\")]
		FC<-datause2\$$treat/datause2\$$cont
		log2FC<-log2(FC)
		datause3<-data.frame(datause2[,c(\"$treat\")],datause2[,c(\"$cont\")],FC,log2FC)
		row.names(datause3)<-row.names(datause2)
		datause4<-merge(datause3, ansd, by = \"row.names\")
		datause4[is.na(datause4)]<-1
		datause4\$aa<-ifelse(datause4\$res2.padj<0.05,\"yes\",\"no\")
		colnames(datause4)<-c(\"gene_id\",\"$treat\",\"$cont\",\"ratio($Temp2[0]\/$Temp2[1])\",\"log2ratio($Temp2[0]\/$Temp2[1])\",\"p_value\",\"q_value\",\"significant\")
		
		
		################
		datause6<-datause4[which((datause4\$$treat>0.3 | datause4\$$cont>0.3)),]
		
		keyvals <- ifelse(
		  datause6\$\"log2ratio($Temp2[0]\/$Temp2[1])\" < -log2($FC) & datause6\$p_value<=$p, 'royalblue',
		  ifelse(datause6\$\"log2ratio($Temp2[0]\/$Temp2[1])\" > log2($FC) & datause6\$p_value<=$p, 'red',
		         'black'))
		keyvals[is.na(keyvals)] <- 'black'
		names(keyvals)[keyvals == 'red'] <- 'Up signification'
		names(keyvals)[keyvals == 'black'] <- 'No signification'
		names(keyvals)[keyvals == 'royalblue'] <- 'Down signification'
		family<-c(length(which(keyvals == 'red')),length(which(keyvals == 'royalblue')))
		datause6\$p_value[is.na(datause6\$p_value)]<-1
		
		mim_p<- -log10(min(datause6\$p_value[datause6\$p_value>0]))		
		if(mim_p\%\/\%5<1){
			min_p2<-1
			mim_p<-4
		}else{
			min_p2<-mim_p\%\/\%5
		}

		df <- data.frame(x = c(max(datause6\$\"log2ratio($Temp2[0]\/$Temp2[1])\")*(2/3),min(datause6\$\"log2ratio($Temp2[0]\/$Temp2[1])\")*(2/3)), y = quantile(seq(0,-log10(min(datause6\$p_value[datause6\$p_value>0])),1),0.7), family = family,col=c(\"red\",\"royalblue\"))
		png(filename = \"${nameR}_Volcano.png\",width = 3000, height = 2000,res=300)
		g<-EnhancedVolcano(datause6,lab = \"\" ,selectLab = rownames(datause6)[which(names(keyvals) \%in\% c('high', 'low'))],FCcutoff=log2($FC), colCustom = keyvals,x = \"log2ratio($Temp2[0]\/$Temp2[1])\",y = \"p_value\",pCutoff =$p)+ scale_y_continuous(trans = \"log1p\",breaks = round(c(0, -log10($p),seq(2,mim_p,min_p2)),2))+geom_text(aes(x, y, label = family),data =df,size = 8,color=c(\"red\",\"royalblue\"))
		print(g)
		dev.off()	
		
		################
		
		datause5<-datause4[which(datause4\$p_value<=$p & ( datause4\$\"log2ratio($Temp2[0]\/$Temp2[1])\" > log2($FC) | datause4\$\"log2ratio($Temp2[0]\/$Temp2[1])\" < -log2($FC) ) & (datause4\$$treat>0.3 | datause4\$$cont>0.3)),]
		write.table(datause4,\"${nameR}_DE1.tab\", sep = \"\t\",na = \"1\",col.names = T, row.names = F, quote = F)
		write.table(datause5,\"${nameR}_FC${FC}X1.tab\", sep = \"\t\",na = \"1\",col.names = T, row.names = F, quote = F)
		
		png(filename = \"${nameR}_MA.png\",width = 6000, height = 6000,res=1000)
		g<-ggplot(datause4,aes(x=($treat+$cont)/2, y=log2($treat/$cont), col=ifelse(p_value<=$p,\"yes\",\"no\"))) + geom_point() + scale_x_log10() + ggtitle(\"MA plot\")+ labs(col = \"p_value(<=$p)\")
		print(g)
		dev.off()
		
		png(filename = \"${nameR}_scatter.png\",width = 6000, height = 6000,res=1000)
		g<-ggplot(datause4,aes(x=$treat, y=$cont, col=ifelse(p_value<=$p,\"yes\",\"no\"))) + geom_point() + scale_x_log10(limits = c(0.0001,100000)) + scale_y_log10(limits = c(0.0001,100000))+ ggtitle(\"${treat}_vs_${cont}\")+ labs(col = \"p_value(<=$p)\")
		print(g)
		dev.off()
		
		";
	close RSCRIPT;
	system("Rscript $nameR.R");
	return $nameR;
}

sub StringTie_R_plot { # draw plot 
	($exp)= @_;
	my $loc_R="cummer";
	if(! -d "cummer"){
		mkdir ("cummer");
	}
	my $exp1;
	if($exp==1){
	print "FPKM";
		$exp1="FPKM";
	}else{
	print "TPM";
		$exp1="TPM";
	}
	
	open(RSCRIPT,">Plot.R") or die;
	print RSCRIPT "
	library(\"gplots\")
	library(\"ggplot2\")
	library(\"philentropy\")
	library(\"GGally\")
	library(\"plyr\")
	library(\"factoextra\")
	 
	readdata = read.csv(\"gene_count_matrix.csv\", row.names=1)
	rawdata <- read.delim(\"rawdata.txt\", row.names=1)
	####corr####
	rawdata_C <- rawdata[,c(-1,-2)]
	rawdata_M <- merge(readdata,rawdata_C,by=\"row.names\")
	rawdata_M <- rawdata_M[,-1]
	col_n <- colnames(rawdata_M)
	b <- c(\"X\",\"Y\",\"rho\")
	for (i in 1:(length(col_n)/2)){
		for(j in ((length(col_n)/2)+1):length(col_n)){
			pearson_R <- cor.test(rawdata_M[,i],rawdata_M[,j],method=\"spearman\")
			a <- c(col_n[i],col_n[j],pearson_R\$estimate)
			b <- rbind(b,a)
		}
	}
	write.csv(b,\"spearman.csv\")
	####BOX####
	cc<-as.list(rawdata[,c(-1,-2)])
	df <- ldply (cc, data.frame)
	cc1<-df[which(df[,2]!=0.0001),]
	colnames(cc1)<-c(\"condition\",\"$exp1\")
	p <- ggplot(cc1, aes(x=condition, y=log10($exp1)),color=condition) + 
	  geom_boxplot(aes(fill=condition)) + 
	  theme(axis.text.x=element_text(angle=50,hjust=0.5, vjust=0.5)) 
	png(filename = \"box.png\",width = 6000, height = 6000,res=1000)
	p
	dev.off()
	####dens####
	p <- ggplot(cc1, aes(x=log10($exp1)),color=condition) + 
	  geom_density(aes(fill=condition),alpha = 1/2) + 
	  theme(axis.text.x=element_text(angle=50,hjust=0.5, vjust=0.5))
	png(filename = \"dens.png\",width = 6000, height = 6000,res=1000)
	p
	dev.off()
	####dens####
	cc<-as.list(readdata)
	df <- ldply (cc, data.frame)
	colnames(df)<-c(\"condition\",\"counts\")
	df\$counts<-df\$counts+1
	p <- ggplot(df, aes(x=log2(counts)),color=condition) + 
	  geom_density(aes(fill=condition),alpha = 1/2) + 
	  theme(axis.text.x=element_text(angle=50,hjust=0.5, vjust=0.5))+ labs(x = \"log2(counts+1)\")
	png(filename = \"dens_count.png\",width = 6000, height = 6000,res=1000)
	p
	dev.off()	
	####JSD####
	dd<-t(as.matrix(rawdata[,c(-1,-2)]))
	row.names(dd)
	dd1<-JSD(dd,est.prob = \"empirical\")
	
	if(length(dd1)==1){
	 dd2<-matrix(c(0,dd1,dd1,0),2,2)
	}else{dd2<-dd1}
	
	row.names(dd2)<-row.names(dd)
	colnames(dd2)<-row.names(dd)
	png(filename = \"DistHeat.png\",width = 6000, height = 6000,res=1000)
	aa<-heatmap.2(dd2, scale=\"none\",key=TRUE, symkey=FALSE, density.info=\"none\", trace=\"none\", cexRow=0.5,cexCol=0.7)
	dev.off()
	
	####Scatter####
	aa<-as.matrix(rawdata[,c(-1,-2)])
	aa[which(aa==0.0001)]<-0
	cc<-as.data.frame(log10(aa))
	upperfun <- function(data,mapping){
	  ggplot(data = data, mapping = mapping)+
		geom_point()+
		scale_x_continuous(limits = c(-0.5,4))+
		scale_y_continuous(limits = c(-0.5,4))
	}   
	lowerfun <- function(data,mapping){
	  ggplot(data = data, mapping = mapping)+
		geom_point()+
		scale_x_continuous(limits = c(-0.5,4))+
		scale_y_continuous(limits =c(-0.5,4))
	} 
	png(filename = \"csScatterMatrix.png\",width = 6000, height = 6000,res=800)
	g<-ggpairs(cc,lower = list(continuous = wrap(lowerfun)),upper = list(continuous = wrap(upperfun)))
	print(g)
	dev.off()
	####PCA####
	aa<-rawdata[,c(-1,-2)]
	res.pca <- prcomp(aa, scale = TRUE)
	png(filename = \"pca.png\",width = 6000, height = 6000,res=1000)
	g<-fviz_pca_biplot(res.pca,invisible=\"ind\",alpha.var=1,geom.var=c(\"point\",\"text\"),repel = TRUE)
	print(g)
	dev.off()";
	
	close RSCRIPT;
	system("Rscript Plot.R");
	my @png=<*.png>;
	foreach my $plot(@png){
            system ("mv $plot $loc_R/$plot");
    }
}

sub StringTie_R_RE {
	($ab,$FC,$gene,$PValue,%hash_s)= @_;
	#ab: groupCompareInfo；FC:FC :  rawread: raw read data； rawdata: FPKM or TPM；gene: TS or Gene ；hash_s :sample2GroupInfo
	chomp $ab;
	my @Temp2=split /\t/, $ab;
	my $treat=$Temp2[0];
	my $cont=$Temp2[1];
	my $rowname="colnames(datause4)<-c(\"gene_id\",\"$treat\",\"$cont\",\"ratio($Temp2[0]\/$Temp2[1])\",\"log2ratio($Temp2[0]\/$Temp2[1])\",\"p_value\",\"q_value\",\"significant\")";
	if ($hash_s{$Temp2[0]} | $hash_s{$Temp2[1]}) { 
		$treat = join ',' , @{$hash_s{$Temp2[0]}};$treat=~s/,/\",\"/g;$treat=~s/-/\./g;
		$cont = join ',' , @{$hash_s{$Temp2[1]}};$cont=~s/,/\",\"/g;$cont=~s/-/\./g;
		$rowname="colnames(datause4)<-c(\"gene_id\",\"$treat\",\"$Temp2[0]\",\"SD_$Temp2[0]\",\"$cont\",\"$Temp2[1]\",\"SD_$Temp2[1]\",\"ratio($Temp2[0]\/$Temp2[1])\",\"log2ratio($Temp2[0]\/$Temp2[1])\",\"p_value\",\"q_value\",\"significant\")";
	}
	my $nameR;
	if ($gene==0){
		$nameR="$Temp2[0]_vs_$Temp2[1]";
	}else{
		$nameR="TS_$Temp2[0]_vs_$Temp2[1]";
	}
	my $p=$PValue;

	open(RSCRIPT,">${nameR}_re.R") or die;
		print RSCRIPT
		"
		library(ggplot2)
		library(EnhancedVolcano)
		datause4<- read.delim(\"${nameR}_DE1.tab\")
		$rowname
		datause4\$significant<-ifelse(datause4\$q_value<0.05,\"yes\",\"no\")	
		datause5<-datause4[which(datause4\$p_value<=$p & (datause4\$\"log2ratio($Temp2[0]\/$Temp2[1])\" <=-log2($FC) | log2($FC)<=datause4\$\"log2ratio($Temp2[0]\/$Temp2[1])\") & (datause4\$\"$Temp2[0]\">0.3 | datause4\$\"$Temp2[1]\">0.3)),]
		write.table(datause5,\"${nameR}_FC${FC}X1.tab\", sep = \"\t\",na = \"1\",col.names = T, row.names = F, quote = F)
		write.table(datause4,\"${nameR}_DE1.tab\", sep = \"\t\",na = \"1\",col.names = T, row.names = F, quote = F)
		################
		datause6<-datause4[which((datause4\$\"$Temp2[0]\">0.3 | datause4\$\"$Temp2[1]\">0.3)),]

		keyvals <- ifelse(
		  datause6\$\"log2ratio($Temp2[0]\/$Temp2[1])\" < -log2($FC) & datause6\$p_value<=$p, 'royalblue',
		  ifelse(datause6\$\"log2ratio($Temp2[0]\/$Temp2[1])\" > log2($FC) & datause6\$p_value<=$p, 'red',
		         'black'))
		keyvals[is.na(keyvals)] <- 'black'
		names(keyvals)[keyvals == 'red'] <- 'Up signification'
		names(keyvals)[keyvals == 'black'] <- 'No signification'
		names(keyvals)[keyvals == 'royalblue'] <- 'Down signification'
		family<-c(length(which(keyvals == 'red')),length(which(keyvals == 'royalblue')))
		datause6\$p_value[is.na(datause6\$p_value)]<-1
		df <- data.frame(x = c(max(datause6\$\"log2ratio($Temp2[0]\/$Temp2[1])\")*(2/3),min(datause6\$\"log2ratio($Temp2[0]\/$Temp2[1])\")*(2/3)), y = quantile(seq(0,-log10(min(datause6\$p_value[datause6\$p_value>0])),1),0.7), family = family,col=c(\"red\",\"royalblue\"))
		
		mim_p<- -log10(min(datause6\$p_value[datause6\$p_value>0]))
		if(mim_p\%\/\%5<1){
			min_p2<-1
			mim_p<-4
		}else{
			min_p2<-mim_p\%\/\%5
		}
		png(filename = \"${nameR}_Volcano.png\",width = 3000, height = 2000,res=300)
		g<-EnhancedVolcano(datause6,lab = \"\" ,selectLab = rownames(datause6)[which(names(keyvals) \%in\% c('high', 'low'))],FCcutoff=log2($FC), colCustom = keyvals,x = \"log2ratio($Temp2[0]\/$Temp2[1])\",y = \"p_value\",pCutoff =$p)+ scale_y_continuous(trans = \"log1p\",breaks = round(c(0, -log10($p),seq(2,mim_p,min_p2)),2))+geom_text(aes(x, y, label = family),data =df,size = 8,color=c(\"red\",\"royalblue\"))
		print(g)
		dev.off()	
		
		################		
		
		png(filename = \"${nameR}_MA.png\",width = 6000, height = 6000,res=1000)
		g<-ggplot(datause4,aes(x=(`$Temp2[0]`+`$Temp2[1]`)/2, y=log2(`$Temp2[0]`/`$Temp2[1]`), col=ifelse(p_value<=$p,\"yes\",\"no\"))) + geom_point() + scale_x_log10() + ggtitle(\"MA plot\")+ labs(col = \"p_value(<=$p)\")
		print(g)
		dev.off()
		png(filename = \"${nameR}_scatter.png\",width = 6000, height = 6000,res=1000)
		g<-ggplot(datause4,aes(x=`$Temp2[0]`, y=`$Temp2[1]`, col=ifelse(p_value<=$p,\"yes\",\"no\"))) + geom_point() + scale_x_log10(limits = c(0.0001,100000)) + scale_y_log10(limits = c(0.0001,100000))+ ggtitle(\"$Temp2[0]_vs_$Temp2[1]\")+ labs(col = \"p_value(<=$p)\")
		print(g)
		dev.off()
		
		";
	close RSCRIPT;
	system("Rscript ${nameR}_re.R");
	return $nameR;
}
return 1;

