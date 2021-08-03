sub enrich_ANA{
	#!/usr/bin/perl -w
	use strict;
	use Excel::Writer::XLSX;
	use Encode qw/encode decode/;
	use List::Util qw( min max );
	use Cwd;
	use Excel::Writer::XLSX;

	my $loc="/home/shihjielai/RNAseq_PIP/RNA_E-ReportData";
	chdir("R_dataFiles");
	####
	my @Rscript=<*.r>;
	foreach my $RR(@Rscript){
		open(HTML_Temp, "$loc/Enrichment_template_new.html") or die;
		my @HTML_Temp=<HTML_Temp>;close HTML_Temp;
	    print "manual";
	    system ("Rscript $RR");
	    my $E_Folder=$RR;
	    $E_Folder=~s/\.r//;$E_Folder="Enrich_".$E_Folder;
	    mkdir ("$E_Folder");
	    my @PNG1=<*.png>;
	    foreach my $png1(@PNG1){
			my $file=$E_Folder."\/".$png1;
			system ("mv $png1 $file");
	    }
	    my @PNG1=<*.svg>;
	    foreach my $png1(@PNG1){
			my $file=$E_Folder."\/".$png1;
			system ("mv $png1 $file");
	    }
	    my $E_htm=$E_Folder;$E_htm=~s/Enrich_/EnrichPlots_/;$E_htm=$E_htm.".html";
	    my $namename=$RR;$namename=~s/\.r//;
	    
	    open(E_HTML, ">$E_htm") or die;
		if($E_Folder!~/_GSEA/){
			foreach my $html_line (@HTML_Temp){
				$html_line=~s/FuckYouDamnit/$E_Folder/g;
				$html_line=~s/ILOVESPURS/$namename/g;
				print E_HTML $html_line;
			}
	    }
	    my $MF=$RR;$MF=~s/\.r/\.tab/;$MF="MF_".$MF;
	    my $BP=$RR;$BP=~s/\.r/\.tab/;$BP="BP_".$BP;
	    my $CC=$RR;$CC=~s/\.r/\.tab/;$CC="CC_".$CC;
		my $DO=$RR;$DO=~s/\.r/\.tab/;$DO="DO_".$DO;
	    my $KK=$RR;$KK=~s/\.r/\.tab/;my $KK_Foolder=$KK;$KK="KK_".$KK;$KK_Foolder="KEGG_".$KK;$KK_Foolder=~s/\.tab//;
	    my $Enrich=$RR;$Enrich=~s/\.r/\.xlsx/;$Enrich="Enrich_".$Enrich;
	    ####################################################################
		my $KK_pathview=$KK;$KK_pathview=~s/\.tab/\.pathview/;
		if($KK_pathview=~/GSEA/){$KK_pathview=~s/^KK/All/;$KK_pathview=~s/_GSEA//;}#
		my %EG2GS;
		open(EG2GS,"$KK_pathview") or die;
		while (<EG2GS>) {
			chomp $_;
			my @Temp=split /\t/, $_;
			$EG2GS{$Temp[0]}=$Temp[1]
		}
		####################################################################
		my %EG2GSS;
		open(EG2GSS,"b.tab") or die;
		while (<EG2GSS>) {
			chomp $_;
			my @Temp=split /\t/, $_;
			$EG2GSS{$Temp[0]}=$Temp[1]
		}
		####################################################################
	    my $a_0;$a_0="<a href=\"$Enrich\">$Enrich</a>";
	    my $workbook=Excel::Writer::XLSX->new("$Enrich");
	    my $header_format_1 = $workbook->add_format(
			                    align   => 'right',
		                        font  => 'Arial'
	    );
	    $header_format_1->set_bold();
	    $header_format_1->set_color('red');
	    $header_format_1-> set_bottom();
	    $header_format_1->set_size('10');
	    my $content_format_1 = $workbook->add_format(
								font  => 'Arial'
	    );    
	    my $html_name=$E_Folder.".html";
		#$html_name=~s/\.tab/_pathway.html/;
		my @html2_R;
	    
		if($RR!~/_GSEA/){
			open(html2, "$loc/Enrichment_html_new.html");@html2_R=<html2>;close html2;
			open(html1,">$html_name");
		}else{
			my $RR1=$RR;
			$RR1=~s/_GSEA//;$RR1=~s/\.r//;
			open(html2, "Enrich_$RR1.html");@html2_R=<html2>;close html2;
			open(html1,">Enrich_$RR1.html");
		}
	    open(CC, "$CC");my @CC_R=<CC>;close CC;
	    open(BP, "$BP");my @BP_R=<BP>;close BP;
	    open(MF, "$MF");my @MF_R=<MF>;close MF;
	    open(KK, "$KK");my @KK_R=<KK>;close KK;
	    my $worksheet_MF = $workbook->add_worksheet("Enrichment_MF");
		my $a="<tr>\n";
	    for(my $i=0;$i<=$#MF_R;$i++){
			chomp $MF_R[$i];$MF_R[$i]=~s/\"//g;
			my @Temp=split /\t/,$MF_R[$i];
			if ($i==0) {
				for(my $j=0;$j<=$#Temp;$j++){
					$worksheet_MF->write($i, $j, $Temp[$j],$header_format_1);
					$a=$a."<th>$Temp[$j]</th>\n";
				} 
				$a=$a."</tr>\n";
			}else{
				if($i<=20){
					$a=$a."<tr>\n";}
				for(my $j=0;$j<=$#Temp-1;$j++){
					my $b="";
					$worksheet_MF->write($i, $j, $Temp[$j+1],$content_format_1);
					if( $Temp[$j+1]=~/^-?[0-9]+\.[0-9]+e?/){
						$b=sprintf("%.4f", $Temp[$j+1]);}
					elsif($j>3 && $Temp[$j+1]=~/^[\w]+\/?[\w]+/){
						my @Temp1=split /\//, $Temp[$j+1];
						foreach my $PathwayAdded(@Temp1){
							$EG2GS{$EG2GSS{$PathwayAdded}}=sprintf("%.6f",$EG2GS{$EG2GSS{$PathwayAdded}});
							if ($EG2GS{$EG2GSS{$PathwayAdded}}>0) {
								$b=$b."<a title=\"$EG2GS{$EG2GSS{$PathwayAdded}}\" style=\"color:red;\">$PathwayAdded</a>/";
							}else{
								$b=$b."<a title=\"$EG2GS{$EG2GSS{$PathwayAdded}}\" style=\"color:blue;\">$PathwayAdded</a>/";
							}     
						}	
					}
					else {
						$b=$Temp[$j+1];}
						if($i<=20){
						$a=$a. "<td>".$b."</td>\n";
					}
				} 
				if($i<=20){
					$a=$a."</tr>";
				}	
			}      
	    }    
	    my $worksheet_BP = $workbook->add_worksheet("Enrichment_BP");
		my $a_1="<tr>\n";
	    for(my $i=0;$i<=$#BP_R;$i++){
		chomp $BP_R[$i];$BP_R[$i]=~s/\"//g;
		my @Temp=split /\t/,$BP_R[$i];
		if ($i==0) {
		    for(my $j=0;$j<=$#Temp;$j++){
		        $worksheet_BP->write($i, $j, $Temp[$j],$header_format_1);
				$a_1=$a_1."<th>$Temp[$j]</th>\n";
		    }
			$a_1=$a_1."</tr>\n";
		}else{
			if($i<=20){
				$a_1=$a_1."<tr>\n";}
		    for(my $j=0;$j<=$#Temp-1;$j++){
				my $b="";
		        $worksheet_BP->write($i, $j, $Temp[$j+1],$content_format_1);
				if( $Temp[$j+1]=~/^-?[0-9]+\.[0-9]+e?/){
					$b= sprintf("%.4f", $Temp[$j+1]);
				}elsif($j>3 && $Temp[$j+1]=~/^[\w]+\/?[\w]+/){
					my @Temp1=split /\//, $Temp[$j+1];
					foreach my $PathwayAdded(@Temp1){
		                $EG2GS{$EG2GSS{$PathwayAdded}}=sprintf("%.6f",$EG2GS{$EG2GSS{$PathwayAdded}});
						if ($EG2GS{$EG2GSS{$PathwayAdded}}>0) {
		                    $b=$b."<a title=\"$EG2GS{$EG2GSS{$PathwayAdded}}\" style=\"color:red;\">$PathwayAdded</a>/";
		                }else{
		                    $b=$b."<a title=\"$EG2GS{$EG2GSS{$PathwayAdded}}\" style=\"color:blue;\">$PathwayAdded</a>/";
		                }     
		            }	
				}else {
					$b= $Temp[$j+1];}
				if($i<=20){
					$a_1=$a_1. "<td>".$b."</td>\n";
				}
		    }
				if($i<=20){			
					$a_1=$a_1. "</tr>";
				}	  	
			}          
	    }    
	    my $worksheet_CC = $workbook->add_worksheet("Enrichment_CC");
		my $a_2="<tr>\n";
	    for(my $i=0;$i<=$#CC_R;$i++){
			chomp $CC_R[$i];$CC_R[$i]=~s/\"//g;
			my @Temp=split /\t/,$CC_R[$i];
			if ($i==0) {
				for(my $j=0;$j<=$#Temp;$j++){
					$worksheet_CC->write($i, $j, $Temp[$j],$header_format_1);
					$a_2=$a_2."<th>$Temp[$j]</th>\n";
				}
				$a_2=$a_2."</tr>\n";
			}else{
				if($i<=20){
					$a_2=$a_2."<tr>\n";}
				for(my $j=0;$j<=$#Temp-1;$j++){
					my $b="";
					$worksheet_CC->write($i, $j, $Temp[$j+1],$content_format_1);
					if( $Temp[$j+1]=~/^-?[0-9]+\.[0-9]+e?/){
						$b= sprintf("%.4f", $Temp[$j+1]);}
					elsif($j>3 && $Temp[$j+1]=~/^[\w]+\/?[\w]+/){
						my @Temp1=split /\//, $Temp[$j+1];
						foreach my $PathwayAdded(@Temp1){
							$EG2GS{$EG2GSS{$PathwayAdded}}=sprintf("%.6f",$EG2GS{$EG2GSS{$PathwayAdded}});
							if ($EG2GS{$EG2GSS{$PathwayAdded}}>0) {
								$b=$b."<a title=\"$EG2GS{$EG2GSS{$PathwayAdded}}\" style=\"color:red;\">$PathwayAdded</a>/";
							}else{
								$b=$b."<a title=\"$EG2GS{$EG2GSS{$PathwayAdded}}\" style=\"color:blue;\">$PathwayAdded</a>/";
							}     
						}	
					}
					else {
						$b= $Temp[$j+1];}
					if($i<20){
						$a_2=$a_2."<td>".$b."</td>\n";
					}
				} 
				if($i<=20){
					$a_2=$a_2."</tr>";
				}	
			}                        
	    }    

	    my $worksheet_KEGG = $workbook->add_worksheet("Enrichment_KEGG");
		my $a_3="<tr>\n";
	    my $url_format = $workbook->add_format(
							color     => 'blue',
							underline => 1,
							font  => 'Arial',
	    );
	    my $PTotal_name=$KK;
	    $PTotal_name=~s/\.tab/_pathway_total.txt/;
	    open(PTotal,">$PTotal_name");
	    my $KEGG_Org;
	    for(my $i=0;$i<=$#KK_R;$i++){
			chomp $KK_R[$i];$KK_R[$i]=~s/\"//g;
			my @Temp=split /\t/,$KK_R[$i];
			if ($i==0){
				for(my $j=0;$j<=$#Temp;$j++){
					$worksheet_KEGG->write($i, $j, $Temp[$j],$header_format_1);
					$a_3=$a_3."<th>$Temp[$j]</th>\n";
				}
				$a_3=$a_3."</tr>\n";
			}else{
				if($i<=20){
					$a_3=$a_3."<tr>\n";
				}
				for(my $j=0;$j<=$#Temp-1;$j++){
					my $b="";
					if ($j==0) {
						my $Flie_location=$KK_Foolder."\/".$Temp[$j+1].".welgene.png";
						$worksheet_KEGG->write_url($i, $j, "$Flie_location",$url_format,$Temp[$j],$Temp[$j+1]);
						$KEGG_Org=substr $Temp[$j+1],0,3;
						if($i<20){
							$a_3=$a_3."<td><a href=\".\/".$KK_Foolder."\/".$Temp[$j+1].".welgene.png\" target=\"_parent\" title=".$Temp[$j+1].">".$Temp[$j+1]."</a></td>\n";}
						print PTotal "$Temp[$j+1]\n";
					}else{
						$worksheet_KEGG->write($i, $j, $Temp[$j+1],$content_format_1);
						if( $Temp[$j+1]=~/^-?[0-9]+\.[0-9]+e?/){
							$b= sprintf("%.4f", $Temp[$j+1]);}
						elsif($j>3 && $Temp[$j+1]=~/^[\w]+\/?[\w]+/){
							my @Temp1=split /\//, $Temp[$j+1];
							foreach my $PathwayAdded(@Temp1){
								$EG2GS{$EG2GSS{$PathwayAdded}}=sprintf("%.6f",$EG2GS{$EG2GSS{$PathwayAdded}});
								if ($EG2GS{$EG2GSS{$PathwayAdded}}>0) {
									$b=$b."<a title=\"$EG2GS{$EG2GSS{$PathwayAdded}}\" style=\"color:red;\">$PathwayAdded</a>/";
								}else{
									$b=$b."<a title=\"$EG2GS{$EG2GSS{$PathwayAdded}}\" style=\"color:blue;\">$PathwayAdded</a>/";
								}     
							}	
						}
						else {
							$b= $Temp[$j+1];
						}
						if($i<20){
							$a_3=$a_3."<td>".$b."</td>\n";
						}
					}   
				}   
				if($i<=20){
					$a_3=$a_3."</tr>";
				}	
			}                        
	    }
		my $a_4="<!---->";
		if(-e $DO){
			open(DO, "$DO");my @CC_R=<DO>;close DO;
			my $worksheet_DO = $workbook->add_worksheet("Enrichment_DO");
			$a_4="<tr>\n";
			for(my $i=0;$i<=$#CC_R;$i++)	{
				chomp $CC_R[$i];$CC_R[$i]=~s/\"//g;
				my @Temp=split /\t/,$CC_R[$i];
				if ($i==0) {
					for(my $j=0;$j<=$#Temp;$j++){
						$worksheet_DO->write($i, $j, $Temp[$j],$header_format_1);
						$a_4=$a_4."<th>$Temp[$j]</th>\n";
					} 
					$a_4=$a_4."</tr>\n";
				}else{
					if($i<=20){
						$a_4=$a_4."<tr>\n";}
					for(my $j=0;$j<=$#Temp-1;$j++){
						my $b="";
						$worksheet_DO->write($i, $j, $Temp[$j+1],$content_format_1);
						if( $Temp[$j+1]=~/^-?[0-9]+\.[0-9]+e?/){
							$b=sprintf("%.4f", $Temp[$j+1]);}
						elsif($j>3 && $Temp[$j+1]=~/^[\w]+\/?[\w]+/){
							my @Temp1=split /\//, $Temp[$j+1];
							foreach my $PathwayAdded(@Temp1){
								$EG2GS{$EG2GSS{$PathwayAdded}}=sprintf("%.6f",$EG2GS{$EG2GSS{$PathwayAdded}});
								if ($EG2GS{$EG2GSS{$PathwayAdded}}>0) {
									$b=$b."<a title=\"$EG2GS{$EG2GSS{$PathwayAdded}}\" style=\"color:red;\">$PathwayAdded</a>/";
								}else{
									$b=$b."<a title=\"$EG2GS{$EG2GSS{$PathwayAdded}}\" style=\"color:blue;\">$PathwayAdded</a>/";
								}     
							}	
						}
						else {
							$b=$Temp[$j+1];}
						if($i<=20){
							$a_4=$a_4. "<td>".$b."</td>\n";
						}
					} 
					if($i<=20){
						$a_4=$a_4."</tr>";
					}	
				}      
			}
		}	
		foreach my $html_line1 (@html2_R){
			if($RR!~/GSEA/){
				$html_line1=~s/All_DATA/$a_0/g;
				$html_line1=~s/MF_DATA/$a/g;
				$html_line1=~s/BP_DATA/$a_1/g;
				$html_line1=~s/CC_DATA/$a_2/g;
				$html_line1=~s/KK_DATA/$a_3/g;
				$html_line1=~s/DO_DATA/$a_4/g;
			}else{
				$html_line1=~s/All_G_DATA/$a_0/g;
				$html_line1=~s/MF_G_DATA/$a/g;
				$html_line1=~s/BP_G_DATA/$a_1/g;
				$html_line1=~s/CC_G_DATA/$a_2/g;
				$html_line1=~s/KK_G_DATA/$a_3/g;
				$html_line1=~s/DO_G_DATA/$a_4/g;
			}
			print html1 $html_line1;
	    }
		close html1;
	    close PTotal;
		
	    my $KK_R=$KK_pathview;$KK_R=~s/\.pathview/\.r/;
	    open(Pathway_List, "$PTotal_name") or die;
	    open(RSCRIPT,">$KK_R") or die;
	    print RSCRIPT "library(pathview)\ngene.data=read.table(\"$KK_pathview\",sep=\"\\t\",row.names=1)\n";
	
		my $co=0;
	    while (<Pathway_List>) {
			$co=$co+1;
			if($co>22){last;}
			chomp $_;
			print RSCRIPT "tryCatch({pv.out <- pathview(gene.data, pathway.id = \"$_\", species = \"$KEGG_Org\", out.suffix = \"welgene\",limit=list(gene=3,cpd=1))},error=function(e){})\n";        
	    }
	    close RSCRIPT;
	    system("Rscript $KK_R");
	    my @png=<*.png>;
	    mkdir ("$KK_Foolder");
	    foreach my $plot(@png){
		    system ("mv $plot $KK_Foolder\/$plot");
	    }
		my @tab=<*.tab>;
		foreach my $TT(@tab){
			unlink $TT;
		}
		my @tab1=<KK*.r>;
		foreach my $TT(@tab1){
			unlink $TT;
		}
	
	}
	my @xxx=<*.xml>;
	foreach my $xxl(@xxx){
		unlink $xxl;
	}
	chdir("..");}
return 1;
