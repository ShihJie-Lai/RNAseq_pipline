sub Stringtie_gene{ # combine all genes of samples
($select_E)= @_;
# select_E: RPKM or TPM
if($select_E==1){
	print "FPKM\n";
	$select_E=7;
}else{
	print "TPM\n";
	$select_E=8;
}
my $dir1 = getcwd;
my $dir = "ballgown";
opendir DIR, $dir;
my @etc_all_file = readdir DIR;
closedir DIR;
@etc_all_file = grep {!/\./} @etc_all_file;
my %hash_a;
my %hash_v=();
open (ans2,">rawdata.txt");
my $fi;
foreach my $aa(@etc_all_file){
	open(IN, "<", "$dir1/ballgown/${aa}/${aa}_gene_abund.tab") or die "Could not open \n\n";
	my @IN=<IN>;close (IN);
	foreach my $le(@IN){
		chomp $le;
		my @Temp1=split /\t/, $le;
		if ($Temp1[0]=~/^Gene\s/){
			$fi=$Temp1[0];
			$hash_a{$Temp1[0]}="$Temp1[1]\tlocus";
		}else{
			my $ff;
			if($Temp1[0]=~/./){
				$ff=$Temp1[1]
			}else{
				$ff=$Temp1[0]
			}
			$hash_a{$ff}=();
			push @{$hash_a{$ff}},"$Temp1[1]";
			push @{$hash_a{$ff}},"$Temp1[2]:$Temp1[4]-$Temp1[5]($Temp1[3])";
			######## if FPKM or TPM is 0 , the gene is 0.0001 #########
			if($Temp1[$select_E]<=0.000000001){
				######## if the gene have more than 2 vale , combine them #######
				if($hash_v{$aa}{$ff}){
					if($hash_v{$aa}{$ff}==0.0001){
						$hash_v{$aa}{$ff}=0.0001;
					}else{
						$hash_v{$aa}{$ff}=$hash_v{$aa}{$ff}+$Temp1[$select_E];
					}
				}else{
					$hash_v{$aa}{$ff}=0.0001;
				}
				
			}else{
				if($hash_v{$aa}{$ff}){
					if($hash_v{$aa}{$ff}==0.0001){
						$hash_v{$aa}{$ff}=$Temp1[$select_E];
					}else{
						$hash_v{$aa}{$ff}=$hash_v{$aa}{$ff}+$Temp1[$select_E];
					}	
				}else{
					$hash_v{$aa}{$ff}=$Temp1[$select_E];
				}
			}
		}
	}
}

my %hash_va=();
my $fi_1;
foreach my $kk (@etc_all_file){
	foreach my $kk1 (keys (%{$hash_v{$kk}})){
		if ($kk1!~/$fi/){
			if (!$hash_va{$kk1}){
				$hash_va{$kk1}="\t".$hash_v{$kk}{$kk1};	
			}else{
				$hash_va{$kk1}=$hash_va{$kk1}."\t".$hash_v{$kk}{$kk1};
			}
			
		}
	}
	$fi_1=$fi_1."\t${kk}";
}
print ans2 "$fi\t$hash_a{$fi}$fi_1\n";
foreach my $kk (keys (%hash_va)){
	my $aa1=join("\t", @{$hash_a{$kk}}[0..1]);
	print ans2 "$kk\t$aa1$hash_va{$kk}\n";
}
close (ans2);
return %hash_a;

};


sub Stringtie_Transcript{ # combine all Transcripts of samples
my $dir1 = getcwd;
$dir = "ballgown";
opendir DIR, $dir;
@etc_all_file = readdir DIR;
closedir DIR;
@etc_all_file = grep {!/\./} @etc_all_file;
my %hash_a;
my %hash_aa;
my %hash_v=();
open (ans2,">T_rawdata.txt");

my $fi;
my $fi_1;
foreach my $aa(@etc_all_file){
	open(IN, "<", "$dir1/ballgown/${aa}/t_data.ctab") or die "Could not open $name\n\n";
	my @IN=<IN>;close (IN);
	foreach my $le(@IN){
		chomp $le;
		my @Temp1=split /\t/, $le;
		if ($Temp1[0]=~/^t_id/){
			$fi=$Temp1[5];
			$hash_a{$Temp1[5]}="$Temp1[8]\t$Temp1[9]\tlocus\t$Temp1[7]";
			$hash_aa{$Temp1[5]}="$Temp1[8]\tlocus";
		}else{
			$hash_a{$Temp1[5]}=();
			push @{$hash_a{$Temp1[5]}},"$Temp1[8]"; # gene id 
			push @{$hash_a{$Temp1[5]}},"$Temp1[9]"; # symbol
			push @{$hash_a{$Temp1[5]}},"$Temp1[1]:$Temp1[3]-$Temp1[4]($Temp1[2])"; # loc
			push @{$hash_a{$Temp1[5]}},"$Temp1[7]"; # exon number 

			push @{$hash_aa{$Temp1[5]}},"$Temp1[8]"; # gene id 
			push @{$hash_aa{$Temp1[5]}},"$Temp1[1]:$Temp1[3]-$Temp1[4]($Temp1[2])"; # loc
			
			if($Temp1[11]<=0.00000001){
				$hash_v{$aa}{$Temp1[5]}=0.0001;
			}else{
				$hash_v{$aa}{$Temp1[5]}=$Temp1[11];
			}
		}
	}
}
my $hash_va=();
foreach my $kk (@etc_all_file){
	foreach my $kk1 (keys (%{$hash_v{$kk}})){
		if ($kk1 !~/$fi/){
			if(!$hash_va{$kk1}){$hash_va{$kk1}="\t".$hash_v{$kk}{$kk1};
			}else{
				$hash_va{$kk1}=$hash_va{$kk1}."\t".$hash_v{$kk}{$kk1};}
		}
	}
	$fi_1=$fi_1."\t${kk}";
}
print ans2 "$fi\t$hash_a{$fi}$fi_1\n";
foreach my $kk (keys (%hash_va)){
	my $aa1=join("\t", @{$hash_a{$kk}}[0..3]);
	print ans2 "$kk\t$aa1$hash_va{$kk}\n";
}
close (ans2);
return (%hash_aa);

}

sub Stringtie_Transcript2{ # combine all Transcripts of samples
	my $dir1 = getcwd;
	$dir = "ballgown";
	opendir DIR, $dir;
	@etc_all_file = readdir DIR;
	closedir DIR;
	@etc_all_file = grep {!/\./} @etc_all_file;
	my @hash_sa=();
	my @hash_t=();
	my %hash_a=();
	my %hash_aa=();
	my %hash_s=();
	my %hash_ss=();
	my $a1=0;
	foreach my $aa(@etc_all_file){
		open(IN, "<", "$dir1/ballgown/${aa}/${aa}.gtf") or die "Could not open $name\n\n";
		my @IN=<IN>;close (IN);
		push @hash_sa,$aa; #sample
		foreach my $le(@IN){
			chomp $le;
			my @Temp1=split /\t/, $le;
			if (defined($Temp1[2]) && $Temp1[2]=~/transcript/){
				my $TPM;
				my $tra;
				my @Temp2=split /;/, $Temp1[8];
				if($Temp2[1]=~/transcript_id\s\"(.+)\"/){
					$tra=$1;
				}

				if($Temp2[0]=~/gene_id\s\"(.+)\"/){
					$hash_s{$tra}{'gene'}=$1;
				}elsif($Temp2[1]=~/transcript_id\s\"(.+)\"/){
					$tra=$1;
				}

				if($Temp2[2]=~/gene_name\s\"(.*)\"/){
					$hash_s{$tra}{'symbol'}=$1;
				}else{
					$hash_s{$tra}{'symbol'}='-';
				}

				if($Temp2[$#Temp2]=~/TPM\s\"(.+)\"/){
					if($1==0){
						$TPM=0.0001;
					}else{
						$TPM=$1;
					}
					push @{$hash_a{$tra}},$TPM;
					if(defined($hash_aa{$hash_s{$tra}{'gene'}}[$a1])){
						if($hash_aa{$hash_s{$tra}{'gene'}}[$a1]==0.0001){
							if($1==0){
								$hash_aa{$hash_s{$tra}{'gene'}}[$a1]=$hash_aa{$hash_s{$tra}{'gene'}}[$a1];
							}else{
								$hash_aa{$hash_s{$tra}{'gene'}}[$a1]=$1;
							}
						}else{
							$hash_aa{$hash_s{$tra}{'gene'}}[$a1]=$hash_aa{$hash_s{$tra}{'gene'}}[$a1]+$1;
						}
					}else{
						push (@{$hash_aa{$hash_s{$tra}{'gene'}}},$TPM);
					}
				}

				if(!defined($hash_s{$tra}{'s'})){
					$hash_s{$tra}{'s'}=$Temp1[6];
					$hash_s{$tra}{'chr'}=$Temp1[0];
					$hash_s{$tra}{'star'}=$Temp1[3];
					$hash_s{$tra}{'end'}=$Temp1[4];
					$hash_ss{$hash_s{$tra}{'gene'}}{'star'}=$Temp1[3];
					$hash_ss{$hash_s{$tra}{'gene'}}{'end'}=$Temp1[4];
					push @hash_t,$tra;
				}else{
					if($hash_ss{$hash_s{$tra}{'gene'}}{'star'}>$Temp1[3]){
						$hash_ss{$hash_s{$tra}{'gene'}}{'star'}=$Temp1[3]
					}
					if($hash_ss{$hash_s{$tra}{'gene'}}{'end'}<$Temp1[4]){
						$hash_ss{$hash_s{$tra}{'gene'}}{'end'}=$Temp1[4]
					}
				}		
			}
		}
		$a1=$a1+1;
	}
	my $fi=join('	',@hash_sa);
	print("AAA\t$fi\n");
	open (ans2,">T_rawdata.txt");
	open (ans3,">rawdata2.txt");
	print ans2 "Transcript ID\tGene ID\tGene Name\tlocus\t$fi\n";
	print ans3 "Gene ID\tGene Name\tlocus\t$fi\n";
	#my %hash_gl;
	my %hash_tl;
	$hash_tl{'t_name'}="gene_id\tlocus";
	foreach my $kk (sort @hash_t){
		print ans2 "$kk\t$hash_s{$kk}{'gene'}\t$hash_s{$kk}{'symbol'}\t$hash_s{$kk}{'chr'}:$hash_s{$kk}{'star'}-$hash_s{$kk}{'end'}($hash_s{$kk}{'s'})\t".join('	',@{$hash_a{$kk}})."\n";
		push @{$hash_tl{$kk}},"$hash_s{$kk}{'gene'}"; # gene id 
		push @{$hash_tl{$kk}},"$hash_s{$kk}{'chr'}:$hash_s{$kk}{'star'}-$hash_s{$kk}{'end'}($hash_s{$kk}{'s'})"; # loc
		#if(!defined($hash_gl{$hash_s{$kk}{'gene'}})){
		print ans3 "$hash_s{$kk}{'gene'}\t$hash_s{$kk}{'symbol'}\t$hash_s{$kk}{'chr'}:$hash_ss{$hash_s{$kk}{'gene'}}{'star'}-$hash_ss{$hash_s{$kk}{'gene'}}{'end'}($hash_s{$kk}{'s'})\t".join('	',@{$hash_aa{$hash_s{$kk}{'gene'}}})."\n";
			#$hash_gl{$hash_s{$kk}{'gene'}}="$hash_s{$kk}{'chr'}:$hash_ss{$hash_s{$kk}{'gene'}}{'star'}-$hash_ss{$hash_s{$kk}{'gene'}}{'end'}($hash_s{$kk}{'s'})";
		#}
	}
	close ans2;
	close ans3;
	return(%hash_tl);
}




return 1;
