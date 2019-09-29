#!/usr/bin/perl 
package processPathway;
use POSIX;

#### First of all get differentail expression of genes to remove those genes that are significantly downregulated. 
sub analyze {
	my ($path,$Gpath,$Ppath,$p_val,$logFC_threshold,$PCCthreshold,$method,$iter,$GC,$Rscript,$seed,$status,$header,$show,$Up_down,$rnaseq) = @_; ##datexprs1 = normal and datexprs2 = diseased
	print "Performing DC and DE analysis...\n";
	($DE,$LOGFC) = &getFC($path,$Gpath,$Ppath,$p_val,$logFC_threshold,$PCCthreshold,$method,$GC,$Rscript,$seed,$header,$Up_down,$rnaseq);
	%DE = %$DE;
	%LOGFC = %$LOGFC;
	print "Loading Gene Regulatory Network...\n";
	($GRN,$GRN2) = &GetGRN($path,$Gpath);
	%GRN = %$GRN; %GRN2 = %$GRN2;
	@G = keys %GRN;
	print "Loading Gene Network properties...\n";
	open (GS,"$path/GS.ssv") or die ("Unable to open gene strength $!");
	while ($x=<GS>)
	{
	my($temp1, $temp2) = split(" ",$x);
	$temp1=~s/\s|\n//gi;$temp2=~s/\s|\n//gi;
	$Rs{$temp1} = $temp2;
	}
	close GS;
	open (SC,">$path/score.tsv") or die ("Cannot open score file to write $!");
	print SC "Pathways\tDensity\tRewiringScore\tGene count\n";
	#$RWGenes = scalar(keys %Rs);
	
	open (PRO,"$Ppath") or die ("Cannot open pathway file $!");
	@allpathways = <PRO>;
	for ($l = 0 ; $l <=$#allpathways; $l++)
	{
	$x = $allpathways[$l];
	($name,$Gnes) = split(/\|/, $x);
	$name=~s/\s|\n//gi;#$Gnes=~s/\n//gi;
	#@genes = split("\t",$Gnes);
	$na = "P"."$l";
	$annotation{$name} = $na;
	$score = 0; $R1 = 0; $R2 = 0; $edgeC = 0; $actual = 0; $yes = 0; %det = ();$Totexprs = 0;$significant = (); $scoreContrib = 0; $contrib = 0;$justRW=0;#%Gcount = ();
	$fileD = "$path"."/html/$na".".sif.txt"; ## reading the network file for this pathway; we will later rewrite this file to mark regulatory edges
		if (-e $fileD)
		{
		open(SNET, "$fileD") or print ("No file $fileD found $!");
		print "|---$na: $name\n";
		@array = <SNET>;
		close SNET;
		open(WR1, ">$fileD") or print ("unable to write $fileD: $!");
		$total = scalar(@array);
		$RWoNLY = 0; $tfRW = 0; $tfRW_outside = 0; $TFCount = 0;
			for $i (1..$total)
			{
				$x = $array[$i];
				($g1,$g2,$r,$p) = split(" ",$x);
				$g1=~s/\s|\n//gi;$g2=~s/\s|\n//gi;$p=~s/\s|\n//gi;$r=~s/\s|\n//gi;
				$gene2pathways{$g1}{$name} = 1;
				$gene2pathways{$g2}{$name} = 1;
				$ND{$g1}{$g2} = $p;
				$ND{$g2}{$g1} = $p;
				$PD{$name}{$g1}++;
				$PD{$name}{$g2}++;
				$det{$g1}{$g2} = 0; $det{$g2}{$g1} = 0;
				$Rsco += $r;
				if ($p >= 1.30103) ##RW => REWIRED (since p value is log transformed with -log10 so p value of 0.05 becomes 1.30103
				{
					#$significant++;
					$Rsco_SUB += $r; 
					$RPP{$name}{$g1} += $r; $RPP{$name}{$g2} += $r; ## rewiring score of each gene per pathway (rewiring per pathway); higer the score higher the number rewired edges 
					if ((exists $DE{$g1}) && ($GRN{$g1}{$g2} eq "TF-TF"))
					{
						## $g1 - $g2 is regulatory interaction
						print WR1 "$g1 $g2 $p RW REG\n"; ##REGULATORY LINK ## BOTH ARE TFs
						$contrib += ($LOGFC{$g1});
						#$contrib += $Rs{$g1};
						$scoreContrib++;
						$Description{$name}{$g1} = "TF";
						$Description{$name}{$g2} = "TF";
						$tfRW = 1;
						$TFCount++; #$Gcount{$g1}++; $Gcount{$g2}++;
					}
					if ((exists $DE{$g2}) && ($GRN{$g2}{$g1} eq "TF-TF"))
					{
						## $g1 - $g2 is regulatory interaction
						print WR1 "$g1 $g2 $p RW REG\n"; ##REGULATORY LINK ## BOTH ARE TFs
						$contrib += ($LOGFC{$g2});
						#$contrib += $Rs{$g2};
						$scoreContrib++;
						$Description{$name}{$g1} = "TF";
						$Description{$name}{$g2} = "TF";
						$tfRW = 1; $TFCount++; #$Gcount{$g1}++; $Gcount{$g2}++;
					}
					elsif ((exists $DE{$g2}) && ($GRN{$g2}{$g1} eq "TF-TG"))
					{
						## $g2 - $g1 is regulatory interaction
						print WR1 "$g1 $g2 $p RW REG\n"; ##REGULATORY LINK ## BOTH ARE TFs
						#$contrib +=$Rs{$g2};
						$contrib += ($LOGFC{$g2});
						$scoreContrib++;
						if ($Description{$name}{$g2} ne "TF"){$Description{$name}{$g2} = "TG"};
						$Description{$name}{$g1} = "TF";
						$tfRW = 1; $TFCount++; #$Gcount{$g1}++; $Gcount{$g2}++;
					}
					elsif ((exists $DE{$g1}) && ($GRN{$g1}{$g2} eq "TF-TG"))
					{
						## $g2 - $g1 is regulatory interaction
						print WR1 "$g1 $g2 $p RW REG\n"; ##REGULATORY LINK ## BOTH ARE TFs
						$contrib += ($LOGFC{$g1});
						#$contrib += $Rs{$g1};
						$scoreContrib++;
						if ($Description{$name}{$g1} ne "TF"){$Description{$name}{$g1} = "TG"};
						$Description{$name}{$g2} = "TF";
						$tfRW = 1; $TFCount++; #$Gcount{$g1}++; $Gcount{$g2}++;
						
					}
					else
					{
						if (($g1=~/[0-9]|[a-z]/gi) && ($g2=~/[0-9]|[a-z]/gi))
						{
							#$Gcount{$g1}++; $Gcount{$g2}++;
							
							print WR1 "$g1 $g2 $p RW NA\n";
							## NON-REGULATORY BUT REWIRED INTERACTION
							$RWoNLY = 1;
							if (($Description{$name}{$g2} ne "TF") && ($Description{$name}{$g2} ne "TG"))
							{
								$Description{$name}{$g2} = "OTHER";
								$significant++;
								$justRW++;
							}
							elsif (($Description{$name}{$g1} ne "TF") && ($Description{$name}{$g1} ne "TG"))
							{
								$Description{$name}{$g1} = "OTHER";
								$significant++;
								$justRW++;
							}
							elsif (($Description{$name}{$g1} ne "TF") && ($Description{$name}{$g1} ne "TF"))
							{
								$Description{$name}{$g1} = "OTHER";
								$Description{$name}{$g2} = "OTHER";
								$significant++;
								$justRW++;
							}
						}
						
					}
				}
				else
				{
				## NON-REWIRED INTERACTION
					if (($g1=~/[0-9]|[a-z]/gi) && ($g2=~/[0-9]|[a-z]/gi))
					{
					print WR1 "$g1 $g2 $p UNK NA\n";
						if (($Description{$name}{$g2} ne "TF") && ($Description{$name}{$g2} ne "TG"))
						{
							$Description{$name}{$g2} = "OTHER";
						}
						if (($Description{$name}{$g1} ne "TF") && ($Description{$name}{$g1} ne "TG"))
						{
							$Description{$name}{$g1} = "OTHER";
						}
					}
				}
			}
			
			foreach $G (keys %{$Description{$name}}) ## reading all genes one by one for the last processed pathway and attempting to find rewired TF links that occured from outside the pathway 
			{
				@TFs = keys %{$GRN{$G}};
				#print "$#TFs => ";
				foreach $myT (@TFs)
				{

					if ((exists $det{$G}{$myT}) || (exists $det{$myT}{$G})) ## ALREADY IN THE EDGE LIST written BEFORE for this pathway 
					{}
					else
					{
						if (($GRN{$G}{$myT} eq "TF-TF") || ($GRN{$myT}{$G} eq "TF-TF")) ## a given pathway TF is mediated by other pathway TF 
						{
						$Description{$name}{$myT} = "TFO";
						if ($Description{$name}{$G} ne "TF"){$Description{$name}{$G} = "TG"};
						## REGULATORY LINK THAT ARE NOT THE PART OF THIS PATHWAY}
						
						if (exists $GRN2{$myT}{$G})
						{$ND{$myT}{$G} = $GRN2{$myT}{$G};
						print WR1 "$G $myT $GRN2{$myT}{$G} RW REGO\n";}
						else
						{$ND{$G}{$myT} = $GRN2{$G}{$myT};
						print WR1 "$G $myT $GRN2{$G}{$myT} RW REGO\n";}
						#$significant++;
						$scoreContrib++; ## counting the number of rewired interactions mediated by TFs
						$contrib += ($LOGFC{$G});
						#$contrib += $Rs{$G};
						$tfRW_outside = 1;  $TFCount++;
						#$Gcount{$G}++;
						}
						elsif ($GRN{$G}{$myT} eq "TF-TG") ## a given pathway TG is mediated by other pathway TF 
						{
							if (exists $DE{$G})
							{ 
							$Description{$name}{$myT} = "TFO";
							if ($Description{$name}{$G} ne "TF"){$Description{$name}{$G} = "TG"};
							$contrib += ($LOGFC{$G}); ## total pathway upregulation due to TF mediated rewiring 
							#$contrib += $Rs{$G};
							print WR1 "$G $myT $GRN2{$G}{$myT} RW REGO\n";## REGULATORY LINK THAT ARE NOT THE PART OF THIS PATHWAY}
							$scoreContrib++; ## counting the number of rewired interactions mediated by TFs
							$ND{$G}{$myT} = $GRN2{$G}{$myT};
							$tfRW_outside = 1; 
							#$significant++; 
							$TFCount++;
							#$Gcount{$G}++;
							}
						}
					}
				}
				@TFs = ();
			}
			@allG = keys %det;
			foreach $g (@allG)
			{  
				$Totexprs += ($LOGFC{$g});
				$RWscore += $Rs{$g};
			}
			### checking type of rewiring that happened 
			if (($tfRW_outside == 1) && ($tfRW == 1))
			{
				$target{$name} = "TF-TFO"; ## THIS PATHWAY IS REGULATED BY BOTH INTER- AND INTRA-PATHWAY TF mediated rewiring 
				
			}
			elsif (($tfRW_outside == 0) && ($tfRW == 1))
			{
				$target{$name} = "TF"; ## THIS PATHWAY IS REGULATED BY only INTRA-PATHWAY TF mediated rewiring
				
			}
			elsif (($tfRW_outside == 1) && ($tfRW == 0))
			{
				$target{$name} = "TFO"; ## THIS PATHWAY IS REGULATED BY only INTER-PATHWAY TF mediated rewiring; TF does belongs among the pathway genes 
			}
			elsif (($justRW > 5) && ($tfRW_outside == 0) && ($tfRW == 0))
			{
				$target{$name} = "other";
			}
			if ($justRW > 5)
			{
				#############################################
				##USELESS: $TOTAL ; $div2 ; $TFpotential
				#$PC = scalar(keys %Gcount);
				$RewiredTF = $scoreContrib / $significant; 
				
				$GRNfraction = $contrib/$Totexprs; ## MAKES IT INDEPENDENT OF CORRELATION CHANGE 
				$div = $Rsco_SUB/$Rsco; ## TOTAL REWIRED EDGE STRNGTH / TOTAL EDGE STRENGTH
				$Fscore  = (($GRNfraction * $scoreContrib) + $div);
				
				### To test the role of rewiring in detecting the pathway abberation please set##
				#$Fscore  = $GRNfraction ; ## and see if pathways are still being predicted or not 
				
				print SC "$name\t$Fscore\t$GRNfraction\t$div\t$RewiredTF\t$scoreContrib\n";
				if (($GRNfraction > -1) && ($Fscore > -1))
				{
					$Dysregulation{$name} = "$Fscore,$Totexprs,$div,$GRNfraction";
					push(@dysscore,$Fscore);
				}
			}
			@allG = (); %det = (); $RWscore = 0; $Rsco_SUB = 0; $Rsco = 0;
		}
	}

	@C = keys %Dysregulation;
	if (scalar(@C) >= 1)
	{
		@sorted = sort {$a <=> $b} @dysscore;
		$minDys = $sorted[0];$maxDys = $sorted[$#sorted]; $diffDys = $maxDys - $minDys;
		#now calculate the disease enrichment score for each pathway
		print "Calculating Disease potential of each pathway...\n"; 
		(%affinity) = &oncology($path,\%PD,\%ND,$method);
		#Now for each disease enriched pathway calculate the total rewired network centrality ...
		print "Performing Final analysis...\n";
		($FH,$CS,$pcount,$DS,$DP,$out_inside,$outside,$inside,$RWpathfound,$DRfound,$nonDRfound) = &display($path,\%affinity,\%gene2pathways,\%Dysregulation,\%annotation,\%target,\%PD,$minDys,$maxDys,$diffDys,$status);
		print "Creating charts...\n";
		&pie_chart($path,$out_inside,$outside,$inside,$pcount,$nonDRfound);
		&bar_chart($path,$RWpathfound,$DRfound);
		&scatterPlot($path,$DS,$DP);
		print "Printing per pathway report... \n";
		&details($path,\%annotation,$FH,\%exprs,\%LOGFC,$CS,\%RPP,\%Rs,\%Description,\%PD,$status,$show);
		print "Cleaning files... \n";
		#&remove($path,$FH,\%annotation);
			if ($iter > 9)
			{
				$script = "src/iterate.R";
				`$Rscript $script $path $iter $header $Up_down $Ppath`; ## here iter is the number of iterations
				&iterate::Iteration($path,$logFC_threshold,$iter,$p_val,$Gpath);
			}
		print "\n########### \n Done !!! \n########### \n";
	}
	else
	{
		print "Sorry!!! No dysregulated pathway found for given parameters\n";
		exit;
	}
}

sub getFC {
my($setwd,$Gpath,$Ppath,$p_val,$FC,$PCC,$method,$GC,$Rscript,$seed,$header,$Up_down,$rnaseq) = @_;
	$script  = "src/de.R"; 
	$sep=" \"\\t\"";
	#`/usr/local/bin/Rscript $script $setwd $Gpath $Ppath`; ## Executing new R script  	 ###<-------------------------------------------
	print "$Rscript $script $setwd $Gpath $Ppath $method $GC $PCC $seed $header $Up_down $rnaseq\n";
	`$Rscript $script $setwd $Gpath $Ppath $method $GC $PCC $seed $header $Up_down $rnaseq`;
	open(FC, "$setwd"."/DE.tsv") or die ("Error: Cannot open DE results, $!"); ## opening output
	print "Creating DE hash ...\n";
	%DE = ();
	while (my $x = <FC>)
	{
		my @data = split("\t",$x);
		$data[0]=~s/\s|\n//gi; ## geneID
		$data[2]=~s/\s|\n//gi; ## Fold change 
		$data[5]=~s/\s|\n//gi; ## p.value 
		$TOTAL{$data[0]} = $data[2];
		if(($data[2] >= $FC) && ($data[5] <= $p_val))
		{ 
		## remove those genes which are downregulated or whose p-value is greater than defined threshold 
		$DE{$data[0]} = $data[2];
		} ## if adjusted p-value is less than or equal to 0.05 get the fold change
		else
		{$DE{$data[0]} = 0;}
	}
	close FC;
return(\%DE,\%TOTAL);
}

sub GetGRN{
my($path,$Gpath) = @_;
open(FH,"$Gpath") or die ("Unable to open GRN file: $!");;
@arr = <FH>;
%Hsh=();
	foreach $x (@arr)
	{
		($temp1,$temp2,$tfg) = split(" ",$x);
		$temp1=~s/\s|\n//gi;$temp2=~s/\s|\n//gi;$tfg=~s/\s|\n//gi;
		
		if ($tfg eq "TF-TF")
		{
		$Hsh{$temp1}{$temp2} = $tfg;
		}
		else
		{
		$Hsh{$temp2}{$temp1} = $tfg; ## key1 = TG and Key2 = TF and tfg="TF-TG" or "TF-TF"
		}
	}
@arr = ();
close FH;
open(FH,"$path/GRN.ssv");
@arr = <FH>;
%TFG= ();%TFG2= ();
	foreach $x (@arr)
	{
		($temp1,$temp2,$score) = split(" ",$x);
		$temp1=~s/\s|\n//gi;$temp2=~s/\s|\n//gi;$score=~s/\s|\n//gi;
		if (exists $Hsh{$temp2}{$temp1})
		{
			$tfg=$Hsh{$temp2}{$temp1};
			if ($tfg eq "TF-TF")
			{
			$TFG{$temp2}{$temp1} = $tfg; ## tg-tf
			$TFG2{$temp2}{$temp1} = $score; ## tg-tf
			}
			else
			{
			$TFG{$temp2}{$temp1} = $tfg; ## key1 = TG and Key2 = TF and tfg="TF-TG" or "TF-TF"
			$TFG2{$temp2}{$temp1} = $score;
			}
		}
		elsif (exists $Hsh{$temp1}{$temp2}) ## tg-tf
		{
			$tfg=$Hsh{$temp1}{$temp2};
			if ($tfg eq "TF-TF")
			{
			$TFG{$temp1}{$temp2} = $tfg; ## tg-tf
			$TFG2{$temp1}{$temp2} = $score;
			}
			else
			{
			$TFG{$temp1}{$temp2} = $tfg; ## key1 = TG and Key2 = TF and tfg="TF-TG" or "TF-TF"
			$TFG2{$temp1}{$temp2} = $score; 
			}
		}
	}
@arr = (); %Hsh = ();
close FH;
return (\%TFG,\%TFG2)
}

sub oncology {
	my ($path,$Pathw,$nd,$method) = @_;
	my %Pathways = %$Pathw; ## hash of pathway and its gene components 
	#%PV = %$P; ## hash of pathway and its gene components 
	%ND = %$nd;
	my @genes  = ();
	## Now we will predict the disease gene enrichment assosiated with each pathway using random walk and network propogation/gene neighbourhood ###
	## created a disease network before ... Now I am going to use it for predicting the disease gene with rewired connection
	### random walk by restart ####	
	open (RWR, "$path"."/RWR.txt") or die ("cannot open affinity matrix : $!");
	@rwr = <RWR>;
	foreach $r (@rwr)
	{
		@temp = split(" ",$r);
		$ID = $temp[0];
		$affinity = $temp[1];
		$ID=~s/\s|\n//gi;
		$affinity=~s/\s|\n//gi;
		if ($affinity > 0)
		{
		$gene{$ID} = $affinity;
		}
		else
		{
		$gene{$ID} = 0;
		}
	}
	%scoreT = (); $score = 0;$rewired=0;
	
	foreach $p (keys %Pathways)
	{
		@genes  = keys %{$Pathways{$p}};
		$score = 0;
		%scoreT = (); $edges = 0;
		foreach $g (@genes)
		{
			foreach $o (@genes)
			{
				if (($ND{$g}{$o} >= 1.30103) && (($gene{$g} > 0) || ($gene{$o} > 0))) ## if interaction is rewired and one or both of the gene is disease gene; if $g == 1 then disease gene is within pathway .. if $o == 1 then disease gene may be outside
				{
					if (!exists $scoreT{$g}{$o}) ## make sure edge is not counted twice 
					{
						$scoreT{$g}{$o}++;
						$scoreT{$o}{$g}++;
						$score += ($gene{$g} + $gene{$o}); ## weighing each edge if it is rewired and at-least one of the gene is disease gene 
						#$score += $ND{$g}{$o}; ## an option ; instaed of summing up the prediction score .. just add the strength of rewiring for disease edge
					}
				}
				if($ND{$g}{$o} >= 1.30103)
				{
					$edges++;
				}
			}
		}
		if ($score > 0) ## if atleast one pathway gene is rewired and made connections with disease gene..
		{
		#$new{$p} = $score / scalar(@genes); ## total number of oncogenic rewired connections in pathway / total number of rewired connection; i.e. out of how many rewired edges oncogenic connections were generated in this pathway  
		$new{$p} = $score / $edges;
		}
	 ## score means total score for all the rewired genes that gained oncogenic connections (>= 1)
	}
	@valu = values %new;
	@ascending = sort { $a <=> $b } @valu;
	$min = $ascending[0]; $max = $ascending[$#ascending];
	@keys = keys %new;
	
	foreach $k (@keys)
	{
		$Xi = $new{$k};
		$z = ($Xi - $min)/($max - $min);
		$Onco{$k} =  $z;
	}
	
return(%Onco);
}

sub display {
	my ($path,$aff,$pd,$D,$anno,$tar,$pathDet,$minDys,$maxDys,$diffDys,$status) = @_;
	%affinity = %$aff;
	%GP = %$pd; ## genes 2 pathways 
	%Dys = %$D;
	%annotation = %$anno;
	%target = %$tar;
	%PD = %$pathDet;
	$pcount = 0; $out_inside =0; $outside =0; $inside=0;$DRfound=0;$RWpathfound=0;
	### open RWR.txt and centrality.txt ###
	## print => pathwayname \t Dysregulation score \t Disease enrichment \t pathway significance (centrality measures) 
	## other things that may be of ineterset => cumulative fold change \t fraction of rewired edges \t number of connected genes in disease network [this information can be added in downloadable file]
	open(FR,"$path/centrality.txt") or die ("Unable to open centalities $!");
	@cen = <FR>;
	open(N1,">$path/result.tsv") or die ("Unable to open results to write $!");
	print N1 "PathwayName Rank Disease_gene_enrichment Dy_Score DR_score Gene_count\n";
	#&print_start($path,$p_val,$logFC_threshold,$PCCthreshold,$datExprs1,$datExprs2,$selected); ## printing the HTML header section
	foreach $c (@cen)
	{
		($ID,$CenScore,$C1,$C2,$C3,$C4) = split(/\t/,$c);
		$ID=~s/\s|\n//gi; $CenScore=~s/\s|\n//gi; $C1=~s/\s|\n//gi; $C2=~s/\s|\n//gi;$C3=~s/\s|\n//gi;$C4=~s/\s|\n//gi;
		@pathway = keys %{$GP{$ID}}; ## a gene can exists in multiple pathways 
		foreach $x (@pathway)
		{
			$cen{$x} += $CenScore; ## overall centrality score for a pathway
		}
		
		$CentralityScores{$ID} = "$CenScore,$C1,$C2,$C3,$C4";
	}
	@ASC  = values %cen;
	@sorted  = sort{ $a <=> $b } @ASC;
	$min = $sorted[0];
	$max = $sorted[$#sorted];
	@pathways = keys %Dys;
	foreach my $pathway (@pathways)
	{
		if (exists $affinity{$pathway})
		{
			$DiseasePotential = $affinity{$pathway};
		}
		else
		{
			$DiseasePotential = 0;
		}
		if (exists $cen{$pathway})
		{
			$cScore = $cen{$pathway};
			$cenScore = ($cScore - $min)/($max - $min);
		}
		else
		{
			$cenScore = 0;
		}
		
		####################
		## calculate total number of genes in this pathway and get avergae centrality for this pathway  
		####################
		
		($Fscore,$FC,$RW,$TFscore) = split(",",$Dys{$pathway});
		$FinalHash{$pathway} = "$Fscore,$DiseasePotential,$cenScore,$FC,$TFscore";
		$n = $annotation{$pathway};
		$FFscore = ($Fscore - $minDys)/($diffDys);
		if (($FFscore > -1) && ($RW > 0)) ## retaining all pathways with atleast one rewired edge  
		{
		if ($target{$pathway} eq "TF-TFO"){$color = "#DF0101";$out_inside++;} ## green
		elsif ($target{$pathway} eq "TFO"){$color = "#0B3B0B";$outside++;} ## red 
		elsif ($target{$pathway} eq "TF"){$color = "#0000FF";$inside++;} ##  blue
		elsif ($target{$pathway} eq "other"){$color = "#000000"; $nonDRfound++;} ## black
		$GeneCount = scalar(keys %{$PD{$pathway}});
			if (($status eq "F") || (($status eq "T" && ($TFscore > 0.05))))
			{
				if ($pathway=~/dummy/gi)
				{
					next;
				}
				else
				{
				&print_table($n,$path,$pathway,$FFscore,$DiseasePotential,$TFscore,$color,$GeneCount,$RW);
				$pcount++; ## total nnumber of dysregulated pathways
				push(@DS,$TFscore);
				push(@DP,$RW);
				print N1 "$n $pathway $FFscore $DiseasePotential $Fscore $TFscore $GeneCount\n";
				if ($TFscore > 0.0) {$DRfound++;}
				$RWpathfound++;
				}
			}
		}
	}
	close WRN1; close FR;
	&print_stop($path);
return(\%FinalHash,\%CentralityScores,$pcount,\@DS,\@DP,$out_inside,$outside,$inside,$RWpathfound,$DRfound,$nonDRfound);
}

sub print_start {
	($path,$p_val,$logFC_threshold,$PCCthreshold,$datExprs1,$datExprs2,$selected) = @_;

	open(WR, ">$path"."/index.html") or print ("Cannot write html file $!");
	print WR "
	<!DOCTYPE html>
<html lang=\"en\">
    <head>
		<meta charset=\"UTF-8\" />
        <meta http-equiv=\"X-UA-Compatible\" content=\"IE=edge,chrome=1\"> 
        <title>Result</title>
        <link rel=\"shortcut icon\" href=\"favicon.ico\"> 
        <link rel=\"stylesheet\" type=\"text/css\" href=\"src/table.css\" />
		<link href=\"src/demo.css\" rel=\"stylesheet\" type=\"text/css\">
		<link href=\"src/arrow.css\" rel=\"stylesheet\" type=\"text/css\">
		<script src=\"src/tsorter.js\" type=\"text/javascript\"></script>
		<script type=\"text/javascript\">
						function init()
						{
							var sorter1 = tsorter.create('housing_table_1')
						}
					window.onload = init;
						</script>
    </head>
    <body>
    <script type=\"text/javascript\" src=\"src/wz_tooltip.js\"></script>
        <div class=\"container\" >          
			<header>
				<h1><span><i>DPA Result</i></span></h1>
				<h2><b>D</b>ysregulated <b>P</b>athway <b>A</b>nalyzer</h2>
			</header>
			<section class=\"tabs\">
	            <input id=\"tab-1\" type=\"radio\" name=\"radio-set\" class=\"tab-selector-1\" checked=\"checked\" />
		        <label for=\"tab-1\" class=\"tab-label-1\">Input</label>
		
	            <input id=\"tab-2\" type=\"radio\" name=\"radio-set\" class=\"tab-selector-2\" />
		        <label for=\"tab-2\" class=\"tab-label-2\">Results</label>
		
	            <input id=\"tab-3\" type=\"radio\" name=\"radio-set\" class=\"tab-selector-3\" />
		        <label for=\"tab-3\" class=\"tab-label-3\">Details</label>
			
	            <input id=\"tab-4\" type=\"radio\" name=\"radio-set\" class=\"tab-selector-4\" />
		        <label for=\"tab-4\" class=\"tab-label-4\">Contact</label>
<!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------->     
			    <div class=\"clear-shadow\"></div>
		        <div class=\"content\">
			        <div class=\"content-1\">
						<h3>Dataset</h3>
                        <p>Control expression matrix: $datExprs1 <br> Genes: $selected </p>
                        <hr width= 60% align=\"left\" COLOR=\"grey\" SIZE=\"1\">
						<p>Disease expression matrix: $datExprs2 <br> Genes: $selected </p>
						<h3>Options selected:</h3>
						<p>P value threshold: $p_val </p>
						<p>Correlation thershold: $PCCthreshold </p> 
						<p>Fold change threshold: $logFC_threshold </p>
						<p>Correlation method: $method </p>
				    </div>
				    <!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------->
			        <div class=\"content-2\">
						<iframe src=\"html/pie.html\" id=\"chartContainer\" style=\" width:60%; height:300px; border:none; margin:20;overflow:hidden;\" scrolling=\"no\">  </iframe>
						<iframe src=\"html/scatter.html\" id=\"chartContainer\" style=\" width:60%; height:300px; border:none; margin:20;overflow:hidden;\" scrolling=\"no\" align=\"center\" >  </iframe>
						<iframe src=\"html/bar.html\" id=\"chartContainer\" style=\" width:60%; height:300px; border:none; margin:20;overflow:hidden;\" scrolling=\"no\">  </iframe>
						<!--
						<p> Total Number of rewired genes and edges: $RW & $EDGECOUNT</p>
						<p> Density of rewired subnetwork: $Den</p>
						<p> Total Number of TFs involved in dysregulation: $TFcount</p>
						<p> Total Number of TGs involved in dysregulation: $TGcount </p>
						-->
						</div>
						<div class=\"content-3\"> 
						<table id=\"housing_table_1\" class=\"sortable\">
								  <caption> <a href=#> Result table. Click header for column-wise sorting </a></caption>
								  <thead>
									<tr>
									  <th data-tsorter=\"link\"> </th>
									  <th data-tsorter=\"link\" onmouseover=\"Tip('Pathway Name')\" onmouseout=\"UnTip()\">Pathway</th>
									  <th data-tsorter=\"numeric\" onmouseover=\"Tip('Cumulative alteration and differential regulation score')\" onmouseout=\"UnTip()\">Rank</th>
									  <th data-tsorter=\"numeric\" onmouseover=\"Tip('Alteration (gene-gene rewiring) associated with pathway gene-set')\" onmouseout=\"UnTip()\">Dy Score</th>
									  <th data-tsorter=\"numeric\" onmouseover=\"Tip('Differential regulation of pathway gene-set')\" onmouseout=\"UnTip()\">DR Score</th>
									  <th data-tsorter=\"numeric\" onmouseover=\"Tip('Context-specific disease gene enrichment')\" onmouseout=\"UnTip()\">DP Score</th>
									</tr>
								  </thead>
								  <tbody>
	";
	close WR;
}

sub print_table {
($n,$path,$P,$rewiring,$OP,$TS,$color,$GeneCount,$RW) = @_;
	$F =  "html/"."$n".".html";
	
	
	$rewiring1 = sprintf "%.3f", $rewiring;
	$RW1 = sprintf "%.3f", $RW;
	$TS1 = sprintf "%.3f", $TS;
	$OP1 = sprintf "%.3f", $OP;
	$OP = 1 - $OP1;
	open(WR, ">>$path"."/index.html") or print ("Cannot write html file $!");
	print WR "
			<tr class=\"\" onmouseover=\"Tip('Gene Count: $GeneCount')\" onmouseout=\"UnTip()\">
			<td><p style=\"border:none;color:$color;font-size:170%;\"><span>&#9679;</span></p></td>
			 <td class=\"lft\">
			 <a href=\"$F\">$P</a>
			 </td>
			 <td>$rewiring1</td>
			 <td>$RW1</td>
			 <td>$TS1</td>
			 <td>$OP1</a></td>
			</tr>
	";
	close WR;
}

sub print_stop {
	($path) = @_;
	open(WR, ">>$path"."/index.html") or print ("Cannot write html file $!");
	print WR "
	</tbody></table></div>
				    <div class=\"content-4\">
						<h2>Contact</h2>
                        <h3>Queries and complaints</h3>
						<p>Dr. Dinesh Gupta: dinesh\@icgeb.res.in </p>
						<p>Abhinav Kaushik: abhinav\@icgeb.res.in </p>
				    </div>
		        </div>
			</section>
        </div>
    </body>
</html>
	";
}

####### The main table has been printed but the html page for each dysregulated pathway is not available ######

sub details {
	($path,$annotation,$FinalHash,$exprs,$DE,$CS,$RPP,$Rs,$Discription,$PD,$status,$show) = @_;
	%$annotation = %$annotation;
	%PD = %$PD;
	%FinalHash = %$FinalHash;
	my %DE = %$DE;
	%CS = %$CS;
	%RPP = %$RPP;
	%Rs = %$Rs;
	%exprs = %$exprs;
	%Discription = %$Discription;
	my $k = "";
	## get pathway one by one 
	## convert pathway name to ID 
	## now for each ID pick its edge file 
	## run R script to calculate the property of each pathway 
	## then print its network properties 
	## you can also generate graphs for each pathway with highlighed dysregulated edges or colored as per foldchange or coneection with disease edges 
	@pathways = keys %FinalHash; 
	$description = "";
	my $PathCount = ceil(scalar(@pathways)/4); ## for printing progress 
	my $partition = 0; ## for printing progress 
	print "0%"; ## for printing progress
	open(E2G, "database/E2GS") or print ("Unable to open database/E2GS file $! \n");
	while($e = <E2G>)
	{
		($t1,$t2) = split(" " , $e);
		$t1=~s/\s|\n//gi; $t2=~s/\s|\n//gi;
		$Convert{$t1} = $t2;
	}
	close E2G;
	foreach my $pathway (@pathways)
	{
		#if ($partition == $PathCount) ## for printing progress 
		#{
		#	$partition = 0;
		#	$progress += 25; 
		#	print "..$progress%";
		#}
		my ($RW,$OP,$centrality,$Texprs,$TFscore) = split (",",$FinalHash{$pathway});
		if (($status eq "F") || (($status eq "T" && ($TFscore > 0.05))))
		{
		$n = $annotation{$pathway};
		%prop = ();
		%prop = &getNetworkProperties($path,$pathway,$n,\%DE,\%Discription); ## this hash may have genes involved as TF that are not the part of TF
		&convert_sif_to_JS($path,$pathway,$n,\%Convert);
		@genes = keys %{$PD{$pathway}}; ## genes involved in this pathway only 
		$RWnodes = scalar(keys %{$RPP{$pathway}}); ## total rewired genes in pathway 
		$AllNodes = scalar(@genes); ## total genes in pathway 
	
		#print "$path,$n,$pathway,$RW,$OP,$#genes,$centrality,$exprs{$pathway} \n";
		&print_pathway_html($path,$n,$pathway,$RW,$OP,$AllNodes,$RWnodes,$centrality,$Texprs,$show);
		#$handler = "TABLE";
		#open(TABLE, ">$path/html/$n.csv");
		#print TABLE "GeneName,Rewiring Score,Degree,Clustering Coefficient,Eccentricity,Betweeness,Closeness,Page Rank,Fold Change"; ##the five properties Clustering Coefficient,Eccentricity,Betweeness,Closeness,pagerank have both local and global topological scores
		for $ij (0..$#genes)
		{
			$Local_Rewiring = 0; $Global_Rewiring = 0;
			$k = $genes[$ij];
				$CC = 0;$PR = 0;$betweeness = 0;$degree = 0;$Closeness=0;
				($CC,$PR,$betweeness,$degree,$Closeness) = split(/\t/,$prop{$k});
				$CC=~s/\s|\n//gi; $PR=~s/\s|\n//gi; $betweeness=~s/\s|\n//gi; $degree=~s/\s|\n//gi; $Closeness=~s/\s|\n//gi; 
				#print "Prop: $Closeness \n";
				($C,$C1,$C2,$C3,$C4) = split(",",$CS{$k});
				#print "Centre: $C == $C1 == $C2 == $C3 == $C4 \n";
				$FC = $DE{$k};
				#if ($degree eq ""){$degree = 0;}
				#if ($CC eq ""){$CC = 0;}
				#if ($betweeness eq ""){$betweeness = 0;}
				#if ($PR eq ""){$PR = 0;}
				#if ($FC eq ""){$FC = 0;}
				if(exists $RPP{$pathway}{$k})
				{
				$Local_Rewiring = $RPP{$pathway}{$k}; ## rewiring score assosiated with gene wrt this pathway sub-network only 
				}
				else
				{
				$Local_Rewiring = 0;
				}
				if(exists $Rs{$k})
				{
				$Global_Rewiring = $Rs{$k}; ## rewiring score assosiated with gene wrt whole network
				}
				else
				{
				$Global_Rewiring = 0;
				}
				if (($Local_Rewiring != 0) || ($Global_Rewiring != 0))
				{
					if (exists $Convert{$k}){$symbol = $Convert{$k};} 
					else {$symbol =$k;}
					&print_remaining($path,$n,$pathway,$k,$symbol,$degree,$CC,$betweeness,$PR,$Closeness,$FC,$description,$C1,$C2,$C3,$C4,$C,$Local_Rewiring,$Global_Rewiring,$show);
				}
		}
		close $handler;
		&print_end($path,$n);
		#$partition++;
		}
	}
	print "..100%\n";
}

sub getNetworkProperties {
	
my ($path,$name,$n,$DE,$Discription) = @_;
my $out = "$path/html/$n"."_prop.csv";
#my $network = "$path/html/$n".".sif";
my %hash = ();
%special = %$Discription;
my %DE = %$DE;
my @gen = keys %{$special{$name}};
$foldchange = "$path/html/temp.fc";
open(WRTTEMP,">$foldchange") or print ("Unable to print fold chnage value to network graph: $!");

foreach $g (@gen)
{
print WRTTEMP "$g $DE{$g} $special{$name}{$g}\n"; ## GENE TC TF/TG/TFO/OTHER

}
close WRTTEMP;
open PROP ,$out or die "$out $!";
@prop = <PROP>;
foreach $p (@prop)
{
		@arr = split(" ",$p);
		$gene = $arr[0];
		$gene=~s/\s//gi;
		$arr[1]=~s/\s//gi;$arr[2]=~s/\s//gi;$arr[3]=~s/\s//gi;$arr[4]=~s/\s//gi;$arr[5]=~s/\s//gi;$arr[6]=~s/\s//gi;
		$ClusteringCoefficient = $arr[1];
		$Degree = $arr[2];
		$PageRank = $arr[3];
		$Betweeness = $arr[4];
		$Closeness = $arr[5];
		if ($Degree > 0)
		{
		$hash{$gene} = "$ClusteringCoefficient "."\t"."$PageRank"."\t"."$Betweeness"."\t"."$Degree"."\t"."$Closeness"; ## clutering coefficient; page rank ; betweenenss ; degree of each gene
		} 
}
close PROP;
return (%hash);
}

sub print_pathway_html {
	my ($path,$n,$P,$R1,$OP1,$nodes,$RWnodes,$C1,$FC1,$show) = @_;
	$F =  "$path"."/html/"."$n".".html"; ## this is the filename containg details of each pathway
	$csv =  "$n"."_prop.csv";
	$network = "$n"."_js.html";
	$sif = "$n".".sif.txt";

	$rewiring = sprintf "%.3f", $R1;
	$FC = sprintf "%.3f", $FC1;
	$OP = sprintf "%.3f", $OP1;
	$centrality = sprintf "%.3f", $C1;
	$string = "";
	open(WR, ">$F") or print ("Cannot write html file for pathway: $!");
	$string = "<html>
<title> $P </title>
<body onload=\"document.getElementById('if').contentWindow.location.reload();\">
<script type=\"text/javascript\" src=\"../src/wz_tooltip.js\"></script>
	<head>
			<script src=\"../src/sorttable.js\"></script>
			<link rel=\"stylesheet\" type=\"text/css\" href=\"../src/table2.css\" />
			<!-- <link href=\"../src/dem2.css\" rel=\"stylesheet\" type=\"text/css\"> -->
	</head>
			<br>
		<h1> $P </h1>
<table border= 0>
<tr> 
	<td>
		<table border= 0>
			<tr>
				<th>Differential Network Properties</th>
				<th></th>
			</tr>
			<tr>
				<td>
				Rewiring Score: $rewiring
				</td>
				<td>
				<b> Mean Expression: $FC </b>
				</td>
			</tr>
			<tr>
				<td>
				Gene count: $nodes
				</td>
				<td>
				<b>Rewired gene count: $RWnodes</b>
				</td>
			</tr>
			<tr>
				<td>
				Disease Potential: $OP
				</td>
				<td>
				<b>Pathway centrality: $centrality </b>
				</td>
			</tr>
		</table>
	</td>
	<td> <a href=\"#\" onmouseover=\"Tip('Click here to know about different color codes used in ')\" onmouseout=\"UnTip()\" onclick=\"window.open('../src/colorcodes.htm', 'newwindow', 'width=600, height=300'); return false;\">Color codes </a><iframe src=$network  width=\"800\" height=\"520\" id='if' >
		<p>Your browser does not support iframes.</p> 	</iframe>
		<a href=$network target= \"_blank\"><img src=\"../src/newtab.png\" width=\"22\" height=\"22\" border=\"0\"></a>
		<a href=$sif target= \"_blank\"><img src=\"../src/download.png\" width=\"22\" height=\"22\" border=\"0\"></a>
	</td>
</tr>
</table>
<table id=\"test1\" class=\"sortable\" border=\"0\" cellpadding=\"0\" cellspacing=\"0\">
								  <caption><a href=\"$csv\">Download Table</a></caption>
								  <thead>
									<tr>";
if ($show eq "T")
{
$string .= "
									  <th style=\"-moz-user-select: none;\" class=\"sortable-keep fd-column-0\"><a>Gene Symbol</a></th>
									  <th style=\"-moz-user-select: none;\" class=\"sortable-keep fd-column-0\"><a>Rewiring Score (local)</a></th>
									  <th style=\"-moz-user-select: none;\" class=\"sortable-keep fd-column-0\"><a>Rewiring Score (Global)</a></th>
									  <th style=\"-moz-user-select: none;\" class=\"fd-column-2 sortable-text reverseSort\"><a>Degree</a></th>
									  <th style=\"-moz-user-select: none;\" class=\"sortable-keep fd-column-0\"><a>Clustering coefficient</a></th>
									  <th style=\"-moz-user-select: none;\" class=\"sortable-date-dmy fd-column-3\"><a>Betweeness</a></th>
									  <th style=\"-moz-user-select: none;\" class=\"sortable-currency fd-column-4\"><a>Page Rank</a></th>
									  <th style=\"-moz-user-select: none;\" class=\"sortable-currency fd-column-4\"><a>Closeness</a></th>
									  <th style=\"-moz-user-select: none;\" class=\"sortable-currency fd-column-4\"><a>Overall Global centrality</a></th>
									  <th style=\"-moz-user-select: none;\" class=\"sortable-numeric fd-column-5\"><a>Fold Change</a></th>
									</tr>
								  </thead>
								  <tbody>";
}
else
{
	$string .= "
									  <th style=\"-moz-user-select: none;\" class=\"sortable-keep fd-column-0\"><a>Gene Symbol</a></th>
									  <th style=\"-moz-user-select: none;\" class=\"sortable-keep fd-column-0\"><a>Rewiring Score (local)</a></th>
									  <th style=\"-moz-user-select: none;\" class=\"sortable-keep fd-column-0\"><a>Rewiring Score (Global)</a></th>
									  <th style=\"-moz-user-select: none;\" class=\"fd-column-2 sortable-text reverseSort\"><a>Degree</a></th>
									  <th style=\"-moz-user-select: none;\" class=\"sortable-keep fd-column-0\"><a>Clustering coefficient</a></th>
									  <th style=\"-moz-user-select: none;\" class=\"sortable-date-dmy fd-column-3\"><a>Betweeness</a></th>
									  <th style=\"-moz-user-select: none;\" class=\"sortable-currency fd-column-4\"><a>Page Rank</a></th>
									  <th style=\"-moz-user-select: none;\" class=\"sortable-currency fd-column-4\"><a>Closeness</a></th>
									  <th style=\"-moz-user-select: none;\" class=\"sortable-currency fd-column-4\"><a>Overall Global centrality</a></th>
								  </thead>
								  <tbody>";
}

								  
#print "$string";
print WR "$string";
close WR;
}

sub print_remaining {
	my($path,$n,$P,$gene,$symbol,$degree,$CC,$betweeness,$PR,$Closeness,$FC,$description,$C1,$C2,$C3,$C4,$C,$Local_Rewiring,$Global_Rewiring,$show) = @_;
	$LR = sprintf "%.3f", $Local_Rewiring;
	$GR = sprintf "%.3f", $Global_Rewiring;

	$C11 = sprintf "%.3f", $C1;
	$C21 = sprintf "%.3f", $C2;
	$C31 = sprintf "%.3f", $C3;
	$C41 = sprintf "%.3f", $C4;
	$C_1 = sprintf "%.3f", $C;
	$FC_1 = sprintf "%.3f", $FC;
	my $F =  "$path"."/html/"."$n".".html"; ## this is the filename containg details of each pathway
	open(WR, ">>$F") or print ("Cannot write html file for pathway in the middle: $!");
	my $string = "";
	if ($show eq "T")
	{
	$string =  "							<tr class=\"\">
												 <td class=\"lft\"><a href = \"http://www.ncbi.nlm.nih.gov/gene/$gene\" target=\"_blank\">$symbol</a></td>
												 <td>$LR</td>
												 <td>$GR</td>
												 <td>$degree</td>
												 <td>$C11</td> <!-- clustering coefficient -->
												 <td>$C21</td>  <!-- betweeness -->
												 <td>$C41</a></td> <!-- page rank -->
												 <td>$C31</a></td>  <!-- closeness -->
												 <td>$C_1</td> <!-- overall global centrality -->
												 <td>$FC_1</td>
										</tr>
			";
	}
	else
	{
		$string =  "							<tr class=\"\">
												 <td class=\"lft\"><a href = \"http://www.ncbi.nlm.nih.gov/gene/$gene\" target=\"_blank\">$symbol</a></td>
												 <td>$LR</td>
												 <td>$GR</td>
												 <td>$degree</td>
												 <td>$C11</td> <!-- clustering coefficient -->
												 <td>$C21</td>  <!-- betweeness -->
												 <td>$C41</a></td> <!-- page rank -->
												 <td>$C31</a></td>  <!-- closeness -->
												 <td>$C_1</td> <!-- overall global centrality -->
										</tr>
			";
	}
	print WR "$string";
	close WR;
} 

sub print_end {
	($path,$n) = @_;
	my $F =  "$path"."/html/"."$n".".html"; ## this is the filename containg details of each pathway
	open(WR, ">>$F") or print ("Cannot write html file for pathway in the end: $!");
	my $string = "</tbody></table></p></body></html>";
	#print "$string";
	print WR "$string";
	close WR;
}

sub intersect(\@\@) {
	my %e = map { $_ => undef } @{$_[0]};
	return grep { exists( $e{$_} ) } @{$_[1]};
}

sub convert_sif_to_JS {
($path,$pathway,$name,$C)=@_;
%Convert = %$C;
$sif = "$path/html/$name.sif.txt";
$tempfc = "$path/html/temp.fc";
$network = "$path/html/$name"."_js.html";

open(NERT,">$network") or print "Cnnot open net file to write $!";
## give each gene an ID ###
print NERT "
<html><head><meta http-equiv=\"Content-Type\" content=\"text/html; charset=UTF-8\">
<title>$pathway</title>";
print NERT '
<meta name="viewport" content="width=device-width, user-scalable=no, initial-scale=1, maximum-scale=1">
<script src="../src/jquery-2.0.3.min.js"></script>
<script src="../src/cytoscape.min.js"></script>
<style>
body {font-family: helvetica;font-size: 14px;}
#cy {
width: 100%;
height: 100%;
position: absolute;
left: 0;
top: 0;
z-index: 999;
}

h1 {
opacity: 0.5;
font-size: 1em;
}
</style>

<script>
$(function()
{
var cy = window.cy = cytoscape({
container: document.getElementById(\'cy\'),
layout: {
name: \'cose\',
idealEdgeLength: 100,
nodeOverlap: 20
},

style: [
{"selector":"node","style":{"width":"mapData(score, 0, 0.006769776522008331, 20, 60)","height":"mapData(score, 0, 0.006769776522008331, 20, 60)","content":"data(name)","font-size":"9px","text-valign":"center","text-halign":"center","background-color":"#555","text-outline-color":"","text-outline-width":"0px","color":"black","overlay-padding":"6px","z-index":"10"}},
{"selector":"node:selected","style":{"border-width":"6px","border-color":"#AAD8FF","border-opacity":"0.5","background-color":"#77828C","text-outline-color":"#77828C"}},
{"selector":"edge","style":{"curve-style":"haystack","haystack-radius":"0.5","opacity":"0.4","line-color":"#bbb","width":"mapData(weight, 0, 1, 1, 8)","overlay-padding":"3px"}},
{"selector":"node.unhighlighted","style":{"opacity":"0.2"}},
{"selector":"edge.unhighlighted","style":{"opacity":"0.05"}},
{"selector":"node.highlighted","style":{"border-width":"6px","border-color":"#AAD8FF","border-opacity":"0.5","background-color":"#394855","text-outline-color":"#394855","shadow-blur":"12px","shadow-color":"#000","shadow-opacity":"0.8","shadow-offset-x":"0px","shadow-offset-y":"4px"}}, 
{"selector":"node[group=\"TF\"]","style":{"background-color":"#08088A","color":"#FFFFFF"}},
{"selector":"node[group=\"TFO\"]","style":{"background-color":"#FF0000"}},
{"selector":"node[group=\"TG\"]","style":{"background-color":"#00FF00"}},
{"selector":"node[group=\"OTHER\"]","style":{"background-color":"#424242","color":"#FFFFFF"}},
{"selector":"edge[group=\"UNK\"]","style":{"line-color":"#000000"}},
{"selector":"edge[group=\"RW\"]","style":{"line-color":"#FF8000"}}
],
elements: 
[
';

open(FH, "$tempfc") or print ("unable to open TEMPfc file $!");
$id = 0;

while ($x = <FH>)
{
	@temp = split(" ",$x);

	$temp[0]=~s/\s|\n//;$temp[1]=~s/\s|\n//;$temp[2]=~s/\s|\n//; ## gene ; FC ; TF/TG/TFO/OTHER
	if (($temp[0] ne ""))
	{
	if (exists $Convert{$temp[0]}){$symbol = $Convert{$temp[0]};} 
	else {$symbol = $temp[0];}

	print NERT "{\"data\":{\"id\":\"$temp[0]\",\"name\":\"$symbol\",\"score\":$temp[1],\"group\":\"$temp[2]\"},\"group\":\"nodes\"},\n";
	$des{$temp[0]} = $id;
	$id++;
	}
}

close FH;
@array= ();
open(FH, "$sif") or print ("unable to open SIF file $!");
@array = <FH>;
$count=1;
$x= "";

foreach $x (@array)
{
	@temp = split(" ",$x);
	$temp[0]=~s/\s|\n//;$temp[1]=~s/\s|\n//;$temp[3]=~s/\s|\n//; $temp[2]=~s/\s|\n//; #$temp[3] => UNK or RW
	$from = $temp[0]; $to=$temp[1];
	$weight = $temp[2];
	if (($from ne "UNK") || ($from ne "")  || ($from ne "NA"))
	{
		print NERT "{\"data\":{\"source\":\"$from\",\"target\":\"$to\",\"weight\":$weight,\"group\":\"$temp[3]\"},\"group\":\"edges\"},\n";
	}
	$count++;
}

print NERT '
]
});
});
</script>
</head>
<body>
<h1>';print NERT "$pathway"; print NERT '</h1>
<div id="cy"><div style="position: absolute; z-index: 0; overflow: hidden; width: 1919px; height: 895px;"><canvas style="position: absolute; z-index: 3; width: 1919px; height: 895px;" data-id="layer0-selectbox" width="1919" height="895"></canvas><canvas style="position: absolute; z-index: 2; width: 1919px; height: 895px;" data-id="layer1-drag" width="1919" height="895"></canvas><canvas style="position: absolute; z-index: 1; width: 1919px; height: 895px;" data-id="layer2-node" width="1919" height="895"></canvas></div></div>
<div style="font-family: &#39;Helvetica Neue&#39;, Helvetica, sans-serif; font-style: normal; font-size: 12px; font-weight: normal; position: absolute; left: -9999px; top: -9999px; z-index: -1; visibility: hidden; pointer-events: none; padding: 0px; line-height: 1; white-space: normal;">DNA2</div></body></html>
';

close NERT;
}

sub remove {
	my($path,$FH,$ann) = @_;
	my %Final = %$FH;
	my %annotation  = %$ann;
	@names = keys %Final;
	$html = "$path/html";
	unlink("$path/DE.tsv");
	unlink("$path/centrality.txt");
	unlink("$path/centrality.txt");
	unlink("$path/GRN.ssv");
	unlink("$path/GS.ssv");
	unlink("$path/GS.ssv");
	unlink("$path/RWR.txt");
	unlink("$path/TF.txt");
	unlink("$path/TG.txt");
	unlink("$path/score.tsv");
	
	foreach $name (@names)
	{
		$n = $annotation{$name}; 
		$present{$n}++;
	}
	### read directory ###
	opendir (DIR, $html);
	while (my $f = readdir(DIR)) 
	{
		$f=~s/\s//gi;
		if ($f=~/(P[0-9]+).*$/)
		{
		$test = $1;
			if (exists $present{$test}) ## this file reprsent dysreulated pathway 
			{
				next;
			}
			else
			{
			unlink("$html/$f");
			}
		}
	}
}

sub pie_chart {
	my($path,$out_inside,$outside,$inside,$pcount,$nonDRfound) = @_;
	open(PIE,">$path/html/pie.html") or print "Unable to print pie.html $!";
	#my $total = $out_inside + $outside + $inside;
	print PIE "
	<!DOCTYPE HTML>
	<html>
	<head>
	<script src=\"../src/canvasjs.min.js\"></script>
	<script type=\"text/javascript\">

	window.onload = function () {
		var chart = new CanvasJS.Chart(\"chartContainer\", {
			theme: \"theme2\",//theme1
			title:{
				text: 'Pathway Count: $pcount'              
			},
			animationEnabled: false,   // change to true
			data: [              
			{
				type: \"pie\",
				dataPoints: [
					{ label: \"Inside-Outside Rewiring\",  y: $out_inside  },
					{ label: \"Outside Rewiring\", y: $outside  },
					{ label: \"Inside Rewiring\", y: $inside  },
					{ label: \"Non-Differentially regulated\", y: $nonDRfound  },
				]
			}
			]
		});
		chart.render();
}
</script>
</head>
<body>
<div id=\"chartContainer\" style=\"height: 250px\"></div>
</body>
</html>
";
close PIE;
}

sub bar_chart {
	my($path,$out_inside,$outside) = @_;
	open(BAR,">$path/html/bar.html") or print "Unable to print bar.html $!";
	#my $total = $out_inside + $outside + $inside;
	print BAR "
	<!DOCTYPE HTML>
	<html>
	<head>
	<script src=\"../src/canvasjs.min.js\"></script>
	<script type=\"text/javascript\">

	window.onload = function () {
		var chart = new CanvasJS.Chart(\"chartContainer\", {
			theme: \"theme2\",//theme1
			title:{
				text: ''              
			},
			animationEnabled: false,   // change to true
			data: [              
			{
				type: \"bar\",
				dataPoints: [
					{ label: \"Total altered pathways\",  y: $out_inside  },
					{ label: \"Differentially Regulated pathways\", y: $outside  },
				]
			}
			]
		});
		chart.render();
}
</script>
</head>
<body>
<div id=\"chartContainer\" style=\"height: 250px\"></div>
</body>
</html>
";
close BAR;
}

sub scatterPlot {
	my ($path,$DS,$DP) = @_;
	@DisPotential = @$DP;
	@DysScore = @$DS;
	$corrlation = "NA";
	$corrlation  = &correlation(\@DisPotential,\@DysScore);
	open(SCT,">$path/html/scatter.html") or print "Unable to print scatter.html $!";
print SCT "
	<!DOCTYPE HTML>
<html>
<head>  
	<script type=\"text/javascript\">
	window.onload = function () {
		var chart = new CanvasJS.Chart(\"chartContainer\",
		{
			title:{
				text: \"Pathway Differential Regulation vs Rewiring (Correlation: $corrlation)\",      
				fontFamily: \"arial black\",
				fontColor: \"DarkSlateGrey\"
			},
            animationEnabled: true,
			axisX: {
				title:\"Rewiring score\",
				titleFontFamily: \"arial\"

			},
			axisY:{
				title: \"Differential Regulation score\",
				titleFontFamily: \"arial\",
				titleFontSize: 12
			},

			data: [
			{        
				type: \"scatter\",  
				//toolTipContent: \"<span style='\"'color: {color};'\"'><strong>{name}</strong></span> <br/> <strong>Cost/ container</strong> {y} \$<br/> <strong>Dysregulation</strong> {x} \",
				dataPoints: [
";
for $i (0..$#DisPotential)
{
$x = $DisPotential[$i];
$y = $DysScore[$i];
	if (($x=~/[0-9]/g) && ($y=~/[0-9]/g))
	{
	print SCT "{x: $x, y: $y},\n";
	}
}
print SCT "				]
			}
			]
		});

chart.render();
}
</script>
<script type=\"text/javascript\" src=\"../src/canvasjs.min.js\"></script>
</head>
<body>
	<div id=\"chartContainer\" style=\"height: 300px; width: 68%;\" >
	</div>
</body>
</html>
	";
close SCT;
}

sub correlation {

   my ($a1,$a2) = @_;
	my $correl=numer2($a1,$a2);
	my $xcorrel=sprintf("%.4f",$correl);
   return($xcorrel);
}

sub numer2 { ### ss = sum of squared deviations to the mean

   my ($x1,$y1)=@_;
   my @x = @$x1;
   my @y = @$y1;
   my $sum = '0';
   for (my $i=1;$i<scalar(@x);++$i){
	 $sum1 += $x[$i] * $y[$i];
	 $sum2 += $x[$i];
	 $sum3 += $y[$i];
	 $sum4 += $x[$i]**2;
	 $sum5 += $y[$i]**2;
   }
   my $fst = ($#x * $sum1) - ($sum2 * $sum3);
   my $third  = ($#x * $sum4) - ($sum2)**2;
   my $forth  = ($#y * $sum5) - ($sum3)**2;
   my $secon = (sqrt($third)) * (sqrt($forth)); 
   my $cor = $fst/$secon;
   return $cor;
}

1;
