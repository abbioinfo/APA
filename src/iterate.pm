#!/usr/bin/perl 
package iterate;
use POSIX;

sub Iteration {
	my ($pa,$logFC_threshold,$iter,$p_val,$Gpath) = @_;
	%Orignal  = &orignal_result($pa); 
	@PIDs = keys %Orignal;
	
	for $ni (1..$iter)
	{
	print "$ni:";
	$path = "$pa/Iteration/$ni";
	%GRN= ();
	if ((-e "$path/GS.ssv") && (-e "$path/DE.tsv"))
	{
	($GRN,$GRN2) = &GetGRN($path,$Gpath);
	%GRN = %$GRN; %GRN2 = %$GRN2;
	@G = keys %GRN;
	($DE,$LOGFC) = &getFC_i($path,$p_val,$logFC_threshold);
	%DE = %$DE;
	%LOGFC = %$LOGFC;

	#print "Loading Gene Network properties...\n";
	
	open (GS,"$path/GS.ssv") or die ("Unable to open gene strength $!");
	while ($x=<GS>)
	{
	my($temp1, $temp2) = split(" ",$x);
	$temp1=~s/\s|\n//gi;$temp2=~s/\s|\n//gi;
	$Rs{$temp1} = $temp2;
	}
	close GS;
	foreach $name (@PIDs)
	{
	$name=~s/\s|\n//gi;
	$score = 0; $R1 = 0; $R2 = 0; $edgeC = 0; $actual = 0; $yes = 0; %det = ();$Totexprs = 0;$significant = (); $scoreContrib = 0; $contrib = 0;$justRW=0; $Rsco_SUB = 0;$Rsco = 0;
	%$Description = (); %RPP= (); %det = ();
	
	$fileD = "$path"."/$name".".sif";
		if (-e $fileD)
		{
		open(SNET, "$fileD") or print ("No file $fileD found $!");
		@array = ();
		@array = <SNET>;
		close SNET;
		$total = scalar(@array);
		$RWoNLY = 0; $tfRW = 0; $tfRW_outside = 0; $TFCount = 0;
			for $i (1..$total)
			{
				$x = $array[$i];
				($g1,$g2,$p) = split(" ",$x);
				$g1=~s/\s|\n//gi;$g2=~s/\s|\n//gi;$p=~s/\s|\n//gi;
				$det{$g1}{$g2} = 0; $det{$g2}{$g1} = 0;
				$Rsco += $p;
				if ($p >= 1.30103) ##RW => REWIRED (since p value is lod transformed with -log10 so p value of o.o5 becomes 1.30103
				{
					$Rsco_SUB += $p;
					
					$RPP{$name}{$g1} += $p; $RPP{$name}{$g2} += $p; ## rewiring score of each gene per pathway (rewiring per pathway); higer the score higher the number rewired edges 
					if ((exists $DE{$g1}) && ($GRN{$g1}{$g2} eq "TF-TF"))
					{
						## $g1 - $g2 is regulatory interaction
						$contrib += $LOGFC{$g1};
						$scoreContrib++;
						$Description{$name}{$g1} = "TF";
						$Description{$name}{$g2} = "TF";
						#$TFCount++;
					}
					if ((exists $DE{$g2}) && ($GRN{$g2}{$g1} eq "TF-TF"))
					{
						## $g1 - $g2 is regulatory interaction
						$contrib += $LOGFC{$g2};
						$scoreContrib++;
						$Description{$name}{$g1} = "TF";
						$Description{$name}{$g2} = "TF"; 
						#$TFCount++;
					}
					elsif ((exists $DE{$g2}) && ($GRN{$g2}{$g1} eq "TF-TG"))
					{
						## $g2 - $g1 is regulatory interaction
						$contrib += $LOGFC{$g2};
						$scoreContrib++;
						if ($Description{$name}{$g2} ne "TF"){$Description{$name}{$g2} = "TG"};
						$Description{$name}{$g1} = "TF";
						#$TFCount++;
					}
					elsif ((exists $DE{$g1}) && ($GRN{$g1}{$g2} eq "TF-TG"))
					{
						## $g2 - $g1 is regulatory interaction
						$contrib += $LOGFC{$g1};
						$scoreContrib++;
						if ($Description{$name}{$g1} ne "TF"){$Description{$name}{$g1} = "TG"};
						$Description{$name}{$g2} = "TF";
						#$TFCount++;
						
					}
					else
					{
						if (($g1 ne "") && ($g2 ne ""))
						{
							## NON-REGULATORY BUT REWIRED INTERACTION
							$RWoNLY = 1;
							if (($Description{$name}{$g2} ne "TF") && ($Description{$name}{$g2} ne "TG"))
							{
								$Description{$name}{$g2} = "OTHER";
								$significant++;
							}
							elsif (($Description{$name}{$g1} ne "TF") && ($Description{$name}{$g1} ne "TG"))
							{
								$Description{$name}{$g1} = "OTHER";
								$significant++;
							}
							elsif (($Description{$name}{$g1} ne "TF") && ($Description{$name}{$g1} ne "TF"))
							{
								$Description{$name}{$g1} = "OTHER";
								$Description{$name}{$g2} = "OTHER";
								$significant++;
							}
						}
						
					}
				}
			}
			
			foreach $G (keys %{$Description{$name}}) ## reading all genes one by one for the last processed pathway and attempting to finf rewired TF links that occured from outside the pathway 
			{
				@TFs = keys %{$GRN{$G}};
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
						$significant++;
						$scoreContrib++; ## counting the number of rewired interactions mediated by TFs
						$contrib += $LOGFC{$G};
						#$TFCount++;
						}
						elsif ($GRN{$G}{$myT} eq "TF-TG") ## a given pathway TG is mediated by other pathway TF 
						{
							if (exists $DE{$G})
							{ 
							$Description{$name}{$myT} = "TFO";
							if ($Description{$name}{$G} ne "TF"){$Description{$name}{$G} = "TG"};
							$contrib += $LOGFC{$G}; ## total pathway upregulation due to TF mediated rewiring 
							$scoreContrib++; ## counting the number of rewired interactions mediated by TFs
							$significant++; 
							}
						}
					}
				}
				@TFs = ();
			}
			@allG = keys %det;
			foreach $g (@allG)
			{  
				$Totexprs += $LOGFC{$g};
				$RWscore += $Rs{$g};
			}
				## score for this pathway
				$GRNfraction = $contrib/$Totexprs;
				$div = $Rsco_SUB/$Rsco;
				$Fscore  = (($GRNfraction * $scoreContrib) + $div);
				($Orig_RW,$Orig_GRN,$t1,$t2,$t3) = split(" " ,$Orignal{$name});
				$Orig_RW=~s/\s//gi; $Orig_GRN=~s/\s//gi;
				
				$Orig_RW1 = sprintf "%.2f", $Orig_RW;
				$Fscore1 = sprintf "%.2f", $Fscore;
				$Orig_GRN1 = sprintf "%.2f", $Orig_GRN;
				$GRNfraction1 = sprintf "%.2f", $GRNfraction;
				
				if (($Fscore1 > $Orig_RW1) && ($GRNfraction1 > $Orig_GRN1))
				{
					$result{$name}++;
					print ".";
				}
				else
				{
					$result{$name} += 0;
					print "|";
				}
			@allG = (); %det = (); 
			$RWscore = 0; 
			$GRNfraction = 0; 
			$contrib = 0; 
			$Totexprs = 0;
			$Rsco_SUB = 0; 
			$Rsco = 0; 
			$GRNfraction =0; 
			$scoreContrib = 0; 
			$significant = 0;
		}
		else
		{
			$result{$name} += 0;
			print "|";
		}
	}
	print "\n";
	}
}
	&print_start($pa);
	foreach $name (keys %result)
	{
	$pvalue = $result{$name}/$iter;
	($P,$rank,$DySc,$DR,$OP) = split(" " ,$Orignal{$name}); ## pathwayname Rank Dysregulationscore DRscore Oncopotential
		if ($pvalue <= 0.05) ## unadjusted p-value
		{
		$DySc=~s/\s//gi;
		$P=~s/\s//gi;
		$DR=~s/\s//gi;
		$OP=~s/\s//gi;
		
		&print_table($name,$pa,$P,$rank,$DySc,$DR,$pvalue);
		}
		print "$name: $P: $pvalue\n";
	}
	&print_stop($pa);
	print "Iterations completed !!! Please find $pa/Significant.html file
	
	Thank you for using DPA :)\n";
}

sub getFC_i {
my($setwd,$p_val,$FC) = @_;
	#$script  = "src/de.R"; 
	
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

sub GetGRN {
my($path,$Gpath) = @_;
open(FH,$Gpath) or die ("Unable to open GRN file: $!"); ## the orignal GRN backgound network 
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
open(FH,"$path/GRN.ssv"); ## processed GRN created with processpathway.pm 
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

sub orignal_result {
my($path) = @_;
my %hash = ();
open (RES, "$path/result.tsv") or die ("unable to open RESULT.tsv $!");
@data = <RES>;
	for $i (1..$#data)
	{
	@line = split(" ",$data[$i]);
	$line[0] =~s/\s//gi;$line[2] =~s/\s//gi;$line[6] =~s/\s|\n//gi;$line[5] =~s/\s|\n//gi;$line[2] =~s/\s|\n//gi;$line[3] =~s/\s|\n//gi; $line[1] =~s/\s//gi;
	$na = "$line[0]";
	$hash{$na} = "$line[1] $line[2] $line[4] $line[5] $line[3]"; ## pathwayname Rank Dysregulationscore DRscore Oncopotential
	}
return(%hash);
}

sub print_start {
	($path) = @_;
	open(WR, ">$path"."/Significant.html") or print ("Cannot write html file $!");
	open(READ, "$path"."/index.html") or print ("Cannot write html file $!");
	@index = <READ>; close READ;
	for $i (0..69)
	{
	print WR "$index[$i]";
	}
	
	
	print WR "
									  <th data-tsorter=\"link\"> </th>
									  <th data-tsorter=\"link\" onmouseover=\"Tip('Pathway Name')\" onmouseout=\"UnTip()\">Pathway</th>
									  <th data-tsorter=\"numeric\" onmouseover=\"Tip('Cumulative alteration and differential regulation score')\" onmouseout=\"UnTip()\">Rank</th>
									  <th data-tsorter=\"numeric\" onmouseover=\"Tip('Alteration (gene-gene rewiring) associated with pathway gene-set')\" onmouseout=\"UnTip()\">Dy Score</th>
									  <th data-tsorter=\"numeric\" onmouseover=\"Tip('Differential regulation of pathway gene-set')\" onmouseout=\"UnTip()\">DR Score</th>
									  <th data-tsorter=\"numeric\" onmouseover=\"Tip('P value significance')\" onmouseout=\"UnTip()\">P value</th>
									</tr>
								  </thead>
								  <tbody>";
	close WR;
}

sub print_table {
my ($n,$path,$P,$rank,$DySc,$DR,$pvalue) = @_;

	$rank1 = sprintf "%.3f", $rank;
	$DySc1 = sprintf "%.3f", $DySc;
	$DR1 = sprintf "%.3f", $DR;
	$F =  "html/"."$n".".html";
	open(WR, ">>$path"."/Significant.html") or print ("Cannot write html file $!");
	print WR "
			<tr class=\"\" onmouseover=\"Tip('Gene Count: $GeneCount')\" onmouseout=\"UnTip()\">
			<td><p style=\"border:none;color:$color;font-size:170%;\"><span>&#9679;</span></p></td>
			 <td class=\"lft\">
			 <a href=\"$F\">$P</a>
			 </td>
			 <td>$rank1</td>
			 <td>$DySc1</a></td>
			 <td>$DR1</td>
			 <td>$pvalue</td>
			</tr>
	";
	close WR;
}

sub print_stop {
	($path) = @_;
	open(WR, ">>$path"."/Significant.html") or print ("Cannot write html file $!");
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


1;
