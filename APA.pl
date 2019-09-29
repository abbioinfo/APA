use lib qw(src/);
require processPathway;
require iterate;
use Cwd 'abs_path';
use Cwd;


### DEAFULT VALUES ###
$cancer = "";
$control = "";
$path ="";
$s1 = "database/seeds.txt";
$G1 = "database/GRN.ssv";
$pname = "KEGG";
$FC = 0.00;		## Fold Change Threshold
$PCC = 0.1;		## Correlation Threshold
$pvalue = 0.05; ## P-VALUES
$method="RWR"; 	## RWR
$iter = NULL;		## Number of random iteration to perform to calculate p-value
$GC = 10;		## Gene count threshold
$status = "F";
$header = "T";  ## if the expression file contains column/sample names
$show = "F"; ## show DE of genes in sub-network prioritization
$Up_down = "F";
$Rscript = "";
$rnaseq = "FALSE";
########################
$pathway = "database/Pathways/KEGG.txt";

if ($#ARGV < 2)
{
	&help();
}
else
{
&arguments(@ARGV);
}

sub arguments
{
	print "Starting analysis...\n";
	my(@args) = @_;
	$i = 0;
	$err = 0;
	foreach $x (@args)
	{
		$x=~s/\s//gi;
		$y = $args[$i+1];
		$y=~s/\s//gi;
		$y=~s/\/$//gi;
	if ($x=~/^-/)
	{
		if ($x eq "-case")
		{
			
			if (-e $y)
			 {
				$cancer = $y;
				
				$err = 1;
			 }
			 else
			 {
				 print "[Error] Cannot find file $y";
				 &help();
				 exit;
			 }
		}
		
		elsif ($x eq "-control")
		{
			
			if (-e $y)
			{
				
					$control = $y;
					$err++;
			}
			else
			{
				 print "[Error] Cannot find file $y";
				 &help();
				 exit;
			}
		}
		
		elsif ($x eq "-o")
		{
			$y=~s/\/$//;
			if (-e $y) 
			{
				
				$path = $y;
				$err += 1;
			}
			else
			{
				mkdir $y, 0777;
				$path = $y;
				$err++;
			}
		}
		
		########### Optional Arguments ############
		
		elsif ($x eq "-s") ## search for seed ID(s)
		{
			$s1 = $y;
		}
		elsif ($x eq "-rnaseq") ## if input is read count matrices
		{
			$rnaseq = "TRUE";
		}
		elsif ($x eq "-grn")  ## checking gene regulatory network 
		{
			$G1 = $y;
		}
		
		elsif ($x eq "-P")    ## pathway database
		{
		($pname,$pathway) = &pathway_checker($y);
		}
		
		elsif ($x eq "-R")    ## path to Rscript
		{
			$value = "";
			$Rscript  = &Dir_checker($y);
		}
		
		elsif ($x eq "-GC") ## gene count threshold
		{
			
			if ($y > 0)
			{
				$GC = $y;
			}
			else
			{
				if ($GC=~/[0-9]/gi)
				{
					print "Gene count cannot be less than 3 :(\n";
				}
				else
				{
					print "Gene count must be a numeric value :(\n";
				}
			exit;
			}
		}
		
		elsif ($x eq "-m") ## method for disease gene enrichment
		{
			
			if ($y == 1)
			{
				$method = "RWR";
			}
			elsif ($y == 2)
			{
				$method = "Knet";
			}
			else
			{
				print "[Error] Unknown method: $y. Please select either 1 or 2 :(\n";
				exit;
			}
		}
		
		elsif ($x eq "-r") ## correlation threshold
		{
			
			if (($y < 0.90) && ($y > 0.00))
			{
				$PCC = $y;
			}
			else
			{
				print "[Error] Correlation threshold (r) can only be < 0.9 and > 0.0 :(\n";
				exit;
			}
		}
		
		elsif ($x eq "-f") ## foldchange threshold
		{
			
			if ($y > 0)
			{
				$FC = $y;
			}
			else
			{
				print "[Error] Fold change cannot be less than zero :(\n";
				exit;
			}
		}
		
		elsif ($x eq "-p") ## pvalue threshold
		{
			if($y=~/[a-z]/gi)
			{
			print "[Error] Unknown p value theshold $y :(\n";
			exit;
			}
			if (($y >= 0) && ($y < 1))
			{
				$pvalue = $y;
			}
			else
			{
				print "[Error] p value theshold must be ranging between 0-1 :(\n";
				exit;
			}
		}
		
		elsif ($x eq "-TF") ## List of TF ids .. works only with -grn
		{
			if (-e $y) 
			{
			$TransFactor  = $y;
			}
			else
			{
				print "[Error] File $y, containing list of TFs, not found :( \n";
				exit;
			}
		}
		
		elsif ($x eq "-A") ## do you want all the dysregulated pathways in final list or just differentially regulated
		{
			if (($y eq "T") || ($y eq "True") || ($y eq "TRUE"))
			{
			$status  = "T";
			}
			elsif (($y eq "F") || ($y eq "False") || ($y eq "FALSE"))
			{
			$status  = "F";
			}
			else
			{
				print "[Error] Unknown option $y. Please suggest either T or F (True/False)  :( \n";
				exit;
			}
		}
		
		elsif ($x eq "-D") ## do you want show gene DE 
		{
			if (($y eq "T") || ($y eq "True") || ($y eq "TRUE"))
			{
			$show  = "T";
			}
			elsif (($y eq "F") || ($y eq "False") || ($y eq "FALSE"))
			{
			$show  = "F";
			}
			else
			{
				print "[Error] Unknown option $y for -D. Please suggest either T or F (True/False)  :( \n";
				exit;
			}
		}
		
		elsif ($x eq "-H") ## header in expression file 
		{
			if (($y eq "T") || ($y eq "True") || ($y eq "TRUE"))
			{
			$header  = "T";
			}
			elsif (($y eq "F") || ($y eq "False") || ($y eq "FALSE"))
			{
			$header  = "F";
			}
			else
			{
				print "[Error] Unknown option $y. Please suggest either T or F (True/False)  :( \n";
				exit;
			}
		}
		
		elsif ($x eq "-down") ## header in expression file 
		{
			if (($y eq "T") || ($y eq "True") || ($y eq "TRUE"))
			{
			$Up_down  = "T";
			}
			elsif (($y eq "F") || ($y eq "False") || ($y eq "FALSE"))
			{
			$Up_down  = "F";
			}
			else
			{
				print "[Error] Unknown option $y. Please suggest either T or F (True/False)  :( \n";
				exit;
			}
		}
		
		elsif ($x eq "-niter")
		{
			if ($y > 9)
			{
			$iter = $y;
			}
			else
			{
			print "[Error] Number of iterations (niter) must be an integer greater than or equal to 10. Recommeded value 100:( \n";
			exit;
			}
		}
		
		else
		{
			print "[Error] Unknown option $x :(\n";
			&help();
			exit;
		}
	}
	$i++;
	}
	if ($Rscript eq "")
	{
		$Rscript = &get_R(); ## getting the path of Rscript
	}
	if ($err == 3)
	{
			&check_dependencies($Rscript,$path);
			
			##########disease###########
			%e1 = &exprs_checker($cancer,$path,"case.dat",$header);
			print "Disease file processed successfully :D\n";
			##########control###########
			%e2 = &exprs_checker($control,$path,"control.dat",$header);
			print "Control file processed successfully :D\n";
			############################
			%entrez = (%e1,%e2);
			#close LOG;
			$seed = &ID_checker($s1,$path,%entrez);
			if ($seed eq "")
			{
				print "[Error] None of the seed IDs you submitted is present in expression dataset. :(\n";
				exit;
			}
			############################
			$GRN = &GRN_checker($G1,$path,$TransFactor,%entrez);
			if ($GRN eq "")
			{
				print "[Error] None of the GRN IDs you submitted is present in expression dataset. :(\n";
				exit;
			}
			############################
			$html = "$path/html";
			mkdir $html, 0777;
		$absolute_Path = abs_path($0);
		if ($absolute_Path=~/(.*)\/APA.pl$/)
		{
			$AP = $1;
			### Copy javascript and iamges to src folder to $path #######
			mkdir "$path/src", 0777;
			opendir(DIR, "$AP/src/src/") or die ("[Error] Unable to open SRC folder $!");
			while ($file = readdir(DIR))
			{
				$file=~s/\s|\n//gi;
					open (RD, "$AP/src/src/$file") or die ("[Error] Unable to OPEN $file $!");
					@data = <RD>;
					close RD;
					if ($#data > 0)
					{
						open (FH, ">$path/src/$file") or die ("[Error] Unable to write SRC folder $file : $!");
						foreach $x (@data)
						{
						print FH "$x";
						}
						close FH;
					}
					@data = ();
			}
			$imgfile = "$AP/src/img/newtab.png";
			$outputfile = "$path/src/newtab.png";
			open (IMAGE, $imgfile) || print "[Warning] Can't Open $imgfile\n";
			binmode(IMAGE);
			open (OUTPUT, ">$outputfile") || die "Can't Open $outputfile\n";
			binmode(OUTPUT);
			print OUTPUT while (<IMAGE>);
			close(IMAGE);
			close (OUTPUT);
			
			$imgfile = "$AP/src/img/download.png";
			$outputfile = "$path/src/download.png";
			open (IMAGE, $imgfile) || print "[Warning] Can't Open $imgfile\n";
			binmode(IMAGE);
			open (OUTPUT, ">$outputfile") || die "[Error] Can't Open $outputfile\n";
			binmode(OUTPUT);
			print OUTPUT while (<IMAGE>);
			close(IMAGE);
			close (OUTPUT);
			
			$imgfile = "$AP/src/img/color_codes.png";
			$outputfile = "$path/src/color_codes.png";
			open (IMAGE, $imgfile) || print "[Warning] Can't Open $imgfile\n";
			binmode(IMAGE);
			open (OUTPUT, ">$outputfile") || die "[Error] Can't Open $outputfile\n";
			binmode(OUTPUT);
			print OUTPUT while (<IMAGE>);
			close(IMAGE);
			close (OUTPUT);
			
	print "
 **********Input parameters ********
 Main Directory: $path/
 Pathway source: $pname ($pathway)
 RNAseq redcount: $rnaseq
 Gene Regulatory source: $GRN
 Seed ID source: $seed
 Disease gene prediction: $method
 Gene count threshold: $GC
 Correlation threshold: $PCC
 Fold change threshold: $FC
 p value threshold: $pvalue
 Only Differential Regulation: $status
 R path: $Rscript
 Sub-network DE: $show
 Search Down-regulated genes: $Up_down
 Number of iterations: $iter
 ***********************************\n";
 
 open (LOG, ">>$path/Log.txt") or print ("[Warning] Cannot write Log.txt at $path: $!");
 	print LOG "
 **********Input parameters ********
 Main Directory: $path/
 Pathway source: $pname ($pathway)
 RNAseq redcount: $rnaseq
 Gene Regulatory source: $GRN
 Seed ID source: $seed
 Disease gene prediction: $method
 Gene count threshold: $GC
 Correlation threshold: $PCC
 Fold change threshold: $FC
 p value threshold: $pvalue
 Only Differential Regulation: $status
 R path: $Rscript
 Sub-network DE: $show
 Search Down-regulated genes: $Up_down
 Number of iterations: $iter
 ***********************************\n";
 &start($path,$pvalue,$FC,$PCC,$GC,$cancer,$control,$status,$Rscript,$GRN,$pname,$seed,$pathway,$Up_down,$iter,$rnaseq);
		###############################
		&processPathway::analyze($path,$GRN,$pathway,$pvalue,$FC,$PCC,$method,$iter,$GC,$Rscript,$seed,$status,$header,$show,$Up_down,$rnaseq);
		###############################
			my $OS = $^O;
			if ($OS=~/linux/gi)
			{
				if (-e "/usr/bin/google-chrome")
				{
					`usr/bin/google-chrome $path/index.html`;
				}
			}
		}
	}
	else
	{
		if($cancer eq "")
		{
			print "[Error] Missing argument for disease expression file :( !!! Try -case \n";
		}
		if($control eq "")
		{
			print "[Error] Missing argument for control expression file :( !!! Try -control \n";
		}
		if($path eq "")
		{
			print "[Error] Missing argument for output directory :( !!! Try -o \n";
		}
		&help();
	}

}

sub pathway_checker
{
my($file)=@_;
$good  = 0;
if ($file eq "KEGG")
{
	$pname = "KEGG";
	$pathway = "database/Pathways/KEGG.txt";
}
elsif($file eq "NCI")
{
	$pname = "NCI";
	$pathway = "database/Pathways/NCI.txt";
}
elsif($file eq "Panther")
{
	$pname = "Panther";
	$pathway = "database/Pathways/Panther.txt";
}
elsif($file eq "Reactome")
{
	$pname = "Reactome";
	$pathway = "database/Pathways/Reactome.txt";
}
elsif($file eq "MsigDB_C2")
{
	$pname = "MsigDB_C2";
	$pathway = "database/Pathways/MsigDB_C2.txt";
}
elsif($file eq "ALL")
{
	$pname = "ALL";
	$pathway = "database/Pathways/ALL.txt";
}
else
{
	open (WR, "$file") or die ("[Error] unable to open $file: $!");
	$pathway = $file;
	@pathways1 = <WR>;
	close WR;
	$pname = $file;
	foreach $x(@pathways1)
	{
		if ($x=~/\|/)
		{
			@pathway = split(/\|/,$x);
			if ($x=~/\t/)
			{
				@genes = split("\t",$pathway[1]);
				if ($#genes > 10)
				{
					$good++;
				}
			}
			else
			{
				print "Pathway database not in correct format. Not in tab delimited format :( \n";
			}
		}
		else
		{
		print "[Error] Pathway database not in correct format. Not in bar delimited format :( \n";
		exit;
		}
	}
	@pathways1 = ();
	if ($good > 1)
	{
		print "Pathway database selected successfuly with $good good pathways :) \n";
	}
	else
	{
		print " [Error] Pathway database not selected. All bad pathways, with less than 9 genes :( \n";
		exit;
	}
}
return($pname,$pathway);
}

sub Dir_checker
{
	my($temp2) = @_;
	$OS = $^O;
	$temp2=~s/\/$//gi; ## removing the "/" sign at the end
		if ($OS=~/linux/gi)
		{
			 if (-e "$temp2/Rscript" )
			 {
				$R_script = "$temp2/Rscript";
			 }
			 else
			 {
			 print "[Error] Rscript not found at $temp2 :(\n";
			 exit;
			 }
		}
		elsif ($OS=~/darwinkl/gi) ## mac
		{
			my $mpath = ""; ##default path for mac
			if(-e "$temp2/Rscript")
			{
				$R_script = "$temp2/Rscript";
			}
			else
			{
			print "Rscript not found at $temp2/Rscript and $mpath :( . Please install R\n";
			exit;
			}
		}
		elsif ($OS=~/win1/gi)
		{
			 if (-e "$temp2/Rscript.exe" )
			 {
				$R_script = "$temp2/Rscript.exe";
			 }
			 else
			 {
				print "[Error] Rscript not found at $temp2 :( . Please install R\n";
				exit;
			 }
		}
		
		else
		{
			if (-e "$temp2/Rscript")
			{
			$R_script = "$temp2/Rscript";
			}
			else
			{
			print "[Error] Unknown OS. Please specify the directory name having R binaries :(\n Try -R \n";
			exit;
			}
		}
return($R_script);
}

sub ID_checker
{
	my($file,$path,%entrez) = @_;
	my $seed_err = 0;
	open (WR, "$file") or die ("[Error] unable to open $file: $!");
	open (LOG, ">>$path/Log.txt") or print ("[Warning] cannot write Log.txt at $path: $!");
	while($x=<WR>)
	{
		$x=~s/\s|\n//;
		if (!exists $entrez{$x})
		{
			$err++;
			$seed_err++;
			print LOG "Seed_ID_Unknown:$x\n";
		}
		else
		{
			push(@seeds,$x);
			
		}
	}
	close LOG;
	if (scalar(@seeds) > 0)
	{
		if ($seed_err > 0)
		{
		print "Total Seeds IDs $seed_err do not match existing gene database. List of IDs are available in $path/Log.txt. :| \n";
		}
		$seed = $file;
	}
	else
	{
		print "Seed entrez ID not found in entrez gene database :( \n. ";
	}
return($seed);
}

sub GRN_checker
{
	my($file,$path,$TransFactor,%entrez) = @_;
	$count = 0;
	open (WR, "$file") or die ("[Error] unable to open $file: $!");
	open (LOG, ">>$path/Log.txt") or print ("[Warning] Cannot write Log.txt at $path: $!");
	if ($file ne "database/GRN.ssv")
	{
		open (GRN, ">$path/GRN_orig.ssv") or die ("[Error] Cannot write regulatory network at $path: $! :( \n");
		my @dataset = <WR>;
		$sep = &find_separator($dataset[0],$file);
		if ($TransFactor ne "")
		{
			########### its is required that user either write type of reg. interaction edge file. e.g. TF-TG or TG-TF else specify a list of TF IDs with -TF option
			open(DB, "$TransFactor") or die ("[Error] Cannot read TF file $TransFactor : $! :(\n");
			while($tf1 = <DB>)
			{
				$tf1=~s/\s|\n//gi;
				$myt{$tf1}++;
			}
			############
		}
		foreach $x (@dataset)
		{
			@arr = split($sep,$x);
			$arr[0]=~s/\s|\n//gi;$arr[1]=~s/\s|\n//gi;$arr[2]=~s/\s|\n//gi;
			
			if ((exists $entrez{$arr[0]}) && (exists $entrez{$arr[1]}))
			{
				if ($arr[2] eq "TF-TG")
				{
				print GRN "$arr[0] $arr[1] $arr[2]\n";
				$TF{$arr[0]}++; $TG{$arr[1]}++;
				$count++;
				}
				elsif ($arr[2] eq "TG-TF")
				{
				print GRN "$arr[1] $arr[2] TF-TG\n";
				$TF{$arr[0]}++; $TF{$arr[1]}++;
				$count++;
				}
				elsif ($arr[2] eq "TG-TF")
				{
				print GRN "$arr[1] $arr[2] TF-TG\n";
				$TF{$arr[1]}++; $TG{$arr[0]}++;
				$count++;
				}
				else
				{
					if (($#arr < 3) && ($TransFactor ne ""))
					{
						if ((exists $myt{$arr[0]}) && (!exists $myt{$arr[1]})) ## tf-tg
						{
						print GRN "$arr[0] $arr[1] TF-TG\n";
						$TF{$arr[0]}++; $TG{$arr[1]}++;
						$count++;
						}
						elsif ((exists $myt{$arr[1]}) && (!exists $myt{$arr[0]})) ##tg-tf
						{
						print GRN "$arr[1] $arr[0] TF-TG\n";
						$TF{$arr[1]}++; $TG{$arr[0]}++;
						$count++;
						}
						elsif((exists $myt{$arr[1]}) && (exists $myt{$arr[0]})) ## tf-tf
						{
						print GRN "$arr[1] $arr[0] TF-TF\n";
						$TF{$arr[0]}++; $TF{$arr[1]}++;
						$count++;
						}
						else ## i dont know 
						{
						print LOG "No_TF_for_$arr[0]-$arr[1]\n";
						}
					}
					else
					{
					print LOG "Ignored_GRN_Edge:$x";
					}
				}
				
			}
			my $GRN = "$path/GRN_orig.ssv";
		}
		if ($count != 0)
		{
			open (TF, ">$path/TF.txt") or die ("[Error] Cannot write TF.txt at $path: $!");
			open (TG, ">$path/TG.txt") or die ("[Error] Cannot write TG.txt at $path: $!");
			print "Total $count GRN edges found in expression datasets :| \n";
			foreach $k (keys %TF)
			{
				print TF "$k\n";
			}
			foreach $k (keys %TG)
			{
				print TG "$k\n";
			} 
		}
		else
		{
			print "[Error] No edge found in gene expression datasets :( \n";
			exit;
		}
	}
	else
	{
		my $GRN = $file;
		while($x=<WR>)
		{
			@arr = split(" ",$x);
			$arr[0]=~s/\s|\n//gi;$arr[1]=~s/\s|\n//gi;$arr[2]=~s/\s|\n//gi;
			if ((exists $entrez{$arr[0]}) && (exists $entrez{$arr[1]}))
			{
				$count++;
			}
		}
		if ($count != 0)
		{
			print "Total $count GRN edges found in expression datasets :| \n";
		}
		else
		{
			print "No edge found in gene expression datasets :( \n";
			exit;
		}
	}
return($file);
close LOG;
}

sub exprs_checker
{
	my($file,$path,$type,$head) = @_;
	##############
	%e = ();
	%entrez = ();
	open(F,"database/genes") or die("[Error] Unable to open entrez gene list database $!");
	while($x=<F>)
	{
		$x=~s/\s|\n//gi;
		$entrez{$x}++;
	}
	close F;
	##############
	
	open (W, $file) or die ("unable to open $file: $!");
	open (WR1, ">$path/$type") or die ("[Error] Cannot write $type at $path: $!");
	@dataset = <W>; 
	$sep = find_separator($dataset[0],$file);
	@header = split($sep,$dataset[0]);
	@line = split($sep,$dataset[1]);
	if ($head eq "T")
	{
	 if (($#header != $#line) && ($#header-1 != ($#line)) && ($#header != ($#line -1)))
	 {
		print "$#header != $#line \n ";
		print "[Error] Please check the header/sample Name of each column in disease expression file. :(\n Number of column names does not matches with number of expression values\n";
		exit;
	 }
	}
	else
	{
		foreach $mr (@header)
		{
			if ($mr=~/[a-z]|,|\<|\>|\?|\$|\*|\^|\%|\$|\+|\=|\-|\/|\\|\{|\}|\[|\]|_|"|'|:|;|~|\(|\)/gi)
			{
				print "[Error] Header or unknown character present in line 1 of expression dataset. :( \nIf column name is present in expression dataset then use -H T";
				exit;
			}
		}
	}
	$i = 1;
	$notfound = 0;
	$dataset[0]=~s/\n//gi;
	foreach $h (@header)
	{
		$h=~s/\s//gi;
		print WR1 "$h ";
	}
	print WR1 "\n";
	open (LOG, ">>$path/Log.txt") or print ("[Warning] Cannot write Log.txt at $path: $!");
	foreach $x (@dataset)
	{
		if ($i > 1)
		{
			@temp = split($sep, $x);
			if (scalar(@temp) != scalar(@line))
			{
			print "[Error] Please check line $i of $file. the number of 'Columns' are not equal to the number of 'Column Names' :(";
			exit;
			}
			else
			{
			$string = "";
			$temp[0]=~s/\s//gi;
				if (exists $entrez{$temp[0]})
				{
					$e{$temp[0]}++;
					foreach $y (@temp){$y=~s/\s|\n//gi; $y=~s/$sep//gi; $string .="$y ";}
					$string=~s/\s$//gi;
					print WR1 "$string\n";
				}
				else
				{
					$notfound++;
					print LOG "Entrez_ID_Unknown:$temp[0]\n";
				}
			}
		}
	$i++;
	}
	close LOG;
	
	$got = $notfound/scalar(@dataset);
	if (($got > 0.50) && ($notfound > 3000))
	{
	print "[WARNING] More than 50% of entrez IDs are either obsolete or wrong. List of IDs are available in $path/Log.txt. :| \n";
	}
	elsif (($got > 0.50) && ($notfound < 1000))
	{
	print "[Error] Most of entrez IDs are either obsolete or wrong. List of IDs are available in $path/Log.txt. :( \n";
	exit;
	}
	if ($got > 0)
	{
	print "[WARNING] Few entrez IDs are either obsolete or wrong. List of IDs are available in $path/Log.txt. :| \n";
	}
	close W; close WR1;
return(%e);
}

sub help
{
print "
          Altered Pathway Analyzer (version 3.0)
__________________________________________________________________

Usage: perl APA.pl -case [path] -control [path] -o [path] [Options]
_____________________________________________________________________

[Mandatory]
===========
-case	: Path to gene expression dataset of query samples (recommended sample size: 10)
-control : Path to gene expression dataset of control samples (recommended sample size: 10)
-o	: Path to directory where all results will be saved.

[Options]
=========
-rnaseq	: The input gene expression matrices have RNAseq generated read count values.
-s	: Path to text file containing entrez IDs of seed genes [Default: database/seeds.txt (Human)].
-P	: Select one of precomplied pathway gene sets: KEGG [Default]/ NCI / MsigDB_C2 / PANTHER / Reactome / ALL (Human)
								or
	  Path to text file containing pathway name and corresponding gene set in a specified format [Default: database/Pathways.All].
-grn	: Path to text file containing background regulatory edges in space delimited format [Default: database/GRN.ssv (Human)].
-R	: Path to directory having R executables- Rscript [Default: OS specific; not defined for iOS]
-GC	: Gene count threshold [Default: 10]
-TF	: List of transcription factor entrez IDs. Only valid with -grn [Default: database/TF.txt (Human)]
-m	: Method by which disease gene is to be predicted in disease network. 1/2
		1 - Random walk by restart algorithm [Default]
		2 - Knet algorithm [Warning: Slow]
-r	: Pearson correlation coefficient threshold [Default 0.1]
-p	: P value threshold [Default: 0.05]
-H	: [T/F] Whether column name or header is present in expression dataset or not. [Default: T]
-D	: [T/F] Whether to show gene differential expression for sub-network gene prioritization. [Default: F]
-A	: Whether to predict only differential regulation among altered pathways. T/F [Default: F]
-down	: [T/F] Whether to include downregulated genes. T/F [Default: F]
-niter	: [Int] Number of iterations for resampling and p-value calculation. Must be larger than or equal to 10. Recommeded value 100 [Default off]

e.g.
perl APA.pl -rnaseq -case case.htseq.count -control control.htseq.count -o out
perl APA.pl -case sample/cancer.dat -control sample/control.dat -o out/ -P KEGG -m 2 -r 0.3 
perl APA.pl -case sample/cancer.dat -control sample/control.dat -o out/ -P ALL -r 0.3 -s database/seeds.txt
perl APA.pl -case sample/cancer.dat -control sample/control.dat -o out/ -P database/Pathways.All -r 0.3 -s database/seeds.txt -niter 100

";
}

sub get_R
{
	my $OS = $^O;
		if ($OS=~/linux/gi)
		{
			
			 if (-e "/usr/local/bin/Rscript")
			 {
				$R_script = "/usr/local/bin/Rscript";
			 }
			 elsif (-e "/usr/bin/Rscript")
			 {
				$R_script = "/usr/bin/Rscript";
			 }
			 else
			 {
			 print "[Error] Rscript not found at /usr/local/bin/ or /usr/bin/ directory. Please specify the directory name having R binaries :(\n Try -R \n";
			 exit;
			 }
		}
		elsif ($OS=~/win/gi)
		{
			 if (-e "C:/Program files/R/bin/Rscript.exe")
			 {
				$R_script = "C:/Program files/R/bin/Rscript.exe";
			 }
			 else
			 {
				print "[Error] Rscript not found at C:/Program files/R/bin/ directory.  Please specify the directory name having R binaries :(\n Try -R \n";
				exit;
			 }
		}
		elsif ($OS=~/darwinkl/gi) ## mac
		{
			my $mpath = "";
			if (-e "$mpath")
			{
				$R_script = "$mpath";
			}
			else
			{
			print "Rscript not found at $mpath";
			exit;
			}
		}
		else
		{
			print "[Error] Unknown OS. Please specify the directory name having R binaries :(\n Try -R \n";
			exit;
		}
return($R_script);
}

sub find_separator
{
	my ($line,$file) = @_;
	if ($line=~/\t/)
	{
		$sep = "\t";
	}
	elsif ($line=~/,/)
	{
		$sep = ",";
	}
	elsif ($line=~/\s/)
	{
		$sep = " ";
	}
	elsif ($line=~/\|/)
	{
		$sep = "|";
	}
	else
	{
		print "[Error] Unable to detect delimitor in file header $file :( \n";
		exit;
	}
return($sep);
}

sub start
 {
	($path,$pvalue,$FC,$PCC,$GC,$cancer,$control,$status,$Rscript,$GRN,$pname,$seed,$pathway,$ud,$iter,$rnaseq) = @_;

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
				<h1><span><i>APA Result</i></span></h1>
				<h2><b>A</b>ltered <b>P</b>athway <b>A</b>nalyzer</h2>
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
                        <p>Control expression matrix: <b>$control</b></p>
                        <hr width= 60% align=\"left\" COLOR=\"grey\" SIZE=\"1\">
						<p>Disease expression matrix: <b>$cancer</b></p>
						<h3>Options selected:</h3>
						<p>RNAseq read count: <b>$rnaseq</b></p>
						<p>P value threshold: <b>$pvalue</b></p>
						<p>Correlation thershold: <b>$PCC</b></p> 
						<p>Fold change threshold: <b>$FC</b></p>
						<p>Pathway source: <b>$pname </b>($pathway) </p>
						<p>Gene Regulatory source: <b>$GRN</b> </p>
						<p>Seed ID source:<b> $seed</b> </p>
						<p>Gene count threshold: <b>$GC</b></p>
						<p>Only Differential Regulation: <b>$status</b></p>
						<p>R path: <b>$Rscript</b></p>
						<p>Include Down-regulated genes: <b>$ud</b></p>
						<p>Number of iterations:<b> $iter</b></p>
				    </div>
				    <!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------->
			        <div class=\"content-2\">
						<iframe src=\"html/pie.html\" id=\"chartContainer\" style=\" width:60%; height:300px; border:none; margin:20;overflow:hidden;\" scrolling=\"no\">  </iframe>
						<iframe src=\"html/scatter.html\" id=\"chartContainer\" style=\" width:60%; height:300px; border:none; margin:20;overflow:hidden;\" scrolling=\"no\" align=\"center\" >  </iframe>
						<br>
						<iframe src=\"html/bar.html\" id=\"chartContainer\" style=\" width:60%; height:300px; border:none; margin:20;overflow:hidden;\" scrolling=\"no\">  </iframe>
						</div>
						<div class=\"content-3\"> 
						<table id=\"housing_table_1\" class=\"sortable\">
								  <caption> <a href=#> Result table. Click header for column-wise sorting </a></caption>
								  <thead>
									<tr>
									  <th data-tsorter=\"link\"> </th>
									  <th data-tsorter=\"link\" onmouseover=\"Tip('Pathway Name')\" onmouseout=\"UnTip()\">Pathway</th>
									  <th data-tsorter=\"numeric\" onmouseover=\"Tip('Cumulative dysregulation and differential regulation score')\" onmouseout=\"UnTip()\">Rank</th>
									  <th data-tsorter=\"numeric\" onmouseover=\"Tip('Dysregulation (gene-gene rewiring) associated with pathway gene-set')\" onmouseout=\"UnTip()\">Dy Score</th>
									  <th data-tsorter=\"numeric\" onmouseover=\"Tip('Differential regulation of pathway gene-set')\" onmouseout=\"UnTip()\">DR Score</th>
									  <th data-tsorter=\"numeric\" onmouseover=\"Tip('Context-specific disease gene enrichment')\" onmouseout=\"UnTip()\">DP Score</th>
									</tr>
								  </thead>
								  <tbody>
	";
	close WR;
}

sub check_dependencies
{
	unlink("$path/error.log");
	my($Rpath,$path) = @_;
	my $script = "src/testrun.R";
	`$Rpath $script $path`;
	open (RD, "$path/error.log");
	@error = <RD>;
	foreach $e (@error)
	{
		$e=~s/\s|\n//gi;
		if ($e eq "ERROR")
		{
			exit;
		}
	}
}
