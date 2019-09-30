# APA
Altered Pathway Analyzer
http://bioinfo.icgeb.res.in/APA/

Kaushik A, Ali S, Gupta D. Altered Pathway Analyzer: A gene expression dataset analysis tool for identification and prioritization of differentially regulated and network rewired pathways. Sci Rep. 2017 Jan 13;7:40450. doi:10.1038/srep40450

Altered Pathway Analyzer (APA) is a cross-platform and standalone tool for analyzing gene
expression datasets to highlight significantly rewired pathways across case-vs-control
conditions. The Tool is designed to analyze human gene expression datasets (with Entrez ID);
however, the analysis can be performed on gene expression datasets of other species by using
appropriate flags and input files.
APA algorithm is unique prioritization algorithm that also uses gene-regulatory network to
identify transcriptionally dysregulated pathways. It, thus, uses altered Transcription Factor
(TF) and Target Gene (TG) relationship for prioritizing gene circuit rewired pathways. For
Human datasets, it requires Gene expression matrices (RNAseq read count or Normalized
gene expression values) with genes in Entrez ID format. For other species, please refer the
manual


  Prerequisite:

    OS: Windows, Mac or Linux Compilers:
    perl and R (version ≥ 3.1.2)
    Dependencies: R packages- “matrixStats”, “mefa4”, “dnet”, “SANTA”,“limma”, “Biobase”


  Running the sample dataset:

  The APA is packaged with sample dataset in which 17 control samples has wild-type p53 gene
  status and 33 samples have mutated p53 samples. The APA can be executed as:

    perl APA.pl -case sample/case.txt -control sample/control.txt -o sam_output


  By default KEGG [human] pathway database will be used, if user wants to change the pathway
  database then -P flag can be used:

    perl APA.pl -case sample/case.txt -control sample/control.txt -o sam_output -P PANTHER

  In this case, instead of KEGG, PANTHER database for human pathways will be used. For more
  such options, type:

    perl APA.pl

  Note: For RNAseq read count "-rnaseq" argument is required in APA command line.
