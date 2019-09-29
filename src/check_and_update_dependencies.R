## usage: Rscript script_name.R T
args = commandArgs(TRUE)
error=0;
update_dnet = 0;
update_igraph = 0;
update_limma = 0;
update_santa = 0;
update_mefa4 = 0;
update_matrixStats = 0;
update_matrix = 0;

message("Checking dependencies ....")


#############DNET############
if ("dnet" %in% rownames(installed.packages()))
{
suppressMessages(library(dnet));
Funct  = lsf.str("package:dnet")
test = Funct[1:length(Funct)]
	if ("dRWR" %in% test)
	{
	message("..Ok")
	}else
	{
	message("dRWR function from dnet package (recommended version >= 1.0.7) not found")
	error=1;
update_dnet = 1;
	}
}else
{
message("Package dnet (recommended version >= 1.0.7) not found");
error=1;
update_dnet = 1;
}

################ igraph ###########################
if ("igraph" %in% rownames(installed.packages()))
{
suppressMessages(library(igraph));
Funct  = lsf.str("package:igraph")
test = Funct[1:length(Funct)]
depends = c("graph_from_adjacency_matrix","transitivity","degree","closeness","page_rank","betweenness","vcount","strength","delete_edge_attr");
	for(i in 1:length(depends))
	{
		if (depends[i] %in% test)
		{
		message("..Ok")
		}else
		{
		message(paste(depends[i],"..failed. Function from igraph package (recommended version >= 1.0.1) not found",sep=""))
		error=1;
update_igraph = 1;
		}
	}
}else
{
message("Package igraph (recommended version >= 1.0.1) not found");
error=1;
update_igraph = 1;
}

########### matrix Stats #############
if ("matrixStats" %in% rownames(installed.packages()))
{
suppressMessages(library(matrixStats));
Funct  = lsf.str("package:matrixStats")
test = Funct[1:length(Funct)]

	if ("rowSds" %in% test)
	{
	message("..Ok")
	}else
	{
	message("rowSds..failed. Function from matrixStats package (recommended version >= 0.50.1) not found")
	}
}else
{
message("Package matrixStats (recommended version >= 0.50.1) not found");
}
##################################################

#############santa############
if ("SANTA" %in% rownames(installed.packages()))
{
suppressMessages(library(SANTA));
Funct  = lsf.str("package:SANTA")
test = Funct[1:length(Funct)]
	if ("Knode" %in% test)
	{
	message("..Ok")
	}else
	{
	message("Knode..failed. Function from SANTA package (recommended version >= 2.2.0) not found")
	error=1;
update_santa = 1;
	}
}else
{
message("Package SANTA (recommended version >= 2.2.0) not found");
error=1;
update_santa = 1;
}


################# biobase ########################
if ("matrixStats" %in% rownames(installed.packages()))
{
suppressMessages(library(Biobase));
Funct  = lsf.str("package:Biobase")
test = Funct[1:length(Funct)]
	if ("ExpressionSet" %in% test)
	{
	message("..Ok")
	}else
	{
	message("ExpressionSet..failed. Function from Biobase package (recommended version >= 2.26.0) not found")
	error=1;
update_matrixStats = 1;
	}
	if ("featureNames" %in% test)
	{
	message("..Ok")
	}else
	{
	message("featureNames..failed. Function from Biobase package (recommended version >= 2.26.0) not found")
	error=1;
update_matrixStats = 1;
	}
	if ("exprs" %in% test)
	{
	message("..Ok")
	}else
	{
	message("exprs function from Biobase package (recommended version >= 2.26.0) not found")
	error=1;
update_matrixStats = 1;
	}
}else
{
message("Package Biobase (recommended version >= 2.26.0) not found");
error=1;
update_matrixStats = 1;
}

#############mefa4############
if ("mefa4" %in% rownames(installed.packages()))
{
suppressMessages(library(mefa4));
Funct  = lsf.str("package:mefa4")
test = Funct[1:length(Funct)]
	if ("Melt" %in% test)
	{
	message("..Ok")
	}else
	{
	message("Melt function from mefa4 package (recommended version >= 0.3) not found")
	error=1;
update_mefa4 = 1;
	}
}else
{
message("Package mefa4 (recommended version >= 0.3) not found");
error=1;
update_mefa4 = 1;
}


#############Matrix############
if ("Matrix" %in% rownames(installed.packages()))
{
suppressMessages(library(Matrix));
Funct  = lsf.str("package:Matrix")
test = Funct[1:length(Funct)]
	if ("Matrix" %in% test)
	{
	message("..Ok")
	}else
	{
	message("Matrix..failed. Function from Matrix package (recommended version >= 1.2-4) not found")
	error=1;
update_matrix = 1;
	}
}else
{
message("Package Matrix (recommended version >= 1.2-4) not found");
error=1;
update_matrix = 1;
}

################# LIMMA ########################
if ("limma" %in% rownames(installed.packages()))
{
suppressMessages(library(limma));
Funct  = lsf.str("package:limma")
test = Funct[1:length(Funct)]
	if ("lmFit" %in% test)
	{
	message("..Ok")
	}else
	{
	message("lmFit..failed. Function from limma package (recommended version >= 3.22.7)  not found")
	error=1;
update_limma = 1;
	}
	if ("makeContrasts" %in% test)
	{
	message("..Ok")
	}else
	{
	message("makeContrasts..failed. Function from limma package (recommended version >= 3.22.7)  not found")
	error=1;
update_limma = 1;
	}
	if ("contrasts.fit" %in% test)
	{
	message("..Ok")
	}else
	{
	message("contrasts.fit..failed. Function from limma package (recommended version >= 3.22.7)  not found")
	error=1;
update_limma = 1;
	}
	if ("eBayes" %in% test)
	{
	message("..Ok")
	}else
	{
	message("eBayes..failed. Function from limma package (recommended version >= 3.22.7)  not found")
	error=1;
update_limma = 1;
	}
	if ("topTable" %in% test)
	{
	message("..Ok")
	}else
	{
	message("topTable function from limma package (recommended version >= 3.22.7)  not found")
	error=1;
update_limma = 1;
	}
}else
{
message("Package limma (recommended version >= 3.22.7) not found");
error=1;
update_limma = 1;
}


############# Updating packages ###################
if (error == 1)
{
	if (args[1] == "TRUE")
	{
		if (update_matrix == 1)
		{
		install.packages("Matrix");
		}

		if (update_dnet ==1)
		{
		message("Trying to update package dnet");
		#source("https://bioconductor.org/biocLite.R")
		source("http://bioconductor.org/biocLite.R")
		biocLite("supraHex")
		biocLite("Biobase")
		biocLite("graph")
		biocLite("Rgraphviz")
		install.packages("dnet" , repos="http://cran.rstudio.com/")
		}

		if (update_igraph == 1)
		{
		message("Trying to update package igraph");
		source("http://bioconductor.org/biocLite.R")
		biocLite("graph")
		install.packages("igraph", repos="http://cran.rstudio.com/")
		}

		if(update_limma == 1)
		{
		message("Trying to update package limma");
		#source("https://bioconductor.org/biocLite.R")
		source("http://bioconductor.org/biocLite.R")
		biocLite("limma")
		}

		if(update_santa == 1)
		{
		message("Trying to update package SANTA");
		#source("https://bioconductor.org/biocLite.R")
		source("http://bioconductor.org/biocLite.R")
		biocLite("SANTA")
		}

		if(update_mefa4 == 1)
		{
		message("Trying to update package mefa4");
		install.packages("mefa4", repos="http://cran.rstudio.com/");
		}

		if (update_matrixStats == 1)
		{
		message("Trying to update package matrixStats");
		install.packages("matrixStats", repos="http://cran.rstudio.com/");
		}
	}
}else
{
message("All dependencies are satisfied :)");
}
###############################################





