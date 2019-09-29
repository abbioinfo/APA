## usage: Rscript Functions.R
error = 0;

#############santa############
if ("SANTA" %in% rownames(installed.packages()))
{
suppressMessages(library(SANTA));
Funct  = lsf.str("package:SANTA")
test = Funct[1:length(Funct)]
	if ("Knode" %in% test)
	{
	message("Knode..OK")
	}else
	{
	message("Knode..failed. Function from SANTA package not found")
	error = 1;
	}
}else
{
message("Package SANTA (recommended version >= 2.2.0) not found");
error = 1;
}


if ("matrixStats" %in% rownames(installed.packages()))
{
suppressMessages(library(matrixStats));
Funct  = lsf.str("package:matrixStats")
test = Funct[1:length(Funct)]

	if ("rowSds" %in% test)
	{
	message("rowSds..OK")
	}else
	{
	message("rowSds..failed. Function from matrixStats package not found")
	error = 1;
	}
}else
{
message("Package matrixStats (recommended version >= 0.50.1) not found");
error = 1;
}

################# biobase ########################
if ("matrixStats" %in% rownames(installed.packages()))
{
suppressMessages(library(Biobase));
Funct  = lsf.str("package:Biobase")
test = Funct[1:length(Funct)]
	if ("ExpressionSet" %in% test)
	{
	message("ExpressionSet..OK")
	}else
	{
	message("ExpressionSet..failed. Function from Biobase package not found")
	error = 1;
	}
	if ("featureNames" %in% test)
	{
	message("featureNames..OK")
	}else
	{
	message("featureNames..failed. Function from Biobase package not found")
	error = 1;
	}
	if ("exprs" %in% test)
	{
	message("exprs..OK");
	}else
	{
	message("exprs function from Biobase package not found")
	error = 1;
	}
}else
{
message("Package Biobase (recommended version >= 2.26.0) not found");
error = 1;
}

#############mefa4############
if ("mefa4" %in% rownames(installed.packages()))
{
suppressMessages(library(mefa4));
Funct  = lsf.str("package:mefa4")
test = Funct[1:length(Funct)]
	if ("Melt" %in% test)
	{
	message("Melt..OK");
	}else
	{
	message("Melt function from mefa4 package not found")
	error = 1;
	}
}else
{
message("Package mefa4 (recommended version >= 0.3) not found");
error = 1;
}


#############DNET############
if ("dnet" %in% rownames(installed.packages()))
{
suppressMessages(library(dnet));
Funct  = lsf.str("package:dnet")
test = Funct[1:length(Funct)]
	if ("dRWR" %in% test)
	{
	message("dRWR..OK")
	}else
	{
	message("dRWR function from dnet package not found")
	error = 1;
	}
}else
{
message("Package dnet (recommended version >= 1.0.7) not found");
error = 1;
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
		message(paste(depends[i],"..OK",sep=""))
		}else
		{
		message(paste(depends[i],"..failed. Function from igraph package not found",sep=""))
		error = 1;
		}
	}
}else
{
message("Package igraph (recommended version >= 1.0.1) not found");
error = 1;
}

##############Matrix############
if ("Matrix" %in% rownames(installed.packages()))
{
suppressMessages(library(Matrix));
Funct  = lsf.str("package:Matrix")
test = Funct[1:length(Funct)]
	if ("Matrix" %in% test)
	{
	message("Matrix..OK");
	}else
	{
	message("Matrix..failed. Function from Matrix package not found")
	error = 1;
	}
}else
{
message("Package Matrix (recommended version >= 1.2-4) not found");
error = 1;
}

################# LIMMA ########################
if ("limma" %in% rownames(installed.packages()))
{
suppressMessages(library(limma));
Funct  = lsf.str("package:limma")
test = Funct[1:length(Funct)]
	if ("lmFit" %in% test)
	{
	message("lmfit..OK");
	}else
	{
	message("lmFit..failed. Function from limma package not found")
	error = 1;
	}
	if ("makeContrasts" %in% test)
	{
	message("makeContrasts..OK")
	}else
	{
	message("makeContrasts..failed. Function from limma package not found")
	error = 1;
	}
	if ("contrasts.fit" %in% test)
	{
	message("contrasts.fit..OK")
	}else
	{
	message("contrasts.fit..failed. Function from limma package not found")
	error = 1;
	}
	if ("eBayes" %in% test)
	{
	message("eBayes..OK");
	}else
	{
	message("eBayes..failed. Function from limma package not found")
	error = 1;
	}
	if ("topTable" %in% test)
	{
	message("topTable..OK");
	}else
	{
	message("topTable function from limma package not found")
	error = 1;
	}
}else
{
message("Package limma (recommended version >= 3.22.7) not found");
error = 1;
}

if (error == 1)
{
message("[Error] All dependencies are not staisfied. Please install all the dependencies :(")
}else
{
message("\n************************\nAll dependencies are installed. :)\n************************\n")
}






