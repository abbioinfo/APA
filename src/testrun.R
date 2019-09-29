## usage: Rscript Functions.R

args = commandArgs(TRUE)
errorfile = paste(args[1],"error.log",sep="/");
error = 0;
message("Checking dependencies ....")
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
	error = 1;
	}
}else
{
message("Package SANTA (recommended version >= 2.2.0) not found");
error = 2;
}


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
	error = 3;
	}
}else
{
message("Package matrixStats (recommended version >= 0.50.1) not found");
error = 4;
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
	error = 5;
	}
	if ("featureNames" %in% test)
	{
	message("..Ok")
	}else
	{
	message("featureNames..failed. Function from Biobase package (recommended version >= 2.26.0) not found")
	error = 6;
	}
	if ("exprs" %in% test)
	{
	message("..Ok")
	}else
	{
	message("exprs function from Biobase package (recommended version >= 2.26.0) not found")
	error = 7;
	}
}else
{
message("Package Biobase (recommended version >= 2.26.0) not found");
error = 8;
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
	error = 9;
	}
}else
{
message("Package mefa4 (recommended version >= 0.3) not found");
error = 10;
}


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
	error = 11;
	}
}else
{
message("Package dnet (recommended version >= 1.0.7) not found");
error = 12;
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
		error = 13;
		}
	}
}else
{
message("Package igraph (recommended version >= 1.0.1) not found");
error = 14;
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
	error = 15;
	}
}else
{
message("Package Matrix (recommended version >= 1.2-4) not found");
error = 16;
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
	error = 17;
	}
	if ("makeContrasts" %in% test)
	{
	message("..Ok")
	}else
	{
	message("makeContrasts..failed. Function from limma package (recommended version >= 3.22.7)  not found")
	error = 18;
	}
	if ("contrasts.fit" %in% test)
	{
	message("..Ok")
	}else
	{
	message("contrasts.fit..failed. Function from limma package (recommended version >= 3.22.7)  not found")
	error = 19;
	}
	if ("eBayes" %in% test)
	{
	message("..Ok")
	}else
	{
	message("eBayes..failed. Function from limma package (recommended version >= 3.22.7)  not found")
	error = 20;
	}
	if ("topTable" %in% test)
	{
	message("..Ok")
	}else
	{
	message("topTable function from limma package (recommended version >= 3.22.7)  not found")
	error = 21;
	}
}else
{
message("Package limma (recommended version >= 3.22.7) not found");
error = 22;
}

if (error > 0)
{
message("[Error] All dependencies are not satisfied. Please re-install all the dependencies :(")
write.table(paste("ERROR:",error),file=errorfile , row.names=F, quote=F);
}else
{
message("\n************************\nAll dependencies are installed. :)\n************************\n")
}





