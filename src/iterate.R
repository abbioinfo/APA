#############################################################
######### Warning: Please donot edit the script #############
#############################################################
args = commandArgs(TRUE)
control_file = paste(args[1],"control.dat",sep="/");
disease_file = paste(args[1],"case.dat",sep="/");
result  = paste(args[1],"result.tsv",sep="/");
Totaliteration = as.numeric(args[2]) 
TrueFalse = as.logical(args[3]); ## header: T/F
Up_down = as.logical(args[4]); ## up- and down-regulated: T/F
PPATH = args[5];

TF_path = paste(args[1],"TF.txt",sep="/");
TG_path = paste(args[1],"TG.txt",sep="/");

method = "pearson";
message(control_file);
message(disease_file);

message(c("Total Iterations...",Totaliteration));
message("Loading R packages...")
suppressMessages(library(igraph));
#suppressMessages(library(matrixStats));
suppressMessages(library("reshape2"));
suppressMessages(library(limma));
suppressMessages(library(Biobase));
suppressMessages(library("Matrix"));
suppressMessages(library(mefa4));

options(warn=-1)

RtoZ <- function (x) 
{
    0.5 * log((1 + x)/(1 - x))
}

DE <- function(data1, data2)
{
merged12 = merge(data1,data2,by="row.names")
rn = merged12[,1]
rownames(merged12) = rn;
Eset = ExpressionSet(assayData=as.matrix(merged12[-1]))
ID <- featureNames(Eset)

#Symbol <- getSYMBOL(ID, "org.Hs.eg.db")
tmp <- data.frame(ID=ID, Symbol=ID, stringsAsFactors = F)
# Clean up
tmp[tmp=="NA"] <- NA
fData(Eset) <- tmp
# Clean up
fvarLabels(Eset) <- make.names(fvarLabels(Eset))
sml=NULL
for(i in 1:dim(data1)[[2]]) { sml = append(sml,"G0")}
for(i in 1:dim(data2)[[2]]) { sml = append(sml,"G1")}
ex <- exprs(Eset)
# set up the data and proceed with analysis
fl <- as.factor(sml)
Eset$description <- fl
design <- model.matrix(~ description + 0, Eset)
colnames(design) <- levels(fl)
fit <- lmFit(Eset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
gene_list <- topTable(fit2, coef=1, number=1000000, sort.by="logFC")
write.table (gene_list, file = "DE.tsv", sep = "\t", row.names=F, quote=F);
return (gene_list)
}

EnrichFisher <- function(nN,nK,nG,ng){
#nN The total number of genes in the background distribution.
#nK The number of interesting genes in an interesting gene list.
#nG The total number of genes in the known functional gene set.
#ng The number of interesting genes in a known functional gene set.
       a=nK-ng;
       b=ng;
       c=nN-nK-nG+ng;
       d=nG-ng;
       P_value=fisher.test(matrix(c(a,b,c,d),2,2,byrow=T))$p.value;
       return(P_value);
}


message("Reading case datasets...")
df1= as.data.frame(read.table(file=disease_file, header=TrueFalse, row.names=1))
message("Reading control datasets...")
df2= as.data.frame(read.table(file=control_file, header=TrueFalse, row.names=1))
message("Reading result files...")
message(result);
df3= as.data.frame(read.table(file=result, header=T, row.names=1))
message("Getting pathway list...");
RIDs = rownames(df3) 
message("Getting pathway Database...");
pathways = as.matrix(read.table(PPATH, sep="|", header=F, row.names=NULL)); ## check "|"
message("Getting TF list...");
if (file.exists(TF_path))
{
message("Picking user defined regulatory edges");
TF = scan (TF_path, what="character");
}else
{
TF = scan ("database/TF.txt", what="character");
}

if (file.exists(TG_path))
{
message("Picking user defined regulatory edges");
TG = scan (TG_path, what="character");
}else
{
TG = scan ("database/TG.txt", what="character");
}


#P.names = rownames(pathways);
rn = c(1:dim(pathways)[1]);
rn = paste("P",rn,sep="")
rownames(pathways) = rn;
pathways = pathways[RIDs,] ## now only those pathways will be considered which are orignally found to be dysregulated ; his will speed up the analysis
dir.create(paste(args[1],"Iteration", sep="/"));
setwd(paste(args[1],"Iteration", sep="/"));

## creating a random matrix ###
for (iter in 1:Totaliteration)
{
dir.create(as.character(iter));
setwd(as.character(iter));

message(c("Randomizing datasets...",iter))
range.rand  = as.integer((50 *dim(df1)[[2]]) / 100) + 1;
Normal.rand <- sample(1:dim(df1)[[2]], range.rand, replace=T);
df1.temp = df1[,Normal.rand];
remaining = setdiff(c(1:dim(df1)[[2]]), Normal.rand)
df1.temp2 = df1[,remaining];

range.rand  = as.integer((50 *dim(df2)[[2]]) / 100) + 1;
Normal.rand <- sample(1:dim(df2)[[2]], range.rand, replace=T);
df2.temp = df2[,Normal.rand];
remaining = setdiff(c(1:dim(df2)[[2]]), Normal.rand)
df2.temp2 = df2[,remaining];
###########################################################
merged12 = merge(df1.temp,df2.temp,by="row.names")
rownames(merged12) = merged12[,1]
SAMPLE1 = merged12[-1]

merged12 = merge(df1.temp2,df2.temp2,by="row.names")
rownames(merged12) = merged12[,1]
SAMPLE2 = merged12[-1]
#########################################################
r1  = rownames(SAMPLE1);
r2  = rownames(SAMPLE2);
genes = intersect(r1,r2);
S1 = SAMPLE1[genes,]; 
S2 = SAMPLE2[genes,];


if (Up_down)
{
genelist = DE(S2,S1);
temp1  = subset(genelist, !(logFC < 0.0))
temp2  = genelist[genelist$logFC > 0.5 & genelist$adj.P.Val <0.05,]$ID
genes = c(temp1$ID,temp2);
}else
{
genelist = DE(S1,S2)
temp1  = subset(genelist, !(logFC < 0.0))
genes = temp1$ID
}


message("Creating and processing networks")

datExprs1 = SAMPLE1[genes,];
datExprs2 = SAMPLE2[genes,];

mat1 = cor(t(datExprs1), method=method);
mat2 = cor(t(datExprs2), method=method);

################ new ###########################################
mat1[mat1 == 1] <- 0.98; ## if correlation is exactlty one the make it little less otherwise Fishers Z tranformation generate Inf and thus final p value cannot be calculated for any of the value 
mat2[mat2 == 1] <- 0.98;
mat1[mat1 < 0.1] <- 0; ## removing non-significant or inversely proportional edges 
mat2[mat2 < 0.1] <- 0;

Z1 =RtoZ(mat1)
Z2 =RtoZ(mat2)

selectedGenes = rownames(mat1);

rm(mat1,mat2);


n1 = dim(df1)[[2]];
n2 = dim(df2)[[2]];

RW <- abs((Z1 - Z2)/sqrt(1/(n1 - 3) + 1/(n2 - 3))); 
rm(Z1,Z2,datExprs1,datExprs2,temp1)
m = mean(RW)
sdev = sd(RW)
Z = (RW - m)/sdev

message("Calculating P-value...")
P = 2*pnorm(-abs(Z))
rm(Z)
rownames(P) = selectedGenes;
colnames(P) = selectedGenes;


message("Pathway analysis...")
P_val_vector = vector();
Pathway_vector = vector();

		for (i in 1:dim(pathways)[[1]])
		{
		#NAME = paste("P",i,sep="");
		NAME = rownames(pathways)[i];
		geneList = as.character(unlist(strsplit(pathways[i,2], split="\t")))
		found = intersect(selectedGenes,geneList); #### number of genes that exists in pathway as well as dataset using intersect()
		message(paste("|--",NAME,sep=""))

			if (length(found) >9)
			{
				#### get the subnetwork from RW matrix ###
				subnet = RW[found,found];
				temp = -1 * log10(P[found,found])
				temp[temp == Inf] <- 16 ## all the P.values with zero will be converted to Inf .. so convert them to highest value.. although it is expected that very little number of entries have zero p value 
				subnet.pvalue = subnet;
				subnet.pvalue[subnet.pvalue > 0] <- 1;
				subnet.pvalue = subnet.pvalue * temp; ## multiplying 1 with log transformed p value  
				graphP = graph.adjacency(subnet.pvalue, mode=c("undirected"), weighted=TRUE, diag=FALSE)
				graph = graph.adjacency(subnet, mode=c("undirected"), weighted=TRUE, diag=FALSE)
				graph = graph.adjacency(temp, mode=c("undirected"), weighted=TRUE, diag=FALSE)
				##########write both the edges in single file Pn.sif 
				graph = set.edge.attribute(graph,"P",index=E(graph), E(graphP)$weight);
				total = length(selectedGenes);
				Term = length(geneList);
				inTerm  = vcount(graph); 
				common  = length(intersect(as.character(V(graph)$name),geneList));
				p_val = EnrichFisher(total,inTerm,Term,common);
				if (p_val <= 0.05) ## if there are significant number of connected genes among pathway's gene list
				{
				write.table(cbind(get.edgelist(graph),E(graph)$weight,E(graph)$P),file = paste(NAME,".sif",sep=""), sep=" ", row.names=F, quote=F)
				}
				P_val_vector[length(P_val_vector)+1]  = p_val ;
				Pathway_vector [length(Pathway_vector)+1] = NAME;
			}
		}

### removing pathways with non-significant number of connected genes ####
p_adj = p.adjust(P_val_vector);

for (i in 1:length(p_adj))
{
NAME = Pathway_vector[i];
	if (p_adj > 0.05)
	{
	file.remove(paste(NAME,".sif",sep=""));
	file.remove(paste(NAME,"_prop.csv",sep=""));
	}
}

###### get RW and p.value of all TF-TF and TF-TG interaction (only those with p <= 0.05) and write it in single text file  

########### This is the most Time consuming step: but most important ##############
#######################################################################
message("Calculating P-value...");
P  = (-1 * log10(P))
P[P < 1.30103]<-0; ### ; 1.30103 means p value is less than 0.05 as we have log transformed the p value
P[P==Inf]<- 16;
message("GS analysis...")
#SP = Matrix(P, sparse = TRUE) ## sparse matrix 
#rm(P)
#message(dim(SP))
ig1 = graph.adjacency(P, mode=c("undirected"), weighted=TRUE, diag=FALSE)

GS = strength(ig1, loop=F)
write.table(GS, file="GS.ssv", row.names=T,quote=F,sep=" "); ##
rm(ig1)
message("GRN analysis...")
TF_C = intersect(TF,selectedGenes)
TG_C = intersect(TG,selectedGenes)
GRN  = P[TF_C,TG_C];

write.table(Melt(GRN),"GRN.ssv" ,sep= " ", quote=F, row.names=F)  ## melting the sparse matrix; MELTING WITH MEFA4

setwd("..");
}

q("no");
