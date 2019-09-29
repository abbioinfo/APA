#############################################################
######### Warning: Please donot edit the script #############
#############################################################
args = commandArgs(TRUE)
control_file = paste(args[1],"control.dat",sep="/");
disease_file = paste(args[1],"cancer.dat",sep="/");
method = as.character(args[2]) ## bicor to be added 
Totaliteration = as.numeric(args[3]) ## bicor to be added 

message(c("...",Totaliteration));
message("Loading R packages...")
suppressMessages(library(igraph));
suppressMessages(library(matrixStats));
suppressMessages(library("reshape2"));
suppressMessages(library(limma));
suppressMessages(library(Biobase));
suppressMessages(library("Matrix"));
suppressMessages(library(reshape2));
suppressMessages(library(mefa4));
suppressMessages(library(picante));

options(warn=-1)

RtoZ <- function (x) 
{
    0.5 * log((1 + x)/(1 - x))
}

DE <- function(data1, data2, iter)
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

message("Reading datasets...")
df1= as.data.frame(read.table(file=disease_file, sep= " ",header=T, row.names=1))
X = transform(as.matrix(df1), SD=rowSds(as.matrix(df1), na.rm=TRUE)); ## GETTING THE SD OF EACH ROW I.E. GENE
names = rownames(df1); ## GETTING THE IDS OF ALL GENES 
SAMPLE1 = X[names[X$SD > 0.000001],]; ## GETTING THE GENEIDS AND THEIR EXPRESSION VALUES IN ALL SAMPLES IF THE SD OF EXPRESSION VALUES IS > 0
SAMPLE1$SD <- NULL
r1  = rownames(SAMPLE1);
r1 = rownames(df1)

################################## Start ###############################################

df2= as.data.frame(read.table(file=control_file, sep= " ",header=T, row.names=1))
message("Removing non-varible genes in each dataset...")
X = transform(as.matrix(df2), SD=rowSds(as.matrix(df2), na.rm=TRUE)); ## GETTING THE SD OF EACH ROW I.E. GENE
names = rownames(df2); ## GETTING THE IDS OF ALL GENES 
SAMPLE2 = X[names[X$SD > 0.000001],]; ## GETTING THE GENEIDS AND THEIR EXPRESSION VALUES IN ALL SAMPLES IF THE SD OF EXPRESSION VALUES IS > 0
SAMPLE2$SD <- NULL
r2  = rownames(SAMPLE2);
genes = intersect(r1,r2);
S1 = SAMPLE1[genes,]; 
S2 = SAMPLE2[genes,]; 
genelist = DE(S1,S2,iter)
 
temp1  = subset(genelist, !(logFC < 0.0))
genes = temp1$ID

message("Creating and processing networks")

datExprs1 = df1[genes,];
datExprs2 = df2[genes,];

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





######reading the result.tsv file and getting the list of pathways for which p value is to be calculated ####
### first column aka rowname is the serial number of the pathway used in pathways.all 
Result= as.data.frame(read.table(file=paste(args[1],"result.tsv",sep="/"), sep= "\t",header=T, row.names=1))
RIDs = rownames(Result) ## these are the pathway IDs which are found to be dysregulated ; for which dysregulation score is to be calculated 

#### Read pathway file and create a vector of genes ###

#change path of database

pathways = as.matrix(read.table("../database/Pathways.All", sep="|", header=F, row.names=NULL)); ## check "|"
TF = scan ("../database/TF.txt", what="character");
TG = scan ("../database/TG.txt", what="character");




rn = c(1:dim(pathways)[1]);
rn = paste("P",rn,sep="")
rownames(pathways) = rn;
pathways = pathways[RIDs,] ## now only those pathways will be considered which are orignally found to be dysregulated ; his will speed up the analysis

dir.create(paste(args[1],"Iteration", sep="/"));
setwd(paste(args[1],"Iteration", sep="/"));
## creating a random matrix ###
for (iter in 1:Totaliteration)
{
message("Calculating P-value...")

P = 2*pnorm(-abs(Z))
#P.2 = randomizeMatrix(P.1,null.model = "frequency",iterations = dim(df2)[1]*dim(df2)[2])
#P = randomizeMatrix(t(P.2),null.model = "richness",iterations = dim(df2)[1]*dim(df2)[2])
rm(P.2)
rownames(P) = selectedGenes;
colnames(P) = selectedGenes;

dir.create(as.character(iter));
setwd(as.character(iter));

message("Pathway analysis...")
P_val_vector = vector();
Pathway_vector = vector();
		for (i in 1:dim(pathways)[[1]])
		{
		NAME = paste("P",i,sep="");
		geneList = as.character(unlist(strsplit(pathways[i,2], split="\t")))
		found = intersect(selectedGenes,geneList); #### number of genes that exists in pathway as well as dataset using intersect()
		message(paste("|--",NAME,sep=""))

			if (length(found) >9)
			{
				#### get the subnetwork from RW matrix ###
				#subnet = RW[found,found];
				temp = -1 * log10(P[found,found])
				temp[temp == Inf] <- 16 ## all the P.values with zero will be converted to Inf .. so convert them to highest value.. although it is expected that very little number of entries have zero p value 
				#subnet.pvalue = subnet;
				#subnet.pvalue[subnet.pvalue > 0] <- 1;
				#subnet.pvalue = subnet.pvalue * temp; ## multiplying 1 with log transformed p value  
				#graphP = graph.adjacency(subnet.pvalue, mode=c("undirected"), weighted=TRUE, diag=FALSE)
				#graph = graph.adjacency(subnet, mode=c("undirected"), weighted=TRUE, diag=FALSE)
				graph = graph.adjacency(temp, mode=c("undirected"), weighted=TRUE, diag=FALSE)
				##########write both the edges in single file Pn.sif 
				#graph = set.edge.attribute(graph,"P",index=E(graph), E(graphP)$weight);
				total = length(selectedGenes);
				Term = length(geneList);
				inTerm  = vcount(graph); 
				common  = length(intersect(as.character(V(graph)$name),geneList));
				p_val = EnrichFisher(total,inTerm,Term,common);
				if (p_val <= 0.05) ## if there are significant number of connected genes among pathway's gene list
				{
				write.table(cbind(get.edgelist(graph),E(graph)$weight),file = paste(NAME,".sif",sep=""), sep=" ", row.names=F, quote=F)
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




setwd("../");

}

q("no");
