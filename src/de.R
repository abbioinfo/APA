#############################################################
######### Warning: Please donot edit the script #############
#############################################################
args = commandArgs(TRUE)
control_file = paste(args[1],"control.dat",sep="/");
disease_file = paste(args[1],"case.dat",sep="/");
DEfile = paste(args[1],"DE.tsv",sep="/");
Netfile = paste(args[1],"GS.ssv",sep="/"); ## gene strength in rewired network 
out_file = paste(args[1],"RWR.txt",sep="/");
#node_file = paste("database","seeds.txt",sep="/");
#Rsession = paste(args[1],"Rsession.Rdata",sep="/");
result2 = paste(args[1],"centrality.txt",sep="/");
strInfo  = paste(args[1],"str.Info",sep="/");
igraphSession  = paste(args[1],"igraph.Rdata",sep="/");

GPATH  = args[2]; ## PATH TO BACKGROUND GENE REGULATORY NETWORK
PPATH = args[3]; ## PATH TO BACKGROUND PATHWAY GENE SET
SPATH = args[7]; ## PATH TO SEED GENE FILE

mtype  = as.character(args[4]); ## method to predict disease gene in case network
GC = as.integer(args[5]); ## gene count threshold
PCC = as.integer(args[6]); ##PCC threshold
TrueFalse = as.logical(args[8]); ## header: T/F
Up_down = as.logical(args[9]); ## up- and down-regulated: T/F
rnaseq = as.logical(args[10]); ## up- and down-regulated: T/F

message("Loading R packages...")

suppressMessages(library(matrixStats));
suppressMessages(library(Biobase));
suppressMessages(library(mefa4));
suppressMessages(library(dnet));
suppressMessages(library(limma));
suppressMessages(library(SANTA));

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
write.table (gene_list, file = DEfile, sep = "\t", row.names=F, quote=F);
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
message(disease_file)
message(control_file)

df1= as.data.frame(read.table(file=disease_file,header=TrueFalse, row.names=1))
df2= as.data.frame(read.table(file=control_file,header=TrueFalse, row.names=1))

if (rnaseq)
{
df1 = as.data.frame(voom(df1, normalize.method = "quantile", plot=F)$E)
df2 = as.data.frame(voom(df2, normalize.method = "quantile", plot=F)$E)
}

message("Removing non-variable genes in each dataset...")
X = transform(as.matrix(df1), SD=rowSds(as.matrix(df1), na.rm=TRUE)); ## GETTING THE SD OF EACH ROW I.E. GENE
names = rownames(df1); ## GETTING THE IDS OF ALL GENES 
SAMPLE1 = X[names[X$SD > 0.000001],]; ## GETTING THE GENEIDS AND THEIR EXPRESSION VALUES IN ALL SAMPLES IF THE SD OF EXPRESSION VALUES IS > 0
SAMPLE1$SD <- NULL
r1 = rownames(SAMPLE1);


X = transform(as.matrix(df2), SD=rowSds(as.matrix(df2), na.rm=TRUE)); ## GETTING THE SD OF EACH ROW I.E. GENE
names = rownames(df2); ## GETTING THE IDS OF ALL GENES 
SAMPLE2 = X[names[X$SD > 0.000001],]; ## GETTING THE GENEIDS AND THEIR EXPRESSION VALUES IN ALL SAMPLES IF THE SD OF EXPRESSION VALUES IS > 0
SAMPLE2$SD <- NULL
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
genelist = DE(S2,S1)
temp1  = subset(genelist, !(logFC < 0.0)) 
genes = temp1$ID 
}


message("Creating and processing networks")

datExprs1 = df1[genes,];
datExprs2 = df2[genes,];

mat1 = cor(t(datExprs1), method="pearson");
mat2 = cor(t(datExprs2), method="pearson");

################ new ###########################################
mat1[mat1 == 1] <- 0.98; ## if correlation is exactlty one the make it little less otherwise Fishers Z tranformation generate Inf and thus final p value cannot be calculated for any of the value 
mat2[mat2 == 1] <- 0.98;
mat1[mat1 < PCC] <- 0; ## removing non-significant or inversely proportional edges 
mat2[mat2 < PCC] <- 0;

selectedGenes = rownames(mat1);

Z1 =RtoZ(mat1)
Z2 =RtoZ(mat2)
rm(mat1,mat2);


n1 = dim(df1)[[2]];
n2 = dim(df2)[[2]];

RW <- abs((Z1 - Z2)/sqrt(1/(n1 - 3) + 1/(n2 - 3))); 
rm(Z1,Z2,df1,df2,datExprs1,datExprs2,temp1)
m = mean(RW)
sdev = sd(RW)
Z = (RW - m)/sdev
message("Calculating P-value...")
P =  2*pnorm(-abs(Z))
## like RW create a matrix of p.value 

rm(Z)
message("Pathway analysis...")
#### Read pathway file and create a vector of genes ###
pathways = as.matrix(read.table(PPATH, sep="|", header=T, row.names=NULL)); ## check "|"
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
		subnet = RW[found,found];
		
		temp = -1 * log10(P[found,found])
		temp[temp == Inf] <- 16 ## all the P.values with zero will be converted to Inf .. so convert them to highest value.. although it is expected that very little number of entries have zero p value 
		subnet.pvalue = subnet;
		subnet.pvalue[subnet.pvalue > 0] <- 1;
		subnet.pvalue = subnet.pvalue * temp; ## multiplying 1 with log transformed p value  
		
		graphP = graph_from_adjacency_matrix(subnet.pvalue, mode=c("undirected"), weighted=TRUE, diag=FALSE)
		graph = graph_from_adjacency_matrix(subnet, mode=c("undirected"), weighted=TRUE, diag=FALSE)
		c.l = transitivity(graph, type='local',weights=T);
		p = page_rank(graph, directed=FALSE)$vector;
		b = betweenness(graph);
		d = degree(graph);
		C3 = closeness(graph);
		FC = genelist[found,3]
		Dise = as.data.frame(cbind(c.l,d,p,b,C3,FC));
		
		colnames(Dise) = c("Transitivity","Degree","Page Rank", "Betweenes", "Closeness","Fold Change");
		total = length(selectedGenes);
		Term = length(geneList);
		inTerm  = vcount(graph);
		if (inTerm >= GC) ## more than desired number of connected genes in a pathway network
		{
		common  = length(intersect(as.character(V(graph)$name),geneList));
		p_val = EnrichFisher(total,inTerm,Term,common);
			if (p_val <= 0.05) ## if there are significant number of connected genes among pathway's gene list
			{
				out = paste(args[1],"html",paste(NAME,"_prop.csv",sep=""),sep="/");
				write.table(file = out, Dise, quote=F, sep = " "); ##########calculate the network properties of this network in Pn_prop.csv file
				##########write both the edges in single file Pn.sif 
				graph = set.edge.attribute(graph,"P",index=E(graph), E(graphP)$weight);
				out = paste(args[1],"html",paste(NAME,".sif.txt",sep=""),sep="/");
				write.table(cbind(get.edgelist(graph),E(graph)$weight,E(graph)$P),file = out, sep=" ", row.names=F, quote=F)
			}
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
	file.remove(paste(args[1],"html",paste(NAME,".sif",sep="")));
	file.remove(paste(args[1],"html",paste(NAME,"_prop.csv",sep="")));
	}
}
###### get RW and p.value of all TF-TF and TF-TG interaction (only those with p <= 0.05) and write it in single text file  

############ This is the most Time consuming step: but most important ##############
########################################################################

if ((!file.exists(paste(args[1],"GRN.ssv",sep="/"))) || (!file.exists(Netfile)) || (!file.exists(result2)) || (!file.exists(out_file) )) ## if one of this file do not exists the do following analysis
{
	message("Calculating P-value...");
	P  = (-1 * log10(P))
	P[P < 1.30103]<-0; ### ; 1.30103 means p value is less than 0.05 as we have log transformed the p value
	P[P==Inf]<- 6;
	message("GRN analysis...")
	SP = Matrix(P, sparse = TRUE) ## sparse matrix 
	rm(P)


	if (!file.exists(paste(args[1],"GRN.ssv",sep="/")) )
	{
		if (GPATH != "database/GRN.ssv")
		{
		TF_path = paste(args[1],"TF.txt",sep="/");
		TG_path = paste(args[1],"TG.txt",sep="/");
		TF = scan (TF_path , what="character");
		TF_C = intersect(TF,selectedGenes)
		TG = scan (TG_path, what="character");
		TG_C = intersect(TG,selectedGenes)
		}else
		{
		TF = scan ("database/TF.txt", what="character");
		TF_C = intersect(TF,selectedGenes)
		TG = scan ("database/TG.txt", what="character");
		TG_C = intersect(TG,selectedGenes)
		}
	GRN  = SP[TF_C,TG_C];
	write.table(Melt(GRN),paste(args[1],"GRN.ssv",sep="/") ,sep= " ", quote=F, row.names=F)  ## melting the sparse matrix; MELTING WITH MEFA4
	rm(GRN,TF,TF_C,TG,TG_C)
	}else
	{
	message("Skipping GRN reconstruction. File already exists");
	}
	message("Gene strength calculation...");
	ig1 = graph_from_adjacency_matrix(SP, mode=c("undirected"), weighted=TRUE, diag=FALSE)
	rm(SP);
	message("...");
	if (!file.exists(Netfile))
	{
	GS = strength(ig1, loop=F)
	write.table(GS, file=Netfile, row.names=T,quote=F,sep=" "); ##
	write.table(c(vcount(ig1), ecount(ig1)), file=strInfo, quote=F, row.names=F)
	rm(GS)
	}else{
	message("Skipping Gene strength calculation. File already exists");
	}

	if (!file.exists(result2))
	{
	message("Calculating network properties...")
	###############################################################
	message("|-----Transitivity")
	r1 = transitivity(ig1,type="local", isolates = "zero");
	#r2 = estimate_betweenness(ig1, directed = FALSE, cutoff = 3);
	message("|-----Degree")
	r2 = degree(ig1);
	message("|-----Closeness (Time consuming step)")
	r3 = estimate_closeness(ig1, cutoff = 3);
	message("|-----Page rank")
	r = page_rank(ig1, algo = "prpack" , directed = FALSE);
	r4 = r$vector

	f1 = (r1 - min(r1)) / (max(r1) - min(r1));
	f2 = (r2 - min(r2)) / (max(r2) - min(r2));
	f3 = (r3 - min(r3)) / (max(r3) - min(r3));
	f4 = (r4 - min(r4)) / (max(r4) - min(r4));
	f1[is.na(f1)] <- 0
	f2[is.na(f2)] <- 0
	f3[is.na(f3)] <- 0
	f4[is.na(f4)] <- 0
	##### combine the score ######
	Final = (f1 + f2 + f3 + f4); ## for each gene a common centrality score 
	write.table(cbind(Final,f1,f2,f3,f4), file=result2, row.names=T, quote =  FALSE,sep="\t");
	}else
	{
	message("Skipping centrality calculation. File already exists");
	}

	if (!file.exists(out_file) )
	{
		message("Predicting disease genes...")
		ig1 = delete_edge_attr(ig1,"weight"); ## if not used; the rwr was creating error ...
		s = scan(SPATH,what="characters");
		s = intersect(s,selectedGenes); ## hits with weight 1

	########################

		if (mtype == "Knet")
		{
		###### Method 1 
		########### New Method: SANTA #######
		Vertex = as.vector(V(ig1));
		other = setdiff(Vertex,s); ## non hits with 0 weight 
		names(Vertex) = Vertex; 
		Vertex[s] <- 1; Vertex[other] <- 0;
		ig1 <- set.vertex.attribute(ig1, name="hits", value=Vertex) 
		Score = Knode(ig1, vertex.attr="hits")
		}
		if (mtype == "RWR")
		{
		####### Method 2 
		############ New Method: RWR ########
		names(s) = s;
		prob = 1/length(s)
		s = as.data.frame(s)
		s[] <- prob
		RWR = dRWR(g=ig1,normalise = "laplacian",setSeeds = s, restart = 0.95, verbose = T)
		names<-V(ig1)$name;
		Score = cbind(names,(RWR[,1]))
		}
		write.table(Score, file=out_file, row.names=F,quote=F,sep=" ");
	}else
	{
	message("Skipping disease gene calculation. File already exists");
	}
}else
	{
	message("Skipping other analysis. Files already exists");
	}
save(ig1,file=igraphSession)
q("no");
