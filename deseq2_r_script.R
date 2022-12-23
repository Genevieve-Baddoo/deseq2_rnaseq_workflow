#install BiocManager
#install DESeq2
#where materials are: /Users/genevievebaddoo/Downloads/salmon_gene_expression_dataset (use this directory)

#install/load necessary packages
#BiocManager::install("tximport")
library("tximport") #imports quantification data 
library("readr") #reads data in csv format
#BiocManager::install("tximportData")
library("tximportData") #provides output of running various transcript 
#abundance quantifiers on a set of 6 RNA-seq samples (6 salmon RNA-seq samples (SRRs))
library("DESeq2") #tool for differential expression analysis of RNA-seq data 

dir<-"/Users/genevievebaddoo/Downloads/salmon_gene_expression_dataset" #set directory

samples<-read.table(file.path(dir,"samples.txt"), header=TRUE) #read in samples text file
rownames(samples)<-samples$run #rownames = samples
files<-file.path(dir, samples$run, "quant.sf") #construct of gene-level DESeqDataSet object from Salmon quant.sf files, 
#which are stored in tximportData package
names(files)<-samples$run #specify path to files using appropriate columns of samples,
#read in table that links transcripts to genes for this dataset
tx2gene<-read_csv(file.path(dir, "TEs_tx2gene.csv")) #read csv file
txi<-tximport(files, type="salmon", tx2gene=tx2gene) #import necessary quantification data for
#DESeq2 using the tximport function
ddsTxi<-DESeqDataSetFromTximport(txi, colData=samples, design= ~ condition) #construct DESeqDataSet 
#(dds) from the txi object and sample information 

#pre-filtering: pre filters low count genes
keep <- rowSums(counts(ddsTxi)) >= 10 #keep rows that have at least 10 reads total
dds <- ddsTxi[keep,] #dds = deseq data set. Object used to store read counts 

#relevel reference levels
dds$condition <- relevel(dds$condition, ref = "control")
dds$condition <- droplevels(dds$condition) #remove levels that don't have samples in current dds
dds <- DESeq(dds) #standard differential expression analysis steps
res <- results(dds) #generate results (plot this)
#res
#head(results(dds, tidy=TRUE)) view results table

#log2 fold change: condition mutant vs control

#reset par
par(mfrow=c(1,1))

# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="spnE", xlim=c(-4,6)))

# Add colored points: red if padj<0.05
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

