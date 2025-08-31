#Load the required libraries

library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(rstudioapi)     

# This is to find the correct pathway, only rstudio users can use
setwd(dirname(getActiveDocumentContext()$path))  

#########################################################################################
#Based on FGT_T7_Advanced_Analysis_RNASeq_Using R.pdf
#URL https://www.learn.ed.ac.uk/ultra/courses/_118861_1/outline/file/_11082900_1
#########################################################################################
#Load data (I already have them, so directly load instead of create them like tutorial)
load("mytable_feaures")
adf <-read.csv("SraRunTable.csv",sep=',',row.names=1,fill=T,header=T)

#get table_rnaseq, the same as tutorial
table_rnaseq<-(mytable_feaures)$counts
#tide up sample names, tutorial don't show code for this step
colnames(table_rnaseq) <- gsub("Aligned.sortedByCoord.out.bam$", "", colnames(table_rnaseq))
colnames_rnaseq <-colnames(table_rnaseq)
#The same method as in tutorial(I use 'Sample_short', not 'ShortName')
idx<-match(colnames_rnaseq,rownames(adf))
colnames(table_rnaseq) <-adf$Sample_short[idx]
adf <-adf[idx,]
rownames(adf) <-colnames(table_rnaseq)
head(table_rnaseq)

#I don't have to remove things, so skip "Removal of Unwanted Samples" part

#"Loading the Data into DESeq" 
#I can use everything for dds because I didn't remove 
# the design is only the 'type' column, means gene knocked or not
dds_rnaseq <- DESeqDataSetFromMatrix(countData = table_rnaseq,colData = adf,design = ~ type)

#"Quality Filtering Genes"
#Completely the same method as tutorial
keep <- rowSums(counts(dds_rnaseq)) > 1
dds_rnaseq <- dds_rnaseq[keep,]
nrow(dds_rnaseq)
keep <- rowSums(counts(dds_rnaseq) >= 10) >= 3
dds_rnaseq <- dds_rnaseq[keep,]
nrow(dds_rnaseq)

#"Normalisation of Data"
#Completely the same method as tutorial,I didn't add "non-normalised" steps here.
dds_rnaseq <- estimateSizeFactors(dds_rnaseq)
sizeFactors(dds_rnaseq)
counts_rnaseq <-log2(counts(dds_rnaseq, normalized=TRUE))
fpm_rnaseq <-log(fpm(dds_rnaseq))
vsd_rnaseq <- vst(dds_rnaseq, blind = T)
rld_rnaseq <- rlog(dds_rnaseq, blind = T)
boxplot(counts_rnaseq)
boxplot(fpm_rnaseq)
boxplot(assay(vsd_rnaseq))
boxplot(assay(rld_rnaseq))
#I decided to use rlog normalisation,the same as tutorial did.
#Because rlog shows the most consistent position of the box.
#And the almost completely identical median line.

#"Exploratory Data Plots", PCA
#Completely the same method as tutorial.  
library(scatterplot3d)
library(ggplot2)
pca <- prcomp(t(na.omit(assay(rld_rnaseq))), scale=T)
s3d<-scatterplot3d(pca$x[,1:3], pch=19,color=colData(dds_rnaseq)$Colour)
s3d.coords <- s3d$xyz.convert(pca$x[,1:3])
text(s3d.coords$x, s3d.coords$y, labels = colnames(rld_rnaseq),pos =
       3,offset = 0.5,cex=0.5)
qplot(pca$x[,1],pca$x[,2],xlab="PCA1",ylab="PCA2",color=colData(dds_rnaseq)$Colour)
sampleDists <- dist(t(assay(rld_rnaseq)))
sampleDists
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix)
plotPCA(rld_rnaseq, intgroup = c("type")) #just change '"Sex", "Line"' to "type"

#"Differential Gene Expression"
#Same steps as tutorial
dds_rnaseq <- DESeq(dds_rnaseq)
result_type <- results(dds_rnaseq, contrast = c("type", "cKO", "Control")) #using my design
plotMA(result_type, main = "cKO vs Control", ylim = c(-2, 2))
table(result_type$padj < 0.05)
result_type_selected <- subset(result_type, padj < 0.05)
result_type_selected <- result_type_selected[order(abs(result_type_selected$log2FoldChange), decreasing = TRUE), ]
head(result_type_selected,10)

#"make a heatmap for the top 50"
#nearly the same as tutorial, add title and annotation
top50 <- rownames(result_type_selected)[1:50]
pheatmap(
  assay(rld_rnaseq)[top50, ], 
  scale = "row",          
  show_rownames = FALSE, 
  main = "Top 50 DEGs: cKO vs Control",   # add title
  annotation_col = data.frame(            # add annotation
    Type = colData(rld_rnaseq)$type,
    row.names = colnames(rld_rnaseq)
  )
) 

pheatmap(
    assay(rld_rnaseq)[rownames(result_type_selected), ],
    scale = "row",
    show_rownames = FALSE,
    main = "All Significant DEGs: cKO vs Control",   # add title
    annotation_col = data.frame(                     # add annotation
      Type = colData(rld_rnaseq)$type,
      row.names = colnames(rld_rnaseq)
    )
)

#"Download Ensembl annotation using BiomaRt and rename the samples" 
#########################################################################################
#Based on example_analysis_GSE243520.pdf
#Founded in example_analysis_extraICA_help.zip
#########################################################################################
#the same as "Build Gene Annotation. A serialised is generated in case Ensembl is unavailable"
library(biomaRt) 
if(file.exists("resultAnnot.RData")){ 
  
  load("resultAnnot.RData") 
  print("Loaded existing Ensembl annotation...") 
}else{ 
  
  #UK ensembl is being updated we use a USA mirror, "useast.ensembl.org" 
  ensembl_host <-"https://www.ensembl.org" 
  head(biomaRt::listMarts(host = ensembl_host), 15) 
  head(biomaRt::listAttributes(biomaRt::useDataset(dataset = "mmusculus_gene_ensembl",mart = useMart("ENSEMBL_MART_ENSEMBL",host = ensembl_host))), 40) 
  mart <- biomaRt::useDataset(dataset = "mmusculus_gene_ensembl",mart = useMart("ENSEMBL_MART_ENSEMBL",host = ensembl_host)) 
  resultAnnot <- biomaRt::getBM(values=rownames(dds_rnaseq),attributes = 
                                  c("ensembl_gene_id","external_gene_name","chromosome_name","start_position",
"end_position","description","strand"),filters="ensembl_gene_id",mart=mart) 
  save(resultAnnot,file="resultAnnot.RData") 
} 

#"merge with input data"
#The same as what in FGT_T7_Advanced_Analysis_RNASeq_Using R.pdf
names <- resultAnnot[, 1]  
resultAnnot <- as.data.frame(resultAnnot)
rownames(resultAnnot) = names
idx<-match(rownames(dds_rnaseq),rownames(resultAnnot))
all(rownames(dds_rnaseq) == rownames(resultAnnot))
grr<-resultAnnot[match(rownames(dds_rnaseq),
                       resultAnnot$ensembl_gene_id),]
all(rownames(dds_rnaseq) == rownames(grr))
resultAnnot <-grr
all(rownames(dds_rnaseq) == rownames(resultAnnot))

#make the nice names
nice_names<-
  paste(resultAnnot$ensembl_gene_id,resultAnnot$external_gene_name,
        sep = '_')
resultAnnot$nice_names <-nice_names
head(resultAnnot)
all(rownames(dds_rnaseq) == rownames(resultAnnot))
#check names
rld_rnaseq <- rlog(dds_rnaseq, blind = TRUE)#
idx2 <- match(rownames(result_type_selected)[1:50], rownames(dds_rnaseq))
plotme <-(rld_rnaseq)[rownames(result_type_selected)[1:50],]
rownames(plotme)<-resultAnnot$nice_names[idx2]

#make heatmap with candidate genes
png("heatmap_candidates.png", width = 800, height = 1000)

pheatmap(assay(plotme),scale="row",fontsize_row = 10,cellheight =12,
         cellwidth =12,treeheight_row = 40, treeheight_col = 40)
dev.off()

#Final Results Table
#Write some code to achieve similar results in example_analysis_GSE243520.pdf

#Convert results to dataframe
res_df <- as.data.frame(result_type_selected)  
res_df$ensembl_gene_id <- rownames(res_df)  # add gene_id column

# merge result and annotation
annotated_results <- merge(
  res_df,
  resultAnnot[, c("ensembl_gene_id", "external_gene_name")],  
  by = "ensembl_gene_id",
  all.x = TRUE  # keep all genes, even no annotation
)

# Sort as previously did
annotated_results <- annotated_results[order(abs(annotated_results$log2FoldChange), decreasing = TRUE), ]

# Check results(only necessary columns)
head(annotated_results[, c("ensembl_gene_id", "external_gene_name", "log2FoldChange", "padj")])

#GO enrichment analysis
library(clusterProfiler)
library(org.Mm.eg.db)  
library(enrichplot)

# get the Ensembl ID of selected genes(padj < 0.05).
sig_genes <- rownames(result_type_selected)

# process analysis
ego <- enrichGO(gene = sig_genes,
                OrgDb = org.Mm.eg.db,
                keyType = "ENSEMBL",
                ont = "BP",              # BP: Biological Process
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2,
                readable = TRUE)         

# plot
barplot(ego, showCategory = 10) + ggtitle("GO Enrichment (BP)")

#KEGG enrichment analysis
#Convert ENSEMBL to ENTREZ ID（KEGG using Entrez ID）
kegg_df <- bitr(sig_genes,
                fromType = "ENSEMBL",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)

# run KEGG enrichment analysis
kegg_result <- enrichKEGG(gene = kegg_df$ENTREZID,
                          organism = 'mmu',         
                          pvalueCutoff = 0.05)

# plot
barplot(kegg_result, showCategory = 10, title = "KEGG Pathway Enrichment")

