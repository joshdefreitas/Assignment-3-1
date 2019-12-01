# open a new R script
rm(list = ls())
#  1- for each question insert section using "ctrl+shift+R" and label it (i.e, Q1)
#  2- write comment (but not long comments)


# load required libraries -------------------------------------------------
#  Read sample gene counts with a csv file.
#install.packages("BiocManager") 
#install.packages("pheatmap") 
#install.packages("tidyverse")

library("tidyverse")
library( "DESeq2" )
library("pheatmap")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("org.Mm.eg.db") 
library("dplyr")
library(pathview)
library(gage)
library(gageData)

# Load data ---------------------------------------------------------------
#  code here....
countdata <- read.csv("countdata.csv")
metadata <- read.csv("metaData.csv")

countdata <- countdata %>% remove_rownames %>% column_to_rownames(var="X")
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var="X")

# Q1 ----------------------------------------------------------------------

#Counting genes that have atleast one sample with 0 expression
count = 0
for(j in 1:nrow(countdata))
{
  for(k in 1:ncol(countdata))
  {
    if(countdata[j,k]==0){
      count = count + 1
      break()
    }
  }
}

print(count)

# Q2 ----------------------------------------------------------------------

metadata$experiment.number <- factor(metadata$experiment.number)
metadata$experiment.number <- factor(metadata$num.tech.reps)

metadata$lane.number <- factor(metadata$lane.number)

#Running DESeq2
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = metadata, design = ~ strain, tidy = FALSE, ignoreRank = FALSE )
dds <- DESeq(dds)
res <- results(dds)

sum(res$padj < 0.1,na.rm = TRUE)


# Q3 ----------------------------------------------------------------------

#Changing alpha
alpha <- 0.07
upDownFDRAlpha = results(dds, alpha=alpha)
summary(upDownFDRAlpha)

# Q4 ----------------------------------------------------------------------
resTop<- upDownFDRAlpha[upDownFDRAlpha$log2FoldChange > 0,]
resTop <- resTop[order(resTop$padj),]
resDown<- upDownFDRAlpha[upDownFDRAlpha$log2FoldChange < 0,]
resDown<- resDown[order(resDown$padj),]

print(resTop[1:20,])
print(resDown[1:20,])


# Q5 ----------------------------------------------------------------------
plotMA(res, ylim=c(-2,2))

# Q6 ----------------------------------------------------------------------

plotCounts(dds, gene=which.min(res$pvalue), intgroup = "strain")
meanp <- mean(res$pvalue,na.rm=TRUE)
plotCounts(dds,gene=which.min(meanp), intgroup = "strain")

# Q7 ----------------------------------------------------------------------
gene <- rownames(res[which.max(res$pvalue),])[1]
print(gene)
high_pval_count = counts(dds[gene,],)


# Q8 ----------------------------------------------------------------------
library("AnnotationDbi")
library("org.Hs.eg.db")

print(row.names(res))
library("org.Mm.eg.db")
res$symbol <- mapIds(org.Mm.eg.db,
                     keys=rownames(res), 
                     column="SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Mm.eg.db,
                     keys=row.names(res), 
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
res$name <- mapIds(org.Mm.eg.db,
                   keys=row.names(res), 
                   column="GENENAME",
                   keytype="ENSEMBL",
                   multiVals="first")
write.csv(res,'resultPath.csv')


# Q9 ----------------------------------------------------------------------


data(kegg.sets.mm)
data(sigmet.idx.mm)
kegg.sigmet.mm=kegg.sets.mm[sigmet.idx.mm]
head(kegg.sets.mm, 3)

foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
head(foldchanges)

keggres = gage(foldchanges, gsets=kegg.sets.mm, same.dir=TRUE)

lapply(keggres, head)

summary(keggres)

keggrespathways = data.frame(id=rownames(keggres$less), keggres$less) %>% 
  tbl_df() %>% 
  filter(row_number()<=6) %>% 
  .$id %>% 
  as.character()

keggresids = substr(keggrespathways, start=1, stop=8)

print(keggresids)


# Q10 ---------------------------------------------------------------------
keggrespathwaysup <- data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
  tbl_df() %>% 
  filter(row_number()<=15) %>% 
  .$id %>% 
  as.character()

keggrespathwaysup <- order(keggrespathwaysup)
normt <-normTransform(dds)
df <- as.data.frame(colData(dds)['strain'])
pheatmap(assay(normt)[keggrespathwaysup,], show_rownames=TRUE, annotation_col=df)


