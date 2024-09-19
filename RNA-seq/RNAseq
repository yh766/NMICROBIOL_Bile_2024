#' ---
#' title: "RNAseq_DESeq2_clusterProfiler"
#' author: "Yuko Hasegawa"
#' date: "9/21/2022"
#' output: html_document
#' editor_options: 
#'   chunk_output_type: console
#' ---
#' 
#' # packages
## ---------------------------------------------------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(limma)
library(stringr)
library(rtracklayer)
library(cluster)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library(pals)
library(reshape2)
library(clusterProfiler)
library(msigdbr)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(enrichplot)

#' 
#' # run DESeq2
## ----setup, include=FALSE---------------------------------------------------------------------------------
# read count and create count data
rawData.lm <- read.table("Lm_count.txt", header=T, sep="\t")
rawData.cr <- read.table("Cr_count.txt", header=T, sep="\t")

countData.lm <- rawData.lm[,-c(2:6)]
countData.cr <- rawData.cr[,-c(2:6)]

names.lm <- c("Geneid","SPF1","SPF2","SPF3","Lm1","Lm2","Lm3")
names.cr <- c("Geneid","SPF4","SPF5","SPF6","Cr.n1","Cr.p1","Cr.n2","Cr.n3","Cr.p2","Cr.p3")

colnames(countData.lm) <- names.lm
colnames(countData.cr) <- names.cr

# reorder CR columns
countData.cr <- countData.cr[, c(1:4,5,7:8,6,9:10)]

# combine two dataset
countData <- left_join(countData.lm, countData.cr, by="Geneid")
countData <- countData[,c(1:4,8:10,5:7,11:16)]

# create meta data
metaData <- data.frame(id=colnames(countData)[-1], condition=c(rep("SPF", 6), rep("Lm", 3),  rep("Cr.neg", 3), rep("Cr.pos", 3)), tissue=rep("Liver",15))

# construct DESEQDataset object
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~condition, tidy = TRUE)#tidy; for matrix input: whether the first column of countData is the rownames for the count matrix

# run DESeq
dds <- DESeq(dds)

#' 
#' # lfcShrink
## ---------------------------------------------------------------------------------------------------------
dds$condition <- relevel(dds$condition, ref = "SPF") # set reference sample
dds <- nbinomWaldTest(dds)
resultsNames(dds)

lm.srnk <- lfcShrink(dds, coef="condition_Lm_vs_SPF", type="apeglm")
cr.neg.srnk <- lfcShrink(dds, coef="condition_Cr.neg_vs_SPF", type="apeglm")
cr.pos.srnk <- lfcShrink(dds, coef="condition_Cr.pos_vs_SPF", type="apeglm")

lm.srnk.df <- as.data.frame(lm.srnk)
cr.neg.srnk.df <- as.data.frame(cr.neg.srnk)
cr.pos.srnk.df <- as.data.frame(cr.pos.srnk)

#' 
#' # convert gene ID
## ----echo=FALSE-------------------------------------------------------------------------------------------
gtf <- import("gencode.vM25.annotation.gtf")

gtf.df <- as.data.frame(gtf)
geneID <- subset(gtf.df, type=="gene")[,c("gene_id","gene_type","gene_name")]

list.srnk <- list(lm.srnk.df, cr.neg.srnk.df, cr.pos.srnk.df)
names(list.srnk) <- c("lm.srnk","cr.pos.srnk","cr.neg.srnk")

for (i in 1:length(list.srnk)){
  list.srnk[[i]] <- rownames_to_column(list.srnk[[i]], var="geneid")
  list.srnk[[i]] <- left_join(list.srnk[[i]], geneID, by=c("geneid"="gene_id"))
  list.srnk[[i]] <- list.srnk[[i]][,c(1,7,8,2,3,5,6)]
  colnames(list.srnk[[i]])[-c(1:3)] <- c(paste(colnames(list.srnk[[i]])[-c(1:3)], names(list.srnk[i]), sep="."))
}

#' 
#' # combine data
## ---------------------------------------------------------------------------------------------------------
table.srnk <- list.srnk %>% purrr::reduce(left_join, by = c("geneid","gene_type","gene_name"))

table.srnk <- table.srnk %>% mutate(class=case_when(padj.lm.srnk < 0.05 & padj.cr.pos.srnk < 0.05 ~ "Both", padj.lm.srnk < 0.05 & padj.cr.pos.srnk > 0.05 ~ "Lm DEGs", padj.lm.srnk > 0.05 & padj.cr.pos.srnk < 0.05 ~ "Cr DEGs", TRUE ~ "not DEGs")) # assign class

#' 
#' # pca
## ---------------------------------------------------------------------------------------------------------
vsdata <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsdata, intgroup="condition", returnData=TRUE)
pcaData$condition <- factor(pcaData$condition, levels=c("SPF","Lm", "Cr.neg", "Cr.pos"))
percentVar <- round(100 * attr(pcaData, "percentVar"))
p <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values=c(SPF = "#125ed8", Lm = "#f79406", Cr.pos = "#02935f", Cr.neg = "#99CC99")) + theme_bw() # W4 x H5

p + ggrepel::geom_text_repel(data=pcaData, aes(PC1, PC2, label=name), size = 3, show.legend = FALSE)

#' 
#' # z-score table
## ---------------------------------------------------------------------------------------------------------
# retrieve the normalized counts
vsd <- vst(dds, blind=FALSE)
assay(vsd) -> counts.nom # extract normalized value
data.frame(counts.nom) -> counts.nom

# row Z-score
t(scale(t(counts.nom))) -> counts.nom.z
data.frame(counts.nom.z) -> counts.nom.z

#' 
#' # clustering
## ---------------------------------------------------------------------------------------------------------
# extract DEGs
table.sub <- subset(table.srnk, padj.lm.srnk < 0.05 | padj.cr.pos.srnk < 0.05)

table.sub <- subset(table.sub, abs(log2FoldChange.lm.srnk) > 1 | abs(log2FoldChange.cr.pos.srnk) > 1)
genes <- table.sub$geneid

data <- subset(counts.nom.z, rownames(counts.nom.z) %in% genes)

# elbow plot (optimize the number for clustering)
wss <- (nrow(data)-1)*sum(apply(data,2,var,na.rm=TRUE))
for (i in 2:20) wss[i] <- sum(kmeans(data,
                                     centers=i)$withinss)
plot(1:20, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares") #W4xH4

# briefly check (! cluster numbers are not same to the result of the next clustering)
n = 7 # how many clusters do you want
pheatmap(data, cluster_cols=F, kmeans_k=n)

# assign cluster to genes
data$Pam <- factor(paste0("c", pam(data, 7)$clustering)) # adjust levels after checking the first result

#' 
#' # heatmap
## ---------------------------------------------------------------------------------------------------------
# create annotation
annotation_col = data.frame(conditions = factor(c(rep("SPF",6), rep("Lm",3), rep("Cr.neg",3), rep("Cr.pos",3)), levels=c("SPF","Lm","Cr.pos","Cr.neg")))

rownames(annotation_col) = colnames(counts.nom.z)

pam_col <- data.frame(clusters=data$Pam)
rownames(pam_col) <- rownames(data)

# colors for conditions and clusters
ann_colors = list(
    conditions = c(SPF = "#125ed8", Lm = "#f79406", Cr.pos = "#02935f", Cr.neg = "#99CC99"),
    clusters = c(c1 = "#aca47c", c2 = "#171236", c3 = "#685704", c4 = "#5ca0bc", c5 = "#fcc103", c6 = "#36998d", c7 = "#a661ff")
)

# data range
range <- abs(max(data[order(data$Pam, decreasing = TRUE), !(names(data) %in% "Pam")]))

# reorder clusters
data$Pam <- factor(paste0("c", pam(data, 7)$clustering), levels=c("c1","c5","c4","c3","c7","c6","c2")) # adjust the order of the clusters

# heatmap (1st check, please see the order of the clusters)
pheatmap(data[order(data$Pam, decreasing = TRUE), !(names(data) %in% "Pam")], annotation_col=annotation_col, annotation_row=pam_col, cluster_rows=F, cluster_cols=F, show_rownames=F, breaks = seq(-range, range, length.out = 50), color=colorRampPalette(c("darkslateblue","white","darkgoldenrod1"))(50), column_split=annotation_col$conditions, column_gap = unit(0.05,"in"), annotation_colors = ann_colors, use_raster = F) #H4xW5, Landscape

# change color
pheatmap(data[order(data$Pam, decreasing = TRUE), !(names(data) %in% "Pam")], annotation_col=annotation_col, annotation_row=pam_col, cluster_rows=F, cluster_cols=F, show_rownames=F, breaks = seq(-range, range, length.out = 50), color=rev(brewer.rdbu(50)), column_split=annotation_col$conditions, column_gap = unit(0.05,"in"), annotation_colors = ann_colors, use_raster = F) #H4xW5, Landscape

#' 
#' # Violin plot
## ---------------------------------------------------------------------------------------------------------
# violin plot
stat <- data %>% group_by(Pam) %>% summarize(SPF=rowMeans(across(contains("SPF"))), Lm=rowMeans(across(contains("Lm"))), Cr.pos=rowMeans(across(contains("Cr.p"))), Cr.neg=rowMeans(across(contains("Cr.n"))))

stat.melt <- melt(stat)

ggplot(stat.melt, aes(x=variable, y=value, fill=variable)) + geom_violin() + ylab("Average Z-score") + facet_wrap(~fct_rev(Pam),  ncol=1) + geom_hline(yintercept=0, color = "black", size=1) + scale_fill_manual(values=c(SPF = "#125ed8", Lm = "#f79406", Cr.pos = "#02935f", Cr.neg = "#99CC99"))  + theme_bw()

#' 
#' ## pathway analysis
#' # compareCluster KEGG (clusters)
## ---------------------------------------------------------------------------------------------------------
set.seed(1234)

# create genelist
clusters <- rev(levels(data$Pam))
genelist <- list()

for(i in 1:length(clusters)){
  table <- subset(data, Pam==clusters[i])
  table$entrez = mapIds(org.Mm.eg.db, keys=gsub("\\.[0-9]*","",table$geneid), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
  genelist[[i]] <- na.omit(table$entrez)
}

names(genelist) <- clusters

# run compareCluster
cc.kegg <- compareCluster(geneCluster = genelist, fun = "enrichKEGG", organism="mmu", keyType = "ncbi-geneid", pAdjustMethod = "fdr", qvalueCutoff=0.2, minGSSize = 10, maxGSSize = 500)

cc.kegg <- setReadable(cc.kegg, OrgDb = org.Mm.eg.db, keyType="ENTREZID")

cc.kegg.result <- cc.kegg@compareClusterResult

#' 
#' # dotplot
## ---------------------------------------------------------------------------------------------------------
# KEGG metabolism pathways
dotplot(cc.kegg, showCategory=100) + scale_color_gradientn(colours=rev(brewer.blues(100)), limits=c(0,0.1)) + theme(axis.title.x=element_blank())

pathways <- cc.kegg.result$Description[grepl(paste("^mmu00", "^mmu01", sep="|"), cc.kegg.result$ID)]

dotplot(cc.kegg, showCategory=pathways) + scale_color_gradientn(colours=rev(brewer.blues(100)), limits=c(0,0.1)) + theme(axis.title.x=element_blank()) # W8 x H10

#' 
#' # compareCluster
#' # read DESeq2 results
## ---------------------------------------------------------------------------------------------------------
# create list object
table.lm <- read.delim("DESeq2_Lm.srnk.txt", comment.char="#")

table.cr.pos <- read.delim("DESeq2_Cr_positive.srnk.txt", comment.char="#")

table.cr.neg <- read.delim("DESeq2_Cr_negative.srnk.txt", comment.char="#")

data <- list(table.lm, table.cr.pos, table.cr.neg)
genelist <- list()

for(i in 1:length(data)){
  table <- data[[i]]
  table$entrez = mapIds(org.Mm.eg.db, keys=gsub("\\.[0-9]*","",table$geneid), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
 genelist[[i]] <- na.omit(sort(setNames(table$log2FoldChange, table$entrez), decreasing = T))
}

names(genelist) <- c("lm", "cr.pos", "cr.neg")

#' 
#' # compareCluster KEGG
## ---------------------------------------------------------------------------------------------------------
set.seed(1234)

cc.kegg <- compareCluster(geneCluster = genelist, fun = "gseKEGG", organism="mmu", pAdjustMethod = "fdr", eps = 0, pvalueCutoff=0.1, minGSSize = 10, maxGSSize = 500, seed =TRUE)

cc.kegg <- setReadable(cc.kegg, OrgDb = org.Mm.eg.db, keyType="ENTREZID")

cc.kegg.result <- cc.kegg@compareClusterResult

#' 
#' # compareCluster MSigDb
## ---------------------------------------------------------------------------------------------------------
set.seed(1234)

msig <- msigdbr(species = "Mus musculus", category = "H") %>% dplyr::select(gs_name, entrez_gene) %>% mutate(gs_name=gsub("HALLMARK_","",gs_name))

cc.msig <- compareCluster(geneCluster = genelist, fun = "GSEA", pAdjustMethod = "fdr", eps = 0, pvalueCutoff=0.2, minGSSize = 10, maxGSSize = 500, seed =TRUE, TERM2GENE=msig)

cc.msig <- setReadable(cc.msig, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
cc.msig.result <- cc.msig@compareClusterResult

#' 
#' # dotplot
## ---------------------------------------------------------------------------------------------------------
# MSigDB dotplot
pathways = unique(subset(cc.msig.result, p.adjust < 0.05)$Description)

dotplot(cc.msig, showCategory=pathways, size="setSize") + scale_color_gradientn(colours=rev(brewer.blues(100)), limits=c(0,0.2)) + theme(axis.title.x=element_blank())

dp <- dotplot(cc.msig, showCategory=pathways, size="setSize") + scale_color_gradientn(colours=rev(brewer.blues(100)), limits=c(0,0.2))

ggplot(dp$data, aes(x=Cluster, y=Description, size=setSize, color=p.adjust, fill=p.adjust, shape=.sign)) + geom_point() + theme_bw() + scale_color_gradientn(colours=rev(brewer.blues(100)), limits=c(0,0.2)) + scale_fill_gradientn(colours=rev(brewer.blues(100)), limits=c(0,0.2)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + scale_shape_manual(values=c(21, 25)) # H10 x W8

# KEGG dotplot
pathways = unique(subset(cc.kegg.result, p.adjust < 0.05)$Description)

dotplot(cc.kegg, showCategory=pathways, size="setSize") + scale_color_gradientn(colours=rev(brewer.blues(100)), limits=c(0,0.2)) + theme(axis.title.x=element_blank())

dp <- dotplot(cc.kegg, showCategory=pathways, size="setSize") + scale_color_gradientn(colours=rev(brewer.blues(100)), limits=c(0,0.2))

ggplot(dp$data, aes(x=Cluster, y=Description, size=setSize, color=p.adjust, fill=p.adjust, shape=.sign)) + geom_point() + theme_bw() + scale_color_gradientn(colours=rev(brewer.blues(100)), limits=c(0,0.2)) + scale_fill_gradientn(colours=rev(brewer.blues(100)), limits=c(0,0.2)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + scale_shape_manual(values=c(21, 25)) # H10 x W8

