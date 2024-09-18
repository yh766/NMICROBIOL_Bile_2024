#' ---
#' title: "20220208_bilemetabolome"
#' author: "Yuko Hasegawa"
#' date: "2/8/2022"
#' output: html_document
#' editor_options: 
#'   chunk_output_type: console
#' ---
#' 
#' # packages
## -----------------------------------------------------------------------------------------------
library(tidyverse)
library(reshape2)
library(pheatmap)
library(ComplexHeatmap)
library(viridis)
library(RColorBrewer)
library(ggrepel)
library(plotly)
library(cluster)

#' 
#' # read data
## -----------------------------------------------------------------------------------------------
data <- read.csv("HARV-06-19VW_raw.csv", na.strings=c("NA","NaN", ""))

names <- c("spf1","spf2","spf3","spf4","spf5","spf6","spf7","spf8","lm1","lm2","lm3","lm4","lm5","lm6","lm7","cr1","cr2","cr3","cr4","cr5","cr6","gf1","gf2","gf3","gf4","gf5","gf6","gf7","gf8","gf9")

colnames(data)[2:31] <- names

raw <- data[,2:31]
rownames(raw) <- data[,1]

# change the order (GF group comes first)
raw <- raw[,c(22:30, 1:21)]

#' 
#' # z-score
## -----------------------------------------------------------------------------------------------
# missing values are imputed with the minimum
raw <- t(apply(raw, 1, function(x) 
                          replace(x, is.na(x), min(x[x > 0], na.rm = TRUE))))

# z-score
t(scale(t(raw))) -> zscaledata
data.frame(zscaledata) -> zscaledata

#' 
#' # PCA
## -----------------------------------------------------------------------------------------------
# transpose the matrix
pca <- prcomp(t(zscaledata))

# calculate variances
pca.var <- pca$sdev^2

# create a tibble with a variable indicating the PC number and a variable with the variances
pca.var <- tibble(PC = factor(1:length(pca.var)), variance = pca.var) %>% mutate(pct = round(variance/sum(variance)*100, 1)) %>% mutate(pct.cum = cumsum(pct))

# variance explained by PCs
ggplot(pca.var, aes(x = PC)) + geom_col(aes(y = pct)) + geom_line(aes(y = pct.cum, group = 1)) + geom_point(aes(y = pct.cum)) + labs(x = "Principal component", y = "Fraction variance explained") + theme_bw()

# the PC scores are stored in the "x" value of the prcomp object
pc.scores <- pca$x

# convert to a tibble retaining the sample names as a new column
pc.scores <- pc.scores %>% as_tibble(rownames = "sample") %>% mutate(condition = factor(gsub('[[:digit:]]+', '', sample), levels=c("gf","spf","lm","cr")))

# create a 3D PCA plot
plot_ly(pc.scores, x = ~PC1, y = ~PC2, z = ~PC3) %>% add_markers(color = ~condition, colors = c("orchid4","dodgerblue3", "darkorange2", "seagreen4"), alpha=0.7)

#' 
#' # kmeans and heatmap
## -----------------------------------------------------------------------------------------------
data <- zscaledata

# elbow plot
wss <- (nrow(data)-1)*sum(apply(data,2,var,na.rm=TRUE))
for (i in 2:20) wss[i] <- sum(kmeans(data,
                                     centers=i)$withinss)
plot(1:20, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares") #W4xH4

# heatmap
data$Pam <- factor(paste0("c", pam(data, 7)$clustering), levels=c("c3","c5","c6","c7","c2","c1","c4")) # clustering

annotation_col = data.frame(conditions = factor(c(rep("GF",9), rep("SPF",8), rep("Lm",7), rep("Cr",6)), levels=c("GF","SPF","Lm","Cr")))

rownames(annotation_col) = colnames(data[-31])

pam_col <- data.frame(clusters=data$Pam)
rownames(pam_col) <- rownames(data)

ann_colors = list(
    conditions = c(GF = "#9524ba", SPF = "#125ed8", Lm = "#f79406", Cr = "#02935f"),
    clusters = c(c1 = "#aca47c", c2 = "#171236", c3 = "#685704", c4 = "#5ca0bc", c5 = "#fcc103", c6 = "#36998d", c7 = "#a661ff")
)

pheatmap(data[order(data$Pam, decreasing = TRUE),-31], annotation_col=annotation_col, annotation_row=pam_col, cluster_rows=F, cluster_cols=F, show_rownames=F, breaks = seq(-2, 4, length.out = 50), color=colorRampPalette(c("darkslateblue","white","darkgoldenrod1"))(50), column_split=annotation_col$conditions, column_gap = unit(0.05,"in"), annotation_colors = ann_colors) #H4xW5, Landscape

# color change
pheatmap(data[order(data$Pam, decreasing = TRUE),-31], annotation_col=annotation_col, annotation_row=pam_col, cluster_rows=F, cluster_cols=F, show_rownames=F, border_color = NA, border=FALSE, breaks = seq(-2, 4, length.out = 10), color=inferno(10), column_split=annotation_col$conditions, column_gap = unit(0.01,"in"), annotation_colors = ann_colors) #H4xW5, Landscape

#' 
#' # save heatmap as a png file
## -----------------------------------------------------------------------------------------------
png("heatmap_kmeans.png", res=600, width=7.46, height=5.99, units="in")

pheatmap(data[order(data$Pam, decreasing = TRUE),-31], annotation_col=annotation_col, annotation_row=pam_col, cluster_rows=F, cluster_cols=F, show_rownames=F, border=FALSE, breaks = seq(-2, 4, length.out = 10), color=inferno(10), column_split=annotation_col$conditions, column_gap = unit(0.05,"in"), annotation_colors = ann_colors)

dev.off()

#' 
#' # Violin plot
## -----------------------------------------------------------------------------------------------
# violin plot
stat <- data %>% group_by(Pam) %>% summarize(gf=rowMeans(across(contains("gf"))), spf=rowMeans(across(contains("spf"))), lm=rowMeans(across(contains("lm"))), cr=rowMeans(across(contains("cr"))))

stat.melt <- melt(stat)

ggplot(stat.melt, aes(x=variable, y=value, fill=variable)) + geom_violin() + ylab("Average Z-score") + facet_wrap(~fct_rev(Pam),  ncol=1) + geom_hline(yintercept=0, color = "black", linewidth=1) + scale_fill_manual(values=c("#9524ba", "#125ed8", "#f79406", "#02935f"))  + theme_bw()

