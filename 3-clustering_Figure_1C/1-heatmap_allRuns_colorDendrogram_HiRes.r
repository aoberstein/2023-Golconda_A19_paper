#!/usr/local/bin/Rscript
library(GSA)
library(ggplot2)
library("viridis", quietly = TRUE)
library("ggplot2", quietly = TRUE)
library("dplyr", quietly = TRUE)
library("reshape2", quietly = TRUE)
library(gplots)
library("dendextend")
options("scipen"=100, "digits"=4)

# User defined variables --------------------------------------------------
ttgFile = "../ttg.tab"
metadataFile = "../metadata_merged.tab"
countFile = paste0("../2-limmaGene_GeTMM/" ,"GeTMM_allData_TMM_expValues_cpm.tab")
ttg = read.csv(ttgFile, header = T, sep = '\t', stringsAsFactors = F)

# # import counts and metadata
counts = read.csv(countFile, sep='\t', header= TRUE, stringsAsFactors = FALSE)
metadata = read.csv(metadataFile, sep='\t', header= TRUE, stringsAsFactors = FALSE)

# convert ensembl ids to hgnc gene names
gids = ifelse( ttg$ext_gene[match(counts$gid, ttg$ens_gene)] == "-", counts$gid,  
               ttg$ext_gene[match(counts$gid, ttg$ens_gene)])

run2 = gsub("__.*", "", metadata$run)
treat = paste( run2, metadata$treament,metadata$cell, sep = '|')
lookup = data.frame(metadata$run, treat)
names(lookup)[1] = "run"
counts2=counts
colnames(counts2) <- plyr::mapvalues(colnames(counts),
                                      as.character(lookup$run), 
                                      as.character(lookup$treat))

## change first column GeneIDs to rownames
matInput = data.frame(row.names = counts2[,1], counts2[,2:length(colnames(counts2))])
mat = data.matrix(matInput, rownames.force = TRUE)

# ## filter GFP, puro, and neo genes
gidsToRemove = c("EGFP_A206K", "neomycin", "puromycin")
mat2 = subset(mat, !(row.names(mat) %in% gidsToRemove)) 


## filter rows with low counts
# less than number of samples (columns) times 2
threshold = length(colnames(counts))*1
matFilt = mat2[rowSums(mat2) > threshold,]

##custom colors
col2=HexLogClock=c("#addd8e") #light green
col3=HexLogClock=c("#31a354") #green 
col4=HexLogClock=c("#fee8c8") #pink
col5=HexLogClock=c("#fdbb84") #salmon
col6=HexLogClock=c("#e34a33") #red
col_orange=HexLogClock=c("#FFA800")
col_blue=HexLogClock=c("#2C54E8")

my_palette <- colorRampPalette(c(col_blue, "white", col_orange))(n = 599)
colors = c(seq(-1, -0.11, length = 200), seq(-0.1, 0.1, length = 200), seq(0.11, 1, length = 200))

## run hclust outside of heatmap.2, to allow access to the data
matScaled = t(scale(t(matFilt), center = TRUE, scale = TRUE))
matFilt_log = log2(matFilt+1)
matScaled_log2 = t(scale(t(matFilt_log), center = TRUE, scale = TRUE))
# matHClust = hclust(as.dist(1-cor(t(matScaled), method = "spearman")), method='ward.D2')
# hclust = function(x) hclust(x, method='ward.D2'), distfun=function(x) as.dist(1-cor(t(x), method="spearman"))
png(file=paste0("1-public_EMT_data_all_samples_HEATMAP_THRESH1.png"), height = 8, width = 18, units = "in", res = 900)
hm = heatmap.2(matScaled,
          breaks=colors,
          col=my_palette,
          dendrogram = "both",
          scale = "row",
          trace = "none",
          labRow = c(""),
          density.info="none",
          key.title="",
          margins=c(12,5),
          lwid=c(2,5), #make column of dendrogram and key very small and other colum very big
          # lhei=c(3,5), #make row of key and other dendrogram very small and other row big.
          hclust = function(x) hclust(x, method='ward.D2'), distfun=function(x) as.dist(1-cor(t(x), method="spearman"))
          )
dev.off()


png(file=paste0("1-public_EMT_data_all_samples_HEATMAP__THRESH1_LOG2.png"), height = 8, width = 18, units = "in", res = 900)
hm_log = heatmap.2(matScaled_log2,
               breaks=colors,
               col=my_palette,
               dendrogram = "both",
               scale = "row",
               trace = "none",
               # labRow = c(""),
               density.info="none",
               key.title="",
               margins=c(12,5),
               lwid=c(2,5), #make column of dendrogram and key very small and other colum very big
               # lhei=c(3,5), #make row of key and other dendrogram very small and other row big.
               hclust = function(x) hclust(x, method='ward.D2'), distfun=function(x) as.dist(1-cor(t(x), method="spearman"))
)
dev.off()

colColors1 = gsub(".*RWPE.1.*", col_orange, colnames(matScaled_log2))
colColors2 = gsub(".*MCF10A.*", col6, colColors1)
colColors3 = gsub(".*ARPE.19.*", "black", colColors2)
colColors4 = gsub(".*MRC5.*", col_blue, colColors3)

# cols_branches = c(HexLogClock=c("#cccccc"),
#                   HexLogClock=c("#969696"),
#                   HexLogClock=c("#525252"),
#                   "black"
# 
dendColors = c(col6, col_orange, col_blue, "black" )
dend1 = color_branches(hm_log$colDendrogram, k = 4, col = dendColors)
dend1 = set(dend1, "branches_lwd", 1)
dend2 = set(hm_log$rowDendrogram, "branches_lwd", 0.6)
col_labels = get_leaves_branches_col(dend1)
col_labels = col_labels[order(order.dendrogram(dend1))]

# these are the cell-type colors used in the rest of the paper (e.g. jitter plots)
#colors = c("ARPE-19" = "black", "MCF10A" = col6 , "MRC5" = col_blue  , "RWPE-1" = col_orange)

png(file=paste0("1-public_EMT_data_all_samples_HEATMAP__THRESH1_LOG2_COLORED-3.png"), 
    height = 3, width = 4.5, units = "in", res = 900)
    heatmap.2(matScaled_log2,
                   breaks=colors,
                   col=my_palette,
                   dendrogram = "both",
                   Colv = dend1,
                   Rowv = dend2,
                   scale = "row",
                   labRow = c(""),
                   labCol = c(""),
                   # colCol = col_labels,
                   # colCol = colColors4,
                   ColSideColors = colColors4, 
                   trace = "none",
                   density.info="none",
                   key.title="",
                   # margins=c(12,5),
                   lwid=c(1,5), #make column of dendrogram and key very small and other colum very big
                   lhei=c(1,4), #make row of key and other dendrogram very small and other row big.
                   hclust = function(x) hclust(x, method='ward.D2'), 
                   distfun=function(x) as.dist(1-cor(t(x), method="spearman"))
)
dev.off()

#############################################33333 black dendroggram, colored "sidebar"
dendColors = c("black", "black", "black", "black" )
dend1 = color_branches(hm_log$colDendrogram, k = 4, col = dendColors)
dend1 = set(dend1, "branches_lwd", 1)
dend2 = set(hm_log$rowDendrogram, "branches_lwd", 0.6)
col_labels = get_leaves_branches_col(dend1)
col_labels = col_labels[order(order.dendrogram(dend1))]

# these are the cell-type colors used in the rest of the paper (e.g. jitter plots)
#colors = c("ARPE-19" = "black", "MCF10A" = col6 , "MRC5" = col_blue  , "RWPE-1" = col_orange)

png(file=paste0("1-public_EMT_data_all_samples_HEATMAP_THRESH1_LOG2_COLORED_BLACK_DEND-4.png"), 
    height = 3, width = 4.5, units = "in", res = 900)
heatmap.2(matScaled_log2,
          breaks=colors,
          col=my_palette,
          dendrogram = "both",
          Colv = dend1,
          Rowv = dend2,
          scale = "row",
          labRow = c(""),
          labCol = c(""),
          # colCol = col_labels,
          # colCol = colColors4,
          ColSideColors = colColors4, 
          trace = "none",
          density.info="none",
          key.title="",
          # margins=c(12,5),
          lwid=c(1,5), #make column of dendrogram and key very small and other colum very big
          lhei=c(1,4), #make row of key and other dendrogram very small and other row big.
          hclust = function(x) hclust(x, method='ward.D2'), 
          distfun=function(x) as.dist(1-cor(t(x), method="spearman"))
)
dev.off()

