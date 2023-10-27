#!/home/adam/programs.installed/R-3.6.3/build/bin/Rscript
library("viridis", quietly = TRUE)
library("ggplot2", quietly = TRUE)
library("dplyr", quietly = TRUE)
library("reshape2", quietly = TRUE)
options("scipen"=100, "digits"=4)

# User defined variables --------------------------------------------------
ttgFile = "/home/adam/projects/2020-08-16_HCMV_M5_A19_1315_non-canonical/ttg.tab"
metadataFile = "/home/adam/projects/2022-10-04_EMT-paper/finished/metadata_merged.tab"
countFile = paste0("/mnt/archive/adam_storage/2022-10-04_EMT-paper/finished/2-limmaGene_TMM/"
                   ,"GeTMM_allData_TMM_expValues_cpm.tab")

# Capture command line arguments ------------------------------------------
# args <- commandArgs(trailingOnly = TRUE)
# prefix = unlist(strsplit(as.character(args[1]), ".txt"))
# cat(prefix)
# cat("\n")
# y_limit = as.numeric(args[2])


# # Testing -- mock command line arguments ----------------------------------
infile = "EMT_MET_A19_paper.txt"
prefix = unlist(strsplit(as.character(infile), ".txt"))
cat(prefix)
y_limit = 14

## variables
plot_title = paste0(prefix,"_limma_normalized_log2cpm.pdf")
jitterPlotTitle = paste0(prefix,"_jitterPlot.pdf") 
ttg = read.csv(ttgFile, header = T, sep = '\t', stringsAsFactors = F)
# mock-in ensembl ids for genes with no hgnc (these are custom transcripts)
counts_preFilt = read.csv(countFile, sep='\t', header= TRUE, stringsAsFactors = FALSE)
metadata_preFilt = read.csv(metadataFile, sep='\t', header= TRUE, stringsAsFactors = FALSE)
# metadata = metadata[naturalsort::naturalorder(metadata$run),]

metadata = metadata_preFilt
counts = counts_preFilt

## genes of interest
ids = read.csv( paste0(prefix,".txt"), header=F, stringsAsFactors=F)
names(ids)[1] = "target_id"

## subset count data and append metadata
gids = ifelse(ids$target_id %in% ttg$ens_gene, ids$target_id, ttg$ens_gene[match(ids$target_id, ttg$ext_gene)]  )
conversionTable = cbind(ids, gids)
pdat = subset(counts, counts$gid %in% as.character(gids) )
pdat_long = melt(pdat)
names(pdat_long) = c("target_id", "sample", "cpm")

# condition = metadata$treat[match(pdat_long$sample, metadata$experiment)]
cell = metadata$cell[match(pdat_long$sample, metadata$run)]
run = metadata$run[match(pdat_long$sample, metadata$run)]
run2 = gsub("__.*", "", run)
treatment = metadata$treament[match(pdat_long$sample, metadata$run)]
other = metadata$other[match(pdat_long$sample, metadata$run)]
treat = paste(cell, run2, treatment, other, sep = '|')
pdat2 = cbind(pdat_long, cell, run2, treatment, other, treat)

# create facet labels column
labs = conversionTable$target_id[match(pdat2$target_id, conversionTable$gids)]
plot_data = cbind(labs, pdat2)
plot_data$cell = factor(plot_data$cell, levels = c("RWPE-1", "MCF10A", "ARPE-19",  "MRC5"))
# data4$treatment <- factor(data4$treatment, levels = treatments   )

## compute n valuesfor each group
n_A19 = length(unique(subset(metadata$run, metadata$cell == "ARPE-19")))
n_10A = length(unique(subset(metadata$run, metadata$cell == "MCF10A")))
n_M5 = length(unique(subset(metadata$run, metadata$cell == "MRC5")))
n_RWPE = length(unique(subset(metadata$run, metadata$cell == "RWPE-1")))
s_A19 = length(unique(subset(metadata$study, metadata$cell == "ARPE-19")))
s_10A = length(unique(subset(metadata$study, metadata$cell == "MCF10A")))
s_M5 = length(unique(subset(metadata$study, metadata$cell == "MRC5")))
s_RWPE = length(unique(subset(metadata$study, metadata$cell == "RWPE-1")))

## reorder primers

plot_data$labs = factor(plot_data$labs, levels = 
c("C1orf116",
"CDH1",
"CDH3",
"EPCAM",
"GJB3",
"OVOL1",
"OVOL2",
"ST14",
"MARVELD3",
"CDH2",
"CDH11",
"FBN1",
"FN1",
"VIM",
"ZEB1",
"ZEB2",
"PPIA",
"RPLP0"
))

## jitter plots
## custom colors
col2=HexLogClock=c("#addd8e") #light green
col3=HexLogClock=c("#31a354") #green 
col4=HexLogClock=c("#fee8c8") #pink
col5=HexLogClock=c("#fdbb84") #salmon
col6=HexLogClock=c("#e34a33") #red
col_orange=HexLogClock=c("#FFA800")
col_blue=HexLogClock=c("#2C54E8")
viridis_color = "E"
# colors = rep("black", length(unique(plot_data$cell) ))
# colors = c("ARPE-19" = "black", "MCF10A" = "blue", "MRC5" = "red", "RWPE-1" = "orange")
# colors = factor(c("MCF10A" = col6 , "RWPE-1" = col_orange, "MRC5" = col_blue, "ARPE-19" = "black" )
                # , levels = c("MCF10A", "RWPE-1", "MRC5", "ARPE-19"))
colors = c(col_orange, col6, "black" , col_blue)
# box_shade = 'grey70'
box_shade = 'white'
columns = 9
width = 32
# height = (ceiling(length(levels(pdat_long2$hgnc))/8) ) * 5.5
height = (ceiling(length(unique(plot_data$labs))/columns) ) * 5
# pdf(file = jitterPlotTitle, width=width, height=height)
viridis_color="E"
# colors = rep("black", length(unique(plot_data$cell) ))
plot2 = ggplot(data = plot_data, aes(x = cell, y = log2(cpm), fill = cell), color = cell) +
geom_boxplot(outlier.shape = NA, alpha=0.3, color="black", fill = box_shade,
               width=0.9, size = 0.7) +
  geom_jitter(position = position_jitter(0.15), alpha = 1, aes(fill = cell, color=cell), size=3) +
  facet_wrap(~labs, scale="free", ncol=columns) +
  #"magma" (or "A"), "inferno" (or "B"), "plasma" (or "C"), "viridis" (or "D", the default option) and "cividis" (or "E").
      # scale_color_viridis(discrete = TRUE, direction = 1, option=viridis_color) +
  scale_colour_manual(values = colors, labels = c(
              paste("RWPE-1", paste0("(n = ", n_RWPE, "; s = ", s_RWPE, ")") ),
              paste("MCF10A", paste0("(n = ", n_10A, "; s = ", s_10A, ")") ),
              paste("ARPE-19", paste0("(n = ", n_A19, "; s = ", s_A19, ")") ),
              paste("MRC5", paste0("(n = ", n_M5, "; s = ", s_M5, ")") )
                      )) +
  # scale_colour_manual(values = colors) +
  scale_y_continuous( breaks = seq(-10, 14, by = 2), limits = c(-6, y_limit) ) +
  theme_bw() + 
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme(
    # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=30, color = "black"),
    axis.text.x=element_blank(),
    strip.text=element_text(size=28),
    axis.text.y = element_text(size=28, color = "black"),
    axis.ticks = element_line(size = 1.2, color = "black"),
    axis.ticks.length.y = unit(.3, "cm"),
    axis.ticks.length.x = unit(0, "cm"),
    panel.margin = unit(1.5,"lines"),
    legend.text=element_text(size=30),
    legend.title=element_text(size=0),
    axis.title.y=element_text(size=0),
    axis.title.x=element_text(size=0),
    legend.key.size = unit(1.2, "cm"),
    panel.border = element_rect(size = 1.2),
    panel.grid = element_blank(),
    strip.background = element_blank() ) #element_rect(size = 2, fill = "none") )
# dev.off()
ggsave(plot2, file = jitterPlotTitle, width = width, height = height, unit = "in")

