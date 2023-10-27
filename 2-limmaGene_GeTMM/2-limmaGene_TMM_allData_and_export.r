#!/home/adam/programs.installed/R-3.6.3/build/bin/Rscript
options(show.error.locations = TRUE)
library("limma")
library("edgeR") #for DGEList function
library("RColorBrewer") #for limma filter diagnostic plot
library("tximport")

## PER RUN VARIABLES
tx2geneFile = paste0("../tx2gene_vector_TB40-nonCanonical_TB40_gencode34.transcripts")
pwd = getwd()
kallistoRootDir1 = ".."
metadataFile1 = paste0(kallistoRootDir1, "/metadata_merged.tab")
metadata1 = read.csv(metadataFile1, sep='\t', header= TRUE, stringsAsFactors = FALSE)

## import read counts using tximport
studies = paste(metadata1$study, metadata1$run, sep ="/")
files = file.path(kallistoRootDir1, studies, "abundance.tsv")
names(files) = c(metadata1$run)

# 
tx2gene = read.csv(tx2geneFile, sep = '\t', header = TRUE, stringsAsFactors = FALSE)

txiTX = tximport::tximport(files, type = "kallisto",
                           txIn = TRUE,
                           txOut = TRUE,
                           countsFromAbundance = "no",
                           tx2gene = tx2gene)

txiGene = tximport::tximport(files, type = "kallisto",
                             txIn = TRUE,
                             txOut = FALSE,
                             countsFromAbundance = "no",
                             # countsFromAbundance = "scaledTPM",
                             tx2gene = tx2gene)

save(file = "txiTX.RData", txiTX)
save(file = "txiGene.RData", txiGene)
 # load("txiTX.RData")
 # load("txiGene.RData")
# 
# 
# RPK normalize genes -----------------------------------------------
rawCounts = as.matrix(txiGene$counts, stringsAsFactors = F)

# tximport gives a matrix of average transcript lengths for each gene
lengths = txiGene$length

  # my recipe
  rpk = rawCounts/lengths
  mData = metadata1
  treatment = mData$treament
  runs = mData$run
  sample = mData$experiment
  cell = mData$cell
  prefix = "GeTMM_allData"

  ##add counts to DGEList object
  x = edgeR::DGEList(counts = rpk)

  #add metadata columns to samples in dge object
  x$samples$cell = cell
  x$samples$treat = treatment
  # x$samples$time = time

  ## filter HCMV genes
  # HCMVgenes = ttg$ens_gene[which(grepl("viral",ttg$target_id, ))]
  HCMVgenes = tx2gene$gid[which(grepl("viral", tx2gene$tid ))]
  # matFilt1 = subset(mat, !(row.names(mat) %in% HCMVgenes))
  keep.exprs <- !(row.names(x$counts) %in% HCMVgenes)
  filt1 <- x[keep.exprs,, keep.lib.sizes=FALSE]
  
  # ## mitochondrial genes
  # keep.exprs <- !(grepl("MT-.*", row.names(filt$counts)))
  # filt2 <- filt1[keep.exprs,, keep.lib.sizes=FALSE]
  
  # cpm = edgeR::cpm(x)
  # filtVar = 30
  filtVar = length(mData$run)
  keep.exprs <- rowSums(cpm(filt1$counts)) >= filtVar
  filt <- filt1[keep.exprs,, keep.lib.sizes=FALSE]

  ## normalize filtered counts w/ TMM method
  normTMM <- calcNormFactors(filt, method = "TMM")
  normTMM$samples$norm.factors
  write(file=paste(prefix, "limma_TMM_normFACTORS.txt", sep="_") ,normTMM$samples$norm.factors)

  # RLE
  normRLE <- calcNormFactors(filt, method = "RLE")
  normRLE$samples$norm.factors
  write(file=paste(prefix, "limma_RLE_normFACTORS.txt", sep="_") ,normRLE$samples$norm.factors)

  # TMMwsp
  normTMMwsp <- calcNormFactors(filt, method = "TMMwsp")
  normTMMwsp$samples$norm.factors
  write(file=paste(prefix, "limma_TMMwsp_normFACTORS.txt", sep="_") ,normTMMwsp$samples$norm.factors)

  # upperquartile
  normUpperquartile <- calcNormFactors(filt, method = "upperquartile")
  normUpperquartile$samples$norm.factors
  write(file=paste(prefix, "limma_upperquartile_normFACTORS.txt", sep="_") ,normUpperquartile$samples$norm.factors)

  # transform data for plotting
  raw_cpm = edgeR::cpm(x)
  raw_lcpm = log2(edgeR::cpm(x))

  filt_cpm = edgeR::cpm(filt)
  filt_lcpm = log2(edgeR::cpm(filt))

  normTMM_cpm = edgeR::cpm(normTMM)
  normTMM_lcpm = log2(edgeR::cpm(normTMM))

  normRLE_cpm = edgeR::cpm(normRLE)
  normRLE_lcpm = log2(edgeR::cpm(normRLE))

  normTMMwsp_cpm = edgeR::cpm(normTMMwsp)
  normTMMwsp_lcpm = log2(edgeR::cpm(normTMMwsp))

  normUpperquartile_cpm = edgeR::cpm(normUpperquartile)
  normUpperquartile_lcpm = log2(edgeR::cpm(normUpperquartile))


  ## diagnostic plots
  nsamples <- ncol(x)
  col = seq(1:nsamples)
  pdf(file = paste(prefix, "limmaDENSITY_PLOTS.pdf", sep="_"), width = 12, height = 8)
  par(mfrow = c(2,3))

  # unfiltered
  plot(density(raw_lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.3), xlim=c(-10,15), las=2,main="", xlab="")
  title(main="A. Raw data", xlab="Log-cpm")
  abline(v=0, lty=3)
  for (i in 2:nsamples){
    den <- density(raw_lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright", rownames(x$samples), text.col=col, bty="n")

  # filtered
  plot(density(filt_lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.3), xlim=c(-10,15), las=2,
       main="", xlab="")
  title(main="B. Filtered", xlab="Log-cpm")
  abline(v=0, lty=3)
  for (i in 2:nsamples){
    den <- density(filt_lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright", rownames(x$samples), text.col=col, bty="n")

  # TMM normalized
  plot(density(normTMM_lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.3), xlim=c(-10,15), las=2,
       main="", xlab="")
  title(main="B. TMM", xlab="Log-cpm")
  abline(v=0, lty=3)
  for (i in 2:nsamples){
    den <- density(normTMM_lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright", rownames(x$samples), text.col=col, bty="n")

  # RLE normalized
  plot(density(normRLE_lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.3), xlim=c(-10,15), las=2,
       main="", xlab="")
  title(main="B. RLE", xlab="Log-cpm")
  abline(v=0, lty=3)
  for (i in 2:nsamples){
    den <- density(normRLE_lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright", rownames(x$samples), text.col=col, bty="n")

  # TMMwsp normalized
  plot(density(normTMMwsp_lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.3), xlim=c(-10,15), las=2,
       main="", xlab="")
  title(main="B. TMMwsp", xlab="Log-cpm")
  abline(v=0, lty=3)
  for (i in 2:nsamples){
    den <- density(normTMMwsp_lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright", rownames(x$samples), text.col=col, bty="n")

  # Upperquartile normalized
  plot(density(normUpperquartile_lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.3), xlim=c(-10,15), las=2,
       main="", xlab="")
  title(main="B. Upperquartile", xlab="Log-cpm")
  abline(v=0, lty=3)
  for (i in 2:nsamples){
    den <- density(normUpperquartile_lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright", rownames(x$samples), text.col=col, bty="n")

  dev.off()


  ## boxplots of un-normalized vs normalized data
  pdf(file = paste(prefix, "limmaNORM_BOXPLOTS.pdf"), width = 12, height = 10)
  par(mfrow = c(2,3), mar=c(15, 2, 4, 2))

  boxplot( filt_lcpm, las = 2, col = col, main = "")
  title(main="Pre-normalization", ylab="Log-cpm")

  boxplot( normTMM_lcpm, las = 2, col = col, main = "")
  title(main="TMM", ylab="Log-cpm")

  boxplot( normRLE_lcpm, las = 2, col = col, main = "")
  title(main="RLE", ylab="Log-cpm")

  boxplot( normRLE_lcpm, las = 2, col = col, main = "")

  boxplot( normTMMwsp_lcpm, las = 2, col = col, main = "")
  title(main="TMMwps", ylab="Log-cpm")

  boxplot( normUpperquartile_lcpm, las = 2, col = col, main = "")
  title(main="Upperquartile", ylab="Log-cpm")

  dev.off()

  ## write table of normalized expression values
  outData = data.frame(row.names(normTMM), cpm(normTMM) )
  names(outData)[1] = c("gid")
  write.table(file = paste0(prefix, "_TMM_expValues_cpm.tab"), outData, sep = '\t', row.names = F )

