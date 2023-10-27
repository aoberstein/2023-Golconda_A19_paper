#!/usr/local/bin/Rscript
ptm <- proc.time()
rootdir = getwd()
timestamp=format(Sys.time(), "%Y-%m-%d_%H:%M:%S")
errorLog = paste0(rootdir, "/", timestamp, "_error.log")

infiles = Sys.glob( paste0(rootdir,"/*metadata*") )
# scratchDir = "/mnt/md0/duncan"
scratchDir = rootdir
print(infiles)
print(paste("rootdir:", rootdir))

for (metadataFile in infiles){
  tryCatch({
    # metadata = infiles[1]	#for testing
    infilename = basename(metadataFile) #for errorLog
    metadata = read.csv(metadataFile, header=T, sep='\t', stringsAsFactors=F)
    prefix = metadata$study[1]
    print(paste("metadata prefix:", prefix) )
    print("metadata:")
    print(metadata)
    setwd(scratchDir)
    if ( dir.exists(prefix) == "FALSE" ){
      dir.create(prefix, mode = "0775")
      }
    setwd(prefix)
    Sys.umask("002") #make sure group has all permissions on new files, so samba work is ok.
    print(paste("working dir:", getwd() ) )

    ## get fastqs (.sra files)
   sink(file = paste0(prefix, "_getFastq.log"), split = TRUE)
   fastqToRatios::getFastq(
     object = TRUE,
     metadataObject = metadata,
     sratoolkitFolder = "/opt/sratoolkit.2.10.8-ubuntu64/bin"
   )
   sink()

    ## extract fastqs
   sink(file = paste0(prefix, "_extractFastq.log"), split = TRUE)
   fastqToRatios::extractSRA(
    threads = 30,
	    # fastqDump = "/opt/parallel-fastq-dump-master/parallel-fastq-dump",
    fastqDump = "/opt/parallel-fastq-dump-0.6.6/parallel-fastq-dump",
    python3 = "/usr/bin/python3",
    tmpDir = "/mnt/nvme/tmp"
   )
   sink()

    ## detect phred encoding
   sink(file = paste0(prefix, "_detectPhred.log"), split = TRUE)
   metadata_append_phred = fastqToRatios::detectPhred(
    object = TRUE,
    metadataObject = metadata,
     threads = 30,
    python2 = "/usr/bin/python2"
  )
   write.table(file = "metadata_phred.tab", metadata_append_phred, row.names = F, sep = '\t')
   sink()


   ## detect strandedness
  sink(file = paste0(prefix, "_detectStrandedness.log"), split = TRUE)
  metadata_append_strand = fastqToRatios::detectStrand(
    metadataFile = "metadata_phred.tab",
    object = "FALSE",
    threads = 30,
    hisat_indexes = paste0(
      "/home/adam/database/genomes/rna-seq/reference_genomes/built_indexes/hisat2/2020-04-21_grch38-Merlin"
    ),
    hisat_prefix = "grch38-merlin",
#     python2 = "/opt/Python-2.7.18/build/bin/python",
    python2 = "/usr/bin/python2",
    bed = "/home/adam/database/beds/hg38_RefSeq.bed"
  )
    write.table(
      file = "metadata_phred_strand.tab",
      metadata_append_strand,
      sep = "\t",
      row.names = FALSE
    )
    sink()

# collapse technical replicates (may not be necessary)
   metadata = fastqToRatios::collapseTechReps(
    object = FALSE,
     metadataFile = "metadata_phred_strand.tab"
   )
   print(metadata)
   write.table(
    file = "metadata_phred_strand_rep_collapsed.tab",
    metadata,
    sep = '\t',
     row.names = F
   )



    ## align with Kallisto
    sink(file=paste0(prefix, "_Kallisto_with_reps.log"), split=TRUE)
    fastqToRatios::kallisto(
      fastq_path = ".",
      object = FALSE,
      metadataFile = "metadata_phred_strand_rep_collapsed.tab",
      threads = 30,
      pseudobam = FALSE,
      bootstrap_samples = 30,
      idx = paste0(
        "/home/adam/projects/2020-08-16_HCMV_M5_A19_1315_non-canonical/1-RP_kallisto_database/",
	"vector_TB40-nonCanonical_TB40_gencode34.idx")
    )
    sink()

   # ## cleanup large files
   #file.remove(Sys.glob("*.sra"))
   # file.remove(Sys.glob("*SAMPLE*"))
   # file.remove(Sys.glob("*.fastq"))
   # file.remove(Sys.glob("*.fastq.gz"))
   # file.remove(Sys.glob("*.sam"))

    ## reset working dir to rootdir     
    #setwd(rootdir) 
    cat("loop cycle COMPLETE\n")
   
    }, error=function(e){ 
     #cat("ERROR :",conditionMessage(e), "\n") 
     a=paste0("ERROR processing ", infilename, ":") 
     write( paste(a, toString(e), sep='\n' ), errorLog, append=TRUE )
       } )
}
total_time<-proc.time() - ptm
print(total_time)

