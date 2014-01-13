# Main script/pipeline for generating MI/CLR files from raw series matrix files from GEO
#

source("R/seriesMatrix.R") # Functions for reading GEO series matrix files
source("R/probeMapping.R") # Mapping probes to genes using BLAST result
source("R/platformMerge.R") # Functions for merging platform matrices
source("R/MiClrScripter.R") # Generate script to run MI/CLR program
source("R/slurm.R") # generate header for SLURM batch script
# source("http://bioconductor.org/biocLite.R")
# biocLite("preprocessCore")
require("preprocessCore") # normalize.quantiles


########
# make tables of series files with corresponding platforms and species
# TODO: do this a different way

cat("\nReading headers:")
GSE_files <- system("find input_data/* -name GSE*_series_matrix.txt",intern=T)
platforms <- character(0)
for( GSE_file in GSE_files){
  cat(GSE_file,"... ")
  sample_header <- seriesMatrixSampleHeader(GSE_file)
  platform <- unique(sample_header$Sample_platform_id)
  if(is.null(platform)){ # emergency handling if header format is incompatible
    cat("HEADER FORMAT ERROR ... ")
    platform <- gsub('"','',strsplit(system(paste("grep !Sample_platform_id",GSE_file),intern=T),"\t")[[1]][2])
  }
  cat(platform,"channels:",sample_header$Sample_channel_count[1],"\n")
  platforms <- c(platforms, platform)
}

# these are the platforms to be used
barleyPlatform <- "GPL1340"
ricePlatforms <- c( "GPL14568", "GPL16106", "GPL2025","GPL7344" )
maizePlatforms <- c( "GPL1992", "GPL1993", "GPL4032", "GPL6438", "GPL8162", "GPL8163" )

seriesFileTbl <- data.frame( platform=platforms,file=GSE_files)
seriesFileTbl$species[seriesFileTbl$platform %in% barleyPlatform] <- "barley" 
seriesFileTbl$species[seriesFileTbl$platform %in% ricePlatforms] <- "rice" 
seriesFileTbl$species[seriesFileTbl$platform %in% maizePlatforms] <- "maize"

# ignore series GSE6990 because of invalid file format
ignoreFiles <- c("input_data/barley/whole_series_matrix_files/GSE6990_series_matrix.txt")

# remove files that are not to be used
seriesFileTbl <- seriesFileTbl[ !(is.na(seriesFileTbl$species) | seriesFileTbl$file %in% ignoreFiles) ,]


###############################
# Generate a vector of blast files and name withthe corresponding platforms

blastFile <- c( system("find data/rice/blast_results -name *.blastn",intern=T),
                system("find data/maize/blast_results -name *.blastn",intern=T),
                system("find data/barley/blast_results -name *.blastn",intern=T) )

# get corresponding platform from the filename
names(blastFile) <- substr(blastFileTbl$file,regexpr("GPL[0-9]+",blastFileTbl$file), 
                           regexpr("_[^/]+\\.blastn",blastFileTbl$file)-1)


# Merge, log transform if necessary, and normalize series files (from same platform)
mergeSeriesMatrices <- function(seriesFiles){
  mergedMatrix <- NULL
  for(seriesFile in seriesFiles){
    cat("\nReading series matrix file",seriesFile,"...")
    m <- as.matrix(seriesMatrixRead(seriesFile))
    
    # do log transform if highest value is larger than 100
    if( max(m,na.rm=T) > 100 ){
      cat("log2 transform...")
      m <- log2(m)
    }
    
    
    if( is.null(mergedMatrix) )
      mergedMatrix <- m
    else
      mergedMatrix <- cbind(mergedMatrix ,m)
  }
  #     normalize
  cat("\nApplying quantile normalization on probe values...")
  mergedMatrixNorm <- normalize.quantiles(mergedMatrix)
  rownames(mergedMatrixNorm) <- rownames(mergedMatrix)
  colnames(mergedMatrixNorm) <- colnames(mergedMatrix)
  cat("done\n")
  return(mergedMatrixNorm)
}

probeStats <- list() #use this to store some statistics on the probe mapping

# for each species
for( species in unique(seriesFileTbl$species)){
  speciesGeneMatrix <- null
  # for each platform
  for( platform in unique(seriesFileTbl$platform[ seriesFileTbl$species == species ])){
    
    # merge all series matrices for the platform
    cat("\nMerging series matrices for platform",platform,"\n")
    platformProbeMatrix <- mergeSeriesMatrices( seriesFileTbl$file[seriesFileTbl$platform==platform] )
    
    # map probes to genes
    cat("\nReading probe mapping...")
    probe2genes <- getProbeMapping(blastFile[platform], removeIsoformsFun[[species]])
    cat("\nConverting probe values to gene values...")
    platformGeneMatrix <- convertProbe2GeneValues(platformProbeMatrix, probe2genes)
    cat("done\n")
    
    # collect interresting statistics
    probeStats[[species]] <- list( probesPerGene = table(table(unlist(probe2genes))),
                                   genesPerProbe = table(unlist(lapply(probe2genes, length)))
                                   probesMapped = table(row.names(inM)%in%names(probe2genes)) )
    
    # merge matrices for
    if( is.null(speciesGeneMatrix) ){
      speciesGeneMatrix <- platformGeneMatrix
    } else {
      cat("\nAdding platform matrix",platform,"to",species,"matrix...")
      speciesGeneMatrix <- mergeExpressionMatrix( speciesGeneMatrix, platformGeneMatrix)
      cat("done\n")
    }
  }
  
  # remove rows with too little data
  speciesGeneMatrix <- removeNArows( speciesGeneMatrix )
  
  # normalize species matrix
  cat("\nNormalizing species matrix...")
  speciesGeneMatrix <- normalize.quantiles( speciesGeneMatrix )
  
  # save gene expression matrix
  geneMatrixFile <- file.path("data",species,"geneMatrix.txt")
  cat("done\nSaving to file",outFile)
  write.table( normalize.quantiles( removeNArows( speciesGeneMatrix ) ),
               file=outFile,quote=F,row.names=T,col.names=T,sep="\t")
  
  
  # generate MI/CLR script
  outdir  <- file.path( "data",species)
  script <- generateBatchScript( jobName = paste("mi_",species,sep=""), 
                                 command = createMICLRscript(geneMatrixFile, outdir) )
  cat(script, file = file.path("scripts",paste("miCalc_",species,".sh",sep="")))  
}

save(probeStats, file = "output/probeMappingStats.RData")