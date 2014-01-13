# Mapping probes to genes using BLAST result
#

getProbeMapping <- function(blastFile, removeIsoformsFun){
  blastCols <- unlist(strsplit("qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore"," "))
  bt <- read.table(blastFile, header=F,sep="\t",col.names=blastCols, stringsAsFactors=F)
  
  bt$sseqid_noiso <- removeIsoformsFun(bt$sseqid)
  
  # Find best hits
  # Use the first hit per query (BLAST sorts best hits first) and include additional 
  # hits if they have equal bitscore.
  isBest <- function(bt){
    prevProbe <- ""
    topScore <- 0
    best <- rep(F,dim(bt)[1])
    for(i in 1:dim(bt)[1]){
      if(prevProbe != bt$qseqid[i]){
        prevProbe <- bt$qseqid[i]
        best[i] <- T
        topScore <- bt$bitscore[i]
      } else if( bt$bitscore[i] == topScore ){
        best[i] <- T
      }
    }
    return(best)
  }
  
  bt <- bt[isBest(bt),] # keep only best hits
  # generate list of best hit genes per probe
  probe2genes <- tapply(bt$sseqid_noiso, bt$qseqid, unique)
  
  return(probe2genes)
}

convertProbe2GeneValues <- function(inM, probe2genes){
  genes <- unique(unlist(probe2genes))
  probes <- names(probe2genes)
  
  # get number of probes per gene. (expression values shall be divided by this)
  probesPerGene <- table(unlist(probe2genes))
  
  outM <- matrix(0,nrow=length(genes),ncol=ncol(inM))
  rownames(outM) <- genes
  colnames(outM) <- colnames(inM)
  
  # for each probe:
  #   for each gene matching the probe:
  #     gene expression = probe value / number of probes for that gene
  for(probe in probes){
    probeGenes <- probe2genes[[probe]]
    outM[probeGenes,] <- outM[probeGenes,] + as.matrix(1/probesPerGene[probeGenes])%*%t(inM[probe,])
  }
  return(outM)
}

removeIsoformsFun <- list(  
  rice = function(geneId){  
    # remove isoform numbers from geneId
    geneId_noiso <- sub("[.][^.]*$", "", geneId, perl=TRUE)
    # only the gene id's starting with "LOC_" has isoform numbers. Restore the rest:
    geneId_noiso[!grepl("^LOC_",geneId)] <- geneId[!grepl("^LOC_",geneId)]
    return(geneId_noiso)
  },
  
  maize = function(geneId){  
    # remove isoform numbers from geneId
    geneId_noiso <- sub("_T[0-9]+$","", geneId, perl=TRUE)
    return(geneId_noiso)
  },
  
  barley = function(geneId){  
    # no isoform numbers in barley
    return(geneId)
  }
)