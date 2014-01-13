# Functions for merging platform matrices
#

# Merge two gene expression matrixes
mergeExpressionMatrix <- function(m1,m2){
  # make new matrix the can contain all columns and combined rows of m1 and m2
  combinedRownames <- unique( c(rownames(m1), rownames(m2)) )
  m <- matrix(NA,nrow=length(combinedRownames),ncol=(ncol(m1)+ncol(m2)))
  rownames(m)<- combinedRownames  
  colnames(m)<- c(colnames(m1), colnames(m2))
  
  # sanity checks:
  if( length(unique(c(colnames(m1), colnames(m2)))) != ncol(m) ) stop("duplicate sample names")
  if( length(unique(rownames(m1))) != nrow(m1) ) stop("duplicate genes")      
  if( length(unique(rownames(m2))) != nrow(m2) ) stop("duplicate genes")
  
  m[rownames(m1),colnames(m1)] <- m1
  m[rownames(m2),colnames(m2)] <- m2
  return(m)
}

# Load and merge expression matrixes from a list of files
mergeFiles <- function(GEfiles){
  m <- NULL
  cat("\nMerging matrix files:")
  for( file_name in GEfiles){  
    cat("\n+ ",file_name)
    newM <- as.matrix(read.table(file_name))
    if(is.null(m)){
      m <- newM
    } else {
      m <- mergeExpressionMatrix( m, newM)
    }
  }
  cat("\nMerging complete.\n")
  return(m)
}

# remove genes with too few data points
removeNArows <- function( m, minData = 20 ){
  dataPerRow <- apply(m,1,function(X){sum(!is.na(X))})
  cat( "\nRemoving", sum(dataPerRow < minData),"of",nrow(m),"rows\n" )
  return( m[dataPerRow >= minData, ] )
}
