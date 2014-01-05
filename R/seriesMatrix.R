# Functions for reading GEO series matrix files
#

seriesMatrixSampleHeader <- function(filename){
  f <- file(filename, open="r")
  header <- readLines(f,n=1000)
  header_length <- grep("!series_matrix_table_begin",header)
  header <- header[1:(header_length-1)]
  
  sample_header <- header[grepl("!Sample_",header)]
  sample_header_m <- matrix(unlist(strsplit(sample_header,"\t")),ncol=length(sample_header))
  sample_header_header <- gsub('^!','',sample_header_m[1,])
  sample_header_df <- as.data.frame(gsub('"','',sample_header_m[-1,]), stringsAsFactors=F )
  names(sample_header_df) <- sample_header_header
  close(f)
  
  return(sample_header_df)
}

seriesMatrixRead <- function(filename){
  tbl <- read.table(file=filename,sep="\t",blank.lines.skip=T,comment.char="!",header=T, na.strings=c("Null","null"))
  row.names(tbl) <- tbl$ID_REF
  tbl <- tbl[,names(tbl)!="ID_REF"]
  return(tbl)
}