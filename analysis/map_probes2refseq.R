#
# Map probes to genes using the blast results
#


blastCols <- unlist(strsplit("qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore"," "))
Os_files <- system("find data/rice/blast_results -name *.blastn",intern=T)

allTables <- list()

filename <- Os_files[1]
for(filename in Os_files){
  
  bt <- read.table(filename, header=F,sep="\t",col.names=blastCols, stringsAsFactors=F)
  
  # # Check if target sequence matches multiple genome sequences:
  # table(table(bt$qseqid))
  # # Check if multiple target sequences matches same genome sequence:
  # table(table(bt$sseqid))
  
  # remove isoform numbers from sseqid
  bt$sseqid_noiso <- sub("[.][^.]*$", "", bt$sseqid, perl=TRUE)
  # onbly the ref seq id starting with "LOC_" has isoform numbers. Restore the rest:
  bt$sseqid_noiso[!grepl("^LOC_",bt$sseqid)] <- bt$sseqid[!grepl("^LOC_",bt$sseqid)]
  
  # pick the best hit for each query
  btUnique <-bt[ bt$qseqid != c("xxx", bt$qseqid[1:(length(bt$qseqid)-1)]),]
  
  
  allTables[[filename]] <- btUnique$sseqid_noiso
}

n <- length(Os_files)
platformName <- substr(Os_files,regexpr("GPL[0-9]+",Os_files), regexpr("_rice.all.cds.blastn$",Os_files)-1)

uniGene <- unique(unlist(allTables))
inDataSet <- data.frame(row.names=uniGene)
for( i in 1:n ){
  inDataSet[,platformName[i]] <- uniGene %in% unique( allTables[[i]] )
}

barplot(table( rowSums(inDataSet)),xlab="Number of platforms containing same gene")


library(xtable)
X <- sort(colSums(inDataSet),decreasing=T)
X <- data.frame(row.names=names(X), genes=as.integer(X))
print(xtable(X),type="html")


