#########
#
# Gather metadata from series_matrix_files and generate overview.
#

# find all series_matrix_files
Hv_files <- system("find input_data/barley/whole_series_matrix_files -name GSE*_series_matrix.txt",intern=T)
Os_files <- system("find input_data/rice/whole_series_matrix_files -name GSE*_series_matrix.txt",intern=T)
Zm_files <- system("find input_data/maize/whole_serie_matrix_files -name GSE*_series_matrix.txt",intern=T)

# extract the metadata we need
GSEtbl <- data.frame(file=character(0),platform=character(0),series=character(0), samples=character(0), stringsAsFactors=F)
for(file in c(Hv_files,Os_files,Zm_files)){
  f <- file(file,'r')
  lines <- readLines(f, n=100)
  close(f)
  
  platform = paste(gsub('"','',unlist(lapply(strsplit(lines[grep("^!Series_platform_id",lines)],"\t"),"[",2))),collapse=" ")
  samples = gsub('"','',unlist(lapply(strsplit(lines[grep("^!Series_sample_id",lines)],"\t"),"[",2)))
  series = gsub('"','',unlist(lapply(strsplit(lines[grep("^!Series_geo_accession",lines)],"\t"),"[",2)))
  
  GSEtbl <- rbind(GSEtbl,data.frame(file,platform,series,samples,stringsAsFactors=F))
}


GSEtbl$species <- c(rep("barley",length(Hv_files)),rep("rice",length(Os_files)),rep("maize",length(Zm_files)))

# remove duplicates (There is only one series matrix files per platform per series.
# Series metadata is therefore duplicated if there are several platforms in one series)
GSEtbl <- GSEtbl[match(unique(GSEtbl$series),GSEtbl$series),]
# order
GSEtbl <- GSEtbl[order(GSEtbl$species,GSEtbl$platform,GSEtbl$series),]

# summarize
sapply(c("rice","barley","maize"),function(species){
  nsamples <- length(unlist(strsplit(GSEtbl$samples[GSEtbl$species == species]," ")))
  nplatforms <- length(unique(unlist(strsplit(GSEtbl$platform[GSEtbl$species == species]," "))))
  nseries <- sum(GSEtbl$species == species)
  return(data.frame(nseries,nplatforms,nsamples))
}) -> GSEtblSum

####
#
# Print the tables in html format
#

library(xtable)
print(xtable(GSEtblSum),type="html")


print(xtable(GSEtbl[,c("species","series","platform")]),type="html")
