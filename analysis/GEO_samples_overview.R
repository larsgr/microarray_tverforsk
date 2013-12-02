#########
#
# Gather metadata from series_matrix_files and generate overview.
#

# find all series_matrix_files
Hv_files <- system("find input_data/barley/whole_series_matrix_files -name GSE*_series_matrix.txt",intern=T)
Os_files <- system("find input_data/rice/whole_series_matrix_files -name GSE*_series_matrix.txt",intern=T)
Zm_files <- system("find input_data/maize/whole_serie_matrix_files -name GSE*_series_matrix.txt",intern=T)

# extract the metadata we need
GSEtbl <- data.frame(file=character(0),platform=character(0),series=character(0), 
                     samples=integer(0), sample_platform=character(0), has2ch=logical(0), 
                     stringsAsFactors=F)
for(file in c(Hv_files,Os_files,Zm_files)){
  f <- file(file,'r')
  lines <- readLines(f, n=1000)
  close(f)
  
  platform = paste(gsub('"','',unlist(lapply(strsplit(lines[grep("^!Series_platform_id",lines)],"\t"),"[",2))),collapse=" ")
  series = gsub('"','',unlist(lapply(strsplit(lines[grep("^!Series_geo_accession",lines)],"\t"),"[",2)))
  
  samples = unlist(lapply(strsplit(lines[grep("^!Sample_geo_accession",lines)],"\t"),function(X){ length(X)-1  }))
  sample_platform = paste(gsub('"','',unlist(lapply(strsplit(lines[grep("^!Sample_platform_id",lines)],"\t"),"[",2))),collapse=" ")
  has2ch <- sum(grepl("^!Sample_taxid_ch2",lines))>0
  
  GSEtbl <- rbind(GSEtbl,data.frame(file,platform,series,samples, sample_platform, has2ch, stringsAsFactors=F))
}


GSEtbl$species <- c(rep("barley",length(Hv_files)),rep("rice",length(Os_files)),rep("maize",length(Zm_files)))

# # remove duplicates (There is only one series matrix files per platform per series.
# # Series metadata is therefore duplicated if there are several platforms in one series)
# GSEtbl <- GSEtbl[match(unique(GSEtbl$series),GSEtbl$series),]

# order
GSEtbl <- GSEtbl[order(GSEtbl$species,GSEtbl$platform,GSEtbl$series),]

# summarize
sapply(c("rice","barley","maize"),function(species){
  nsamples <- sum(GSEtbl$samples[GSEtbl$species == species])
  nplatforms <- length(unique(unlist(strsplit(GSEtbl$platform[GSEtbl$species == species]," "))))
  nseries <- length(unique(GSEtbl$series[GSEtbl$species == species]))
  return(data.frame(nseries,nplatforms,nsamples))
}) -> GSEtblSum

# per platform summary
# number of samples (1-channel and 2-channel)
GSEtbl$samples_1ch <- ifelse( GSEtbl$has2ch, 0, GSEtbl$samples)
GSEtbl$samples_2ch <- ifelse( GSEtbl$has2ch, GSEtbl$samples, 0)
platformSummary <- data.frame( species = tapply( GSEtbl$species, GSEtbl$sample_platform, "[",1),
                               nsamples1ch = tapply( GSEtbl$samples_1ch, GSEtbl$sample_platform, sum),
                               nsamples2ch = tapply( GSEtbl$samples_2ch, GSEtbl$sample_platform, sum) )



####
#
# Print the tables in html format
#

library(xtable)
print(xtable(GSEtblSum),type="html")


print(xtable(GSEtbl[,c("species","series","platform")]),type="html")
