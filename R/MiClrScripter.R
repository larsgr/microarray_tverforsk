# Generate script to run MI/CLR program
#

###### Script template
# Fields: %INPUTFILEPATH% - path to gene expression matrix file
#         %OUTPUTPATH% - path to store resulting MI and CLR file 

template_script <- "
# create temporary directory
if [ ! -d %OUTPUTPATH% ]; then
  mkdir %OUTPUTPATH%
fi
if [ ! -d %OUTPUTPATH%/temp ]; then
  mkdir %OUTPUTPATH%/temp
fi

# convert NA to magic number 1e+30
sed 's/NA/1e+30/g' %INPUTFILEPATH% > %OUTPUTPATH%/temp/microarray_file.txt

# get ngenes
NGENES=$(cat %OUTPUTPATH%/temp/microarray_file.txt | wc -l)

# get nexp
NEXP=$(sed -n 1p %OUTPUTPATH%/temp/microarray_file.txt | wc -w)

# remove possible old temporary mi file
if [ -f %OUTPUTPATH%/temp/mi_file.txt ]; then
  rm %OUTPUTPATH%/temp/mi_file.txt
fi

# run mi program
bin/genepair 6 %OUTPUTPATH%/temp/microarray_file.txt %OUTPUTPATH%/temp/mi_file.txt $NGENES $NEXP 0 $NGENES 7 3

# check if correct number of lines in mi_file
if [[ $(cat %OUTPUTPATH%/temp/mi_file.txt | wc -l) != $(expr $NGENES - 1) ]]; then
  echo ERROR! Incorrect number of lines in MI file
  exit 1
fi

# remove possible old mi file
if [ -f %OUTPUTPATH%/mi_file.txt ]; then
  rm %OUTPUTPATH%/mi_file.txt
fi

# move temp mi file to mi file
mv %OUTPUTPATH%/temp/mi_file.txt %OUTPUTPATH%/mi_file.txt

# remove temporary microarray file
rm %OUTPUTPATH%/temp/microarray_file.txt

# run clr program
bin/genepair 2 %OUTPUTPATH%/mi_file.txt $NGENES %OUTPUTPATH%/clr_file.txt 1
"



createMICLRscript <- function(geneExpressionFile, outdir){
  script <- gsub("%INPUTFILEPATH%",geneExpressionFile,template_script)
  script <- gsub("%OUTPUTPATH%",outdir,script)
  return( script )
}
