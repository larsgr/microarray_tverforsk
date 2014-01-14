microarray_tverforsk
====================

This repository contains scripts used as part of a research project to analyse microarray expression data. Read about the project at [grassrewired.wordpress.com](http://grassrewired.wordpress.com)

The scripts are used to merge public expression data from several platforms into a single expression matrix. This expression data is then used to generate MI/CLR matrices.

Since different microarrays use different probes it is necessary to map the probes to a common set of genes. This is done by BLASTing the probe's sequences against the reference genome. This is performed by the shell script **scripts/blast_target_sequences.sh**. As input, you need the reference genome CDS sequences and the probe's target sequences in fasta format (not included in repository). The blast result files are used by the next script **scripts/GEOseriesToCLR.R**. This script takes expression data in the GEO series matrix file format as input and outputs a gene expression matrix for each species and a set of shell scripts that are ready to run on SLURM. These generated batch scripts will run the MI/CLR program.


Folder structure
----------------

* **/input_data** - input data (gene expression)
* **/scripts** - main scripts
* **/data** - generated data
* **/R** - helper functions in R
* **/analysis** - scripts used to analyse the data (typically outputs figures or tables)
* **/output** - figures and tables from the analysis
* **/bin** - binary executables (actually only the MI/CLR program)