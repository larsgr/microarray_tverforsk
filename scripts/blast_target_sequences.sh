#
# generate blast database for each species
#

DATAPATH="/mnt/users/mariansc/microarray_tverforsk"

BLASTPARAM="-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore' -num_threads 12 -perc_identity 95"


# for rice:
mkdir data/rice
mkdir data/rice/db
makeblastdb -in $DATAPATH/rice/annotation_fastas/all.cds -dbtype nucl -out data/rice/db/all.cds

# for maize:
mkdir data/maize
mkdir data/maize/db
makeblastdb -in $DATAPATH/maize/annotation_fastas/ZmB73_5a_WGS_cds.fasta -dbtype nucl -out data/maize/db/WGS_cds

# for barley
mkdir data/barley
mkdir data/barley/db
makeblastdb -in $DATAPATH/barley/annotation_fastas/barley_HighConf_genes_MIPS_23Mar12_CDSSeq.fa -dbtype nucl -out data/barley/db/all.cds

####
# Run blast

# for rice

mkdir data/rice/blast_results

for TARGET_FASTA in $(find $DATAPATH/rice/target_fastas/ -name GPL*)
do
  GPL=$(echo $TARGET_FASTA | sed 's/.*\(GPL[0-9]*\).*/\1/')
  BLASTRES="data/rice/blast_results/"$GPL"_rice.all.cds.blastn"
  blastn -query $TARGET_FASTA \
         -db data/rice/db/all.cds \
         -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore' -num_threads 12 -perc_identity 95 \
         -out $BLASTRES
done

# for maize

mkdir data/maize/blast_results

for TARGET_FASTA in $(find $DATAPATH/maize/target_fastas/ -name GPL*)
do
  GPL=$(echo $TARGET_FASTA | sed 's/.*\(GPL[0-9]*\).*/\1/')
  BLASTRES="data/maize/blast_results/"$GPL"_rice.all.cds.blastn"
  blastn -query $TARGET_FASTA \
         -db data/maize/db/WGS_cds \
         -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore' -num_threads 12 -perc_identity 95 \
         -out $BLASTRES
done

# for barley

mkdir data/barley/blast_results

for TARGET_FASTA in $(find $DATAPATH/barley/target_fastas/ -name GPL*)
do
  GPL=$(echo $TARGET_FASTA | sed 's/.*\(GPL[0-9]*\).*/\1/')
  BLASTRES="data/barley/blast_results/"$GPL"_rice.all.cds.blastn"
  blastn -query $TARGET_FASTA \
         -db data/barley/db/all.cds \
         -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore' -num_threads 12 -perc_identity 95 \
         -out $BLASTRES
done

