#!/bin/sh
#RNA-seq workflow with TOPHAT2 and cufflinks.


# - defining folder structure -

#User must specify the path to the RNA-seq project folder
PROJECTFOLDER=/home/chao/ATAC/ex0079_all/

#User must specify the name of each sample. i.e. for sample1.fastq.gz and sample2.fastq.gz, write:
SAMPLES="12"

#these folders are auto-generated during running the script. Each stage of the processed data will be saved under one these folders. No user intervention needed.
SCRIPTS=$PROJECTFOLDER/scripts/
FASTQFILES=$PROJECTFOLDER/data/
ALIGNMENTFILES=$PROJECTFOLDER/alignments/
ASSEMBLIESFILES=$PROJECTFOLDER/assemblies/
QUANTIFICATION=$PROJECTFOLDER/quantification
DE=$PROJECTFOLDER/DE/


# - reference genome, index and annotation -

#User to specify the path to these files


# - bowtie2 alignment -

mkdir $SCRIPTS


for SAMPLE_ID in $SAMPLES
do
cat > $SCRIPTS/bowtie2_${SAMPLE_ID}.sh <<EOF
bowtie2 -p 30 --very-sensitive  -k 10 -X 800  -x /home/chao/index/bowtie2/genome_38 -1 $FASTQFILES/${SAMPLE_ID}_R1.fastq.gz  -2 $FASTQFILES/${SAMPLE_ID}_R2.fastq.gz -S ${SAMPLE_ID}.sam

EOF

bash $SCRIPTS/bowtie2_${SAMPLE_ID}.sh

done &&

echo "alignment completed"
