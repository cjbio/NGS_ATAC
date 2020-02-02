q#!/bin/sh
#RNA-seq workflow with TOPHAT2 and cufflinks.


# - defining folder structure -

#User must specify the path to the RNA-seq project folder
PROJECTFOLDER=/home/chao/ATAC/ex0079_atac_seq/

#User must specify the name of each sample. i.e. for sample1.fastq.gz and sample2.fastq.gz, write:
SAMPLES="1 4 5 14 15 19 20 24"

#these folders are auto-generated during running the script. Each stage of the processed data will be saved under one these folders. No user intervention needed.
SCRIPTS=$PROJECTFOLDER/scripts/
qFASTQFILES=$PROJECTFOLDER/data/
ALIGNMENTFILES=$PROJECTFOLDER/alignments/
ASSEMBLIESFILES=$PROJECTFOLDER/assemblies/
QUANTIFICATION=$PROJECTFOLDER/quantification
DE=$PROJECTFOLDER/DE/


# - reference genome, index and annotation -

#User to specify the path to these files


# - bowtie2 alignment -




for SAMPLE_ID in $SAMPLES
do
cat > $SCRIPTS/sort_${SAMPLE_ID}.sh <<EOF
sambamba sort -m 10GB -t 32 -o $ALIGNMENTFILES/${SAMPLE_ID}_sorted.bam $ALIGNMENTFILES/${SAMPLE_ID}.bam
EOF

bash $SCRIPTS/sort_${SAMPLE_ID}.sh

done &&

echo "alignment completed"
