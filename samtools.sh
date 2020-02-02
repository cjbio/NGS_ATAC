#!/bin/sh
#RNA-seq workflow with TOPHAT2 and cufflinks.


# - defining folder structure -

#User must specify the path to the RNA-seq project folder
PROJECTFOLDER=/home/chao/ATAC/ex0079_all/

#User must specify the name of each sample. i.e. for sample1.fastq.gz and sample2.fastq.gz, write:
SAMPLES="2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 64 65"

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
cat > $SCRIPTS/samtools_${SAMPLE_ID}.sh <<EOF
newsamtools view -bhS -O BAM -@ 30 $ALIGNMENTFILES/${SAMPLE_ID}.sam > $ALIGNMENTFILES/${SAMPLE_ID}.bam &&
rm $ALIGNMENTFILES/${SAMPLE_ID}.sam
EOF

bash $SCRIPTS/samtools_${SAMPLE_ID}.sh

done &&

echo "alignment completed"
