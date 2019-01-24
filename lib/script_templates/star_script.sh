#!/bin/bash
echo $HOSTNAME

/usr/local/bin/STAR --quantMode GeneCounts --runThreadN 6 --genomeDir GENOME_DIR --readFilesIn INPUT_FILES1 INPUT_FILES2 --outFileNamePrefix OUTPUT_DIR/SAMPLE_NAME_ --readFilesCommand zcat  STAR_PARAMS || (echo "$HOSTNAME SAMPLE_NAME" >> OUTPUT_DIR/samples_ko.txt)

samtools view -@ 6 -bS OUTPUT_DIR/SAMPLE_NAME_Aligned.out.sam | samtools sort -l 9 -m 2G -@ 6  -o OUTPUT_DIR/SAMPLE_NAME_Aligned.out.bam && (echo "SAMPLE_NAME" >> OUTPUT_DIR/samples_ok.txt) || (echo "$HOSTNAME SAMPLE_NAME" >> OUTPUT_DIR/samples_ko.txt)

#Cleanup the temporary files
rm -f OUTPUT_DIR/SAMPLE_NAME_Aligned.out.sam
rm -f  OUTPUT_DIR/SAMPLE_NAME_Aligned.out.bam.tmp.*.bam

