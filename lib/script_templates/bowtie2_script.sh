# Run bowtie2 to remove the RNA and/or globin
/usr/local/bin/bowtie2 -p 6 BOWTIE2_OPTIONS -x BOWTIE2_INDEX  --un-conc-gz OUTPUT_DIR/OUTFILE -1 READ1 -2 READ2 > /dev/null 2> OUTPUT_DIR/ERR_OUTFILE && (echo "SAMPLE_NAME" >> OUTPUT_DIR/samples_ok.txt; rm -f READ1 READ2) || (echo "$HOSTNAME SAMPLE_NAME" >> OUTPUT_DIR/samples_ko.txt)

grep "reads; of these:"  OUTPUT_DIR/ERR_OUTFILE > OUTPUT_DIR/RNA_REMOVAL_REPORT
grep "aligned concordantly" OUTPUT_DIR/ERR_OUTFILE >> OUTPUT_DIR/RNA_REMOVAL_REPORT
grep "overall alignment rate" OUTPUT_DIR/ERR_OUTFILE >> OUTPUT_DIR/RNA_REMOVAL_REPORT
rm -f OUTPUT_DIR/ERR_OUTFILE

mv OUTPUT_DIR/OUTFILE.1 OUTPUT_DIR/OUTFILE_noRNA_R1_001.fastq.gz
mv OUTPUT_DIR/OUTFILE.2 OUTPUT_DIR/OUTFILE_noRNA_R2_001.fastq.gz
