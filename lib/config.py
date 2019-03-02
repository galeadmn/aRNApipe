# -*- coding: utf-8 -*-
import sys
import os

# LIBRARY USED TO SUBMIT JOBS:
# - 'LSF' FOR IBM LSF WORKLOAD MANAGER (it uses 'sys_LSF.py')
# - 'LOCAL' FOR SEQUENTIAL RUN ON SINGLE MACHINE (it uses 'sys_single.py')
# - 'OTHER' FOR LIBRARIES ADAPTED TO OTHER WORKLOAD MANAGERS (it uses 'sys_OTHER.py')
mode = "OTHER"

# PATH TO THE FOLDER "genomes_processed" WHERE THE DIFFERENT GENOME BUILDS ARE STORED
path_db = "/share/"
path_code = "/share/code/"

# FULL PATHS TO BINARIES USED BY aRNApipe (users must change these values to match
# the current locations of the binaries used by aRNApipe in their system).
path_trimgalore = "/usr/local/bin/trim_galore"
path_bowtie2  = "/usr/local/bin/bowtie2"
path_fastqc   = "/usr/local/FastQC/fastqc"
path_star     = "/usr/local/bin/STAR"
path_htseq    = "/usr/bin/htseq-count"
path_samtools = "/usr/local/bin/samtools"
path_cutadapt = "/usr/local/bin/cutadapt"

# STAR options (users can add their own options):
# The keys of this dict are used in the project config files to use the
# referenced STAR arguments within the corresponding dictionary values
star_options = {"default": "",
                "encode":  "--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000"}

# ENVIRONMENT VARIABLES:
# The following system environment variables are changed to add or overwrite
# their current values.
environment = {"JAVA_HOME":      ["/usr/java/jdk1.8.0_60/","add"],
                "PYTHONPATH":     ["/usr/lib64/python2.7/site-packages","overwrite"],
#               "PATH":           ["/gpfs/gpfs1/software/Python-2.7.2/bin","add"],
#               "PATH":           ["/gpfs/gpfs1/software/bedtools2-2.20.0/bin","add"],
#               "PATH":           ["/gpfs/gpfs1/software/samtools-1.2/bin","add"],
#               "LD_LIBRARY_PATH":["/gpfs/gpfs1/software/gcc-4.8.2/usr/lib64","add"],
               "PERL5LIB"       :["/gpfs/gpfs1/software/perl-modules/lib/perl5/5.10.1:/gpfs/gpfs1/software/perl-modules/lib/perl5/5.10.1/lib64/perl5","add"]}

# ANNOTATIONS AND FULL PATH TO THE PIPELINE BASE DIRECTORY (do not change)
path_genome = path_db + "/genomes_processed/#LABEL/STAR_genome"
path_annotation = path_db + "/genomes_processed/#LABEL/genesets.gtf"
path_fasta = path_db + "/genomes_processed/#LABEL/genome.fa"
annots = [path_db + "/genomes_processed/#LABEL/genesets.refFlat",
          path_db + "/genomes_processed/#LABEL/refFlats/protein_coding.refFlat",
          path_db + "/genomes_processed/#LABEL/refFlats/rRNA.refFlat"]
nannots = ["general","protein_coding","ribosomal"]
