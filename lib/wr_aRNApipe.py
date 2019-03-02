# -*- coding: utf-8 -*-
import os
import time
import optparse
import config
import vcrparser
import programs
import shutil

if config.mode == "LSF":
    import sys_LSF as manager
elif config.mode == "LOCAL":
    import sys_single as manager
else:
    import sys_OTHER as manager

################################################################################
## HELP OPTIONS AND DOCUMENTATION
desc = "aRNApipe: RNA-seq framework"
parser = optparse.OptionParser(description=desc)
parser.add_option("-f", "--project_folder", dest="folder", default="",  help="")
parser.add_option("-m", "--write", dest = "m", default="0", help="")
parser.add_option("-b", "--path_base", dest="path_base", default="", help="")
parser.add_option("-t", "--timestamp", dest="timestamp", default="", help="")
################################################################################

(opt, args) = parser.parse_args()
timestamp = opt.timestamp
print "###########################################"
print "### aRNApipe: RNA-seq framework    #########"
print "###########################################"
print "> Analysis started: " + timestamp
print "> Parsing options..."

## PATH BASE AND PROJECT FOLDER CHECK
[path_base, folder] = vcrparser.check_args(opt.path_base, opt.folder)

## PARSE CONFIGURATION FILE
opt.samples = path_base + folder + "/samples.list"
opt.config = path_base + folder + "/config.txt"
[opt.config, var] = vcrparser.config_file(opt.config, path_base, folder, config)
shutil.copy(opt.samples, opt.samples.replace("samples.list", "logs/" + timestamp + "_samples.list"))
shutil.copy(opt.config, opt.config.replace("config.txt","logs/" + timestamp + "_config.txt"))

## CHECK GENOME BUILD
config.path_genome = config.path_genome.replace("#LABEL", var["genome_build"])
if not os.path.exists(config.path_genome):
    exit("path_genome not found. Genome build " + var["genome_build"] + " missing or incomplete.")
config.path_annotation = config.path_annotation.replace("#LABEL", var["genome_build"])
print config.path_annotation
if not os.path.exists(config.path_annotation):
    exit("path_annotation not found. Genome build ")
for i in range(3):
    config.annots[i] = config.annots[i].replace("#LABEL", var["genome_build"])
    if not os.path.exists(config.annots[i]):
        print config.annots[i]  
        exit("annots not found. ")

var["htseq-gene"] = var["htseq-gene"].strip()
var["htseq-exon"] = var["htseq-exon"].strip()
var["htseq-gene-mode"] = var["htseq-gene-mode"].strip()
var["htseq-exon-mode"]

## VERBOSE OPTIONS
print "> GLOBAL VARIABLES:"
print "  - Project folder:  " + folder
print "  - Base path:       " + path_base
print "  - Config file:     " + opt.config
print "  - Genome build:    " + var["genome_build"]
print "  - Bowtie2 index:   " + var["bowtie2_index"]
print "  - Run mode:        " + opt.m
print "  - Stranded RNA:    " + var["strandedness"]
print "> CLUSTER VARIABLES:"
print "  - Queue:           " + var["q"]
print "  - Walltime:        " + var["wt"]
print "> ANALYSIS:"
print "  - TrimGalore:      " + var["trimgalore"]
print "  - Bowtie2:         " + var["bowtie2"]
print "  - FastQC:          " + var["fastqc"]
print "  - STAR:            " + var["star"]
print "  - HTseq (gene):    " + var["htseq-gene"]
print "  - HTseq (exon):    " + var["htseq-exon"]
if var.has_key('sam2sortbam'):
    print "  - sam2sortbam:     " + var["sam2sortbam"]
print "> INDIVIDUAL ANALYSIS SETTINGS:"
print "  - TrimGalore args: " + var["trimgal_args"]
print "  - STAR arguments:  " + var["star_args"]
print "  - STAR 2-pass:     " + var["star2pass"]
print "  - HTseqGene mode:  " + var["htseq-gene-mode"]
print "  - HTseqExon mode:  " + var["htseq-exon-mode"]

samples = vcrparser.get_samples(path_base, folder, opt.samples)

##########################################################
## Creates 'temp' folder for temporary managing files
##########################################################
n = os.listdir(opt.path_base + "/" + opt.folder)
if not ("temp" in n):
    os.mkdir(path_base + "/" + folder + "/temp")
else:
    out = open(path_base + "/" + folder + "/temp/pids.txt", 'w')
    out.close()
    out = open(path_base + "/" + folder + "/temp/pids_scheduler.txt", 'w')
    out.close()

##########################################################
## Strand-specific protocol and HTseq mode
##########################################################
var["strandedness"] = var["strandedness"].strip()
if var["strandedness"] == "yes":
    var["strandedness"] = " --stranded=yes"
elif var["strandedness"] == "reverse":
    var["strandedness"] = " --stranded=reverse"
elif var["strandedness"] == "no":
    var["strandedness"] = " --stranded=no"
else:
    exit("Error: Strandedness value not correct (yes/no/reverse).")
if (int(var["htseq-gene"].split("/")[0]) > 0) or (int(var["htseq-exon"].split("/")[0]) > 0):
    var["htseq-gene-mode"] == "union"
    var["htseq-exon-mode"] = var["htseq-exon-mode"].strip()

    if (not (var["htseq-gene-mode"] in ["union", "intersection-strict", "intersection-nonempty"])) or (not (var["htseq-exon-mode"] in ["union", "intersection-strict", "intersection-nonempty"])):
        exit("Error: HTseq mode value not correct (union/intersection-strict/intersection-nonempty).")

##########################################################
## Check STAR options
##########################################################
if int(var["star"].split("/")[0]) > 0:
    var["star_args"] = var["star_args"].strip()
    if var["star_args"] == "own":
            var["star_args"] = var["star_args_own"].strip()
    else:
        if config.star_options.has_key(var["star_args"]):
            var["star_args"] = config.star_options[var["star_args"]].strip()
        else:
            exit("star_args key not found in 'config.py' file.")
    if var["star2pass"] == "yes":
        var["star2pass"] = "--twopassMode Basic"
    else:
        var["star2pass"] = ""
    args = dict()
    for i in ["star2pass", "star_args"]:
        if var[i] != "":
            j = var[i].split("--")
            for k in j:
                if len(k) > 0:
                    k = k.split(' ')
                    args['--' + k[0]] = (' '.join(k[1:])).replace('  ', '')
    star_params = ""
    if len(args) > 0:
        for i, j in args.iteritems():
            star_params = star_params + " " + i + " " + j
##########################################################
## Starts analysis
##########################################################
procs = list()
bt2 = False
## Pre-alignment analysis
if int(var["trimgalore"].split("/")[0]) > 0:
    samples_v, stats = vcrparser.check_samples(samples, path_base, folder, "trimgalore", opt.m)
    if len(samples_v) > 0:
        print "About to call trimgalore"
        job_id_tg, logs_tg = programs.trimgalore(timestamp, path_base, folder, samples_v, var["trimgalore"], var["wt"], var["q"], var["trimgal_args"])
        w = vcrparser.job_wait(job_id_tg, 20)
        procs.append(job_id_tg)
if int(var["bowtie2"].split("/")[0]) > 0:
    bt2 = True
    samples_v, stats = vcrparser.check_samples(samples, path_base, folder, "bowtie2", opt.m)
    if len(samples_v) > 0:
        print "About to call bowtie2"
        job_id_btw, logs_qc = programs.bowtie2(timestamp, path_base, folder, samples_v, var["bowtie2_args"], var["bowtie2_index"], var["wt"], var["q"])
        w = vcrparser.job_wait(job_id_btw, 20)
if int(var["fastqc"].split("/")[0]) > 0:
    samples_v, stats = vcrparser.check_samples(samples, path_base, folder, "fastqc", opt.m)
    if len(samples_v) > 0:
        print "About to call fastqc"
        job_id_qc, logs_qc = programs.fastqc(timestamp, path_base, folder, samples_v, var["wt"], var["q"], bt2)
        procs.append(logs_qc)
if int(var["star"].split("/")[0]) > 0:
    samples_v, stats = vcrparser.check_samples(samples, path_base, folder, "star", opt.m)
    if len(samples_v) > 0:
        job_id_star, logs_star = programs.star(timestamp, path_base, folder, samples_v, var["wt"], var["q"], config.path_genome, star_params, bt2)
        w = vcrparser.job_wait(job_id_star, 20)
        procs.append(job_id_star)
if int(var["htseq-gene"].split("/")[0]) > 0:
    samples_v, stats = vcrparser.check_samples(samples, path_base, folder, "htseq-gene", opt.m)
    if len(samples_v) > 0:
        job_id_htseqG, logs_htseqG = programs.htseq(timestamp, path_base, folder, samples_v, config.path_annotation, var["wt"], var["q"], "gene", var["strandedness"], var["htseq-gene-mode"], bt2)
        procs.append(logs_htseqG)
if int(var["htseq-exon"].split("/")[0]) > 0:
    samples_v, stats = vcrparser.check_samples(samples, path_base, folder, "htseq-exon", opt.m)
    if len(samples_v) > 0:
        job_id_htseqE, logs_htseqE = programs.htseq(timestamp, path_base, folder, samples_v, config.path_annotation, var["htseq-exon"], var["wt"], var["q"], "exon", var["strandedness"],var["htseq-exon-mode"])
        procs.append(logs_htseqE)
if len(procs) > 0:
    for proc in procs:
        w = vcrparser.job_wait(proc, 10)

timestamp = time.strftime("%y%m%d_%H%M%S")
print "> Analysis finished: " + timestamp
