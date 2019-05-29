# -*- coding: utf-8 -*-
import os
import config
import vcrparser


def write_log(message):
    fh = open("/share/log.txt", 'a')
    print >> fh, message
    fh.close()

# Conditional library load according to 'config.mode'
if config.mode == "LSF":
    import sys_LSF as manager
elif config.mode == "LOCAL":
    import sys_single as manager
else:
    import sys_OTHER as manager

sample_checker = ' && (echo "#SAMPLE" >> #FOLDER/samples_ok.txt) || (echo "$HOSTNAME #SAMPLE" >> #FOLDER/samples_ko.txt)'
slurm_partiton = os.getenv('SLURM_PARTITION', 'aRNA-Seq')

def submit_job_super(pname, ncpu, path, wt, q, nstar, timestamp):
    folder = path.split("/")[-1]
    if pname.startswith("htseq"):
        job_name = "hts" + (pname.split("-")[1][0]).upper()
    else:
        job_name = pname[0:3]
    print "> Submitting " + pname + " job(s) to cluster..."
    job_ids = list()
    logs = list()
    for i in range(nstar):
        logname = path + "/logs/" + pname + "_" + str(i) + ".log"
        logs.append(logname)
        job_ids.append(str(manager.submit_job(wt, ncpu, slurm_partiton, logname, folder + "_" + job_name + "_" + str(i),
                              path + "/results_" + pname + "/script_" + str(i) + ".sh", 1, path)))
    job_ids = "|".join(job_ids)
    logs = "|".join(logs)
    out = open(path + "/temp/pids.txt", 'a')
    print >> out, pname + "\t" + job_ids + "\t" + logs
    out.close()
    return job_ids, logs


def secure_mkdir(path,folder):
    L = os.listdir(path)
    if not (folder in L):
        os.mkdir(path + "/" + folder)
        print "> '" + folder + "' folder created"
    else:
        print "> - '" + folder + "' folder already exists"
    return 1


def create_scripts(sample_count, commands, path_base, folder, output):
    out = list()
    for i in range(sample_count):
        out.append(open(path_base + folder + "/" + output + "/script_" + str(i) + ".sh",'w'))
        print >> out[i], "#!/bin/bash\necho $HOSTNAME\n"
    for k in range(len(commands)):
        print >> out[k % sample_count], commands[k]
    for i in out:
        i.close()


###################################
#
# Inputs
# samples: a dictionary object
#
# Returns
# keys: a list containing sorted items from the dictionary
###################################
def sortbysize(samples):
    sizes = dict()
    for sample, files in sorted(samples.iteritems()):
        if len(files) == 2:
            r = files[1]
        else:
            r = files[2] + files[3]
        if not sizes.has_key(r):
            sizes[r] = [sample]
        else:
            sizes[r].append(sample)
    keys = list()
    for i,j in sorted(sizes.iteritems()):
        for k in j:
            keys.append(k)
    return keys


def rename_tg_output(sample, files, path):
    g = path + "/results_trimgalore/"
    cmds = ['sleep 10']
    if len(files) == 4:
        for i in range(2):
            output = g + files[i].split("/")[-1].replace(".fastq.gz","").replace(".fastq","") + "_val_" + str(i+1) +".fq.gz"
            cmds.append("mv '" + output + "' '" + output.replace("_val_" + str(i+1) +".fq.gz", ".fastq.gz") + "' || true")
    else:
        output = g + files[0].split("/")[-1].replace(".fastq.gz","").replace(".fastq","") + "_trimmed.fq.gz"
        cmds.append("mv '" + output + "' '" + output.replace("_trimmed.fq.gz", ".fastq.gz") + "' || true")
    return "\n".join(cmds)


def trimgalore(timestamp, path_base, folder, samples, configuration, wt, q, extra_args):
    ########################################################################
    ## Trimgalore analysis
    #  samples is a dictionary keyed to the sample ID. It returns a list with four items:
    #   1. First read file
    #   2. Second read file
    #   3. Size of first read file in bytes
    #   4. Size of second read file in bytes
    ########################################################################
    print "## Trim-galore: Quality and adapter trimming"
    print "> Quality and adapter trimming with Trim Galore..."
    output = "results_trimgalore"
    secure_mkdir(path_base + folder, "results_trimgalore")
    output_folder = path_base + folder + "/results_trimgalore"
    print "> Writing jobs for TrimGalore analysis..."
    commands = list()
    ksamp = sortbysize(samples)
    for sample in ksamp:
        files = samples[sample]
        if len(files) == 4:
            args = extra_args + " --paired"
            fnames = files[0] + " " + files[1]
        else:
            args = extra_args
            fnames = files[0]
        if (args != "") and (not args.startswith(" ")):
            args = " " + args
        call = config.path_trimgalore.strip() + " " + args + " --gzip --path_to_cutadapt " + config.path_cutadapt + " -o " + output_folder + " " + fnames
        #call = config.path_trimgalore.strip() + " " + args + "  --path_to_cutadapt " + config.path_cutadapt + " -o " + output_folder + " " + fnames  # NO GZIP fcf
        call = call + sample_checker.replace("#FOLDER", output_folder).replace("#SAMPLE", sample) + "\n" + rename_tg_output(sample, files, path_base + folder)
        commands.append(call)
    create_scripts(len(samples), commands, path_base, folder, output)
    return submit_job_super("trimgalore", 2, path_base + folder, wt, q, len(samples), timestamp)

def bowtie2(timestamp, path_base, folder, samples, bowtie2_options, bowtie2_index, wt, q):
    output = "results_bowtie2"
    secure_mkdir(path_base + folder, output)
    print "## RNA removal with Bowtie2..."
    print "> Writing jobs for Bowtie2 filtering..."
    commands = list()
    ksamp = sortbysize(samples)
    for sample in ksamp:
        files = samples[sample]
        fn1 = ""
        fn2 = ""
        g = path_base + folder + "/results_trimgalore/"
        suf = ""
        if not files[0].split("/")[-1].endswith(".gz"):
            suf = ".gz"
        fn1 = files[0].split("/")[-1] + suf
        base_name = fn1[:len(fn1)-16] 
        fn1 = g + fn1
        if len(files) > 2:
            fn2 = g + files[1].split("/")[-1] + suf
        #with open('/share/code/lib/script_templates/bowtie2_script.sh','r') as script_template:
        with open(config.path_code + '/lib/script_templates/bowtie2_script.sh','r') as script_template:
            command = script_template.read()
            script_template.close()
        command = command.replace("OUTPUT_DIR", path_base+folder+"/results_bowtie2/")
        command = command.replace("READ1", fn1)
        command = command.replace("READ2", fn2)
        command = command.replace("BOWTIE2_INDEX", bowtie2_index)
        command = command.replace("SAMPLE_NAME", sample)
        command = command.replace("BOWTIE2_OPTIONS", bowtie2_options)
        command = command.replace("OUTFILE", base_name)
        command = command.replace("ERR_OUTFILE", base_name+"_err")
        command = command.replace("RNA_REMOVAL_REPORT", sample+".RNA_removal.report.txt")
        commands.append(command)
    create_scripts(len(samples), commands, path_base, folder, output)
    return submit_job_super("bowtie2", 5, path_base + folder, wt, q, len(samples), timestamp)
 
def fastqc(timestamp, path_base, folder, samples, wt, q, bt2):
    ########################################################################
    ## FastQC analysis
    ########################################################################
    print "## QC: FastQC"
    print "> Quality control with fastQC..."
    output = "results_fastqc"
    secure_mkdir(path_base + folder, "results_fastqc")
    output_folder = path_base + folder + "/results_fastqc"
    print "> Writing jobs for fastqc analysis..."
    commands = list()
    ksamp = sortbysize(samples)
    for sample in ksamp:
        files = samples[sample]
        suf = ""
        if not files[0].split("/")[-1].endswith(".gz"):
            suf = ".gz"
        if not bt2:
            g = path_base2+ folder + "/results_trimgalore/"
            fn1 = files[0].split("/")[-1] + suf
            fn2 = files[1].split("/")[-1] + suf
        else:
            g = path_base + folder + "/results_bowtie2/"
            base_name = files[0].split("/")[-1]
            fn1 = base_name[:len(base_name)-16]+"_noRNA_R1_001.fastq.gz"
            fn2= base_name[:len(base_name)-16]+"_noRNA_R2_001.fastq.gz"
        if len(files) == 4:
            fnames = g + fn1 + " " + g + fn2
        else:
            fnames = g + fn1
        call = config.path_fastqc + " --extract -q -o " + output_folder + " " + fnames
        commands.append(call + sample_checker.replace("#FOLDER", output_folder).replace("#SAMPLE", sample))
    create_scripts(len(samples), commands, path_base, folder, output)
    return submit_job_super("fastqc", 5, path_base + folder, wt, q, len(samples), timestamp)

def star(timestamp, path_base, folder, samples, wt, q, path_genome, star_params, bt2):
    output = "results_star"
    secure_mkdir(path_base + folder, output)
    print "## RNAseq alignment with STAR..."
    print "> Writing jobs for STAR alignment..."
    commands = list()
    ksamp = sortbysize(samples)
    for sample in ksamp:
        files = samples[sample]
        fn1 = ""
        fn2 = ""
        suf = ""
        bowtie2_suffix = ""
        if not files[0].split("/")[-1].endswith(".gz"):
            suf = ".gz"
        if bt2:
            bowtie2_suffix = "_noRNA"
            g = path_base + folder + "/results_bowtie2/"
            base_name = files[0].split("/")[-1]
            fn1 = g + base_name[:len(base_name)-16]+"_noRNA_R1_001.fastq.gz"
            fn2 = g + base_name[:len(base_name)-16]+"_noRNA_R2_001.fastq.gz"
        else:
            g = path_base + folder + "/results_trimgalore/"
            fn1 = g + files[0].split("/")[-1] + suf
            fn2 = g + files[1].split("/")[-1] + suf
        with open(config.path_code + '/lib/script_templates/star_script.sh','r') as script_template:
            command = script_template.read()
            script_template.close()
        command = command.replace("OUTPUT_DIR", path_base+folder+"/results_star/")
        command = command.replace("INPUT_FILES1", fn1)
        command = command.replace("INPUT_FILES2", fn2)
        command = command.replace("GENOME_DIR", path_genome)
        command = command.replace("SAMPLE_NAME", sample)
        command = command.replace("STAR_PARAMS", star_params)
        command = command.replace("BOWTIE2_SUFFIX", bowtie2_suffix)
        commands.append(command)
    create_scripts(len(samples), commands, path_base, folder, output)
    return submit_job_super("star", 10, path_base + folder, wt, q, len(samples), timestamp)

def htseq(timestamp, path_base, folder, samples, path_annotation, wt, q, mode, strand, countmode, bt2):
    output_dir = "results_htseq-" + mode
    secure_mkdir(path_base + folder, output_dir)
    print "## HTseq-count"
    print "> Writing jobs for HTseq-count " + mode + " analysis..."
    commands = list()
    ksamp = sortbysize(samples)
    proc_files = os.listdir(path_base + folder + "/results_star/")
    for sample in ksamp:
        if bt2:
            in_file = path_base + folder + "/results_star/" + sample + "_noRNA_Aligned.out.bam"
        else:
            in_file = path_base + folder + "/results_star/" + sample + "_Aligned.out.bam"
        if in_file.split("/")[-1] in proc_files:
            outputf = path_base + folder + "/results_htseq-" + mode + "/" + sample + ".tab"
            if mode == "gene":
                ld1 = config.path_htseq + strand + " -f bam -m " + countmode  + " -q " + in_file + " " + path_annotation
            else:
                ld1 = config.path_htseq + strand + " -f bam -m " + countmode  + " -i exon_id -q " + in_file + " " + path_annotation
            call = ld1 + " > " + outputf
            commands.append(call  + sample_checker.replace("#FOLDER", path_base + folder + "/" + output_dir).replace("#SAMPLE", sample))
        else:
            print "Warning: [HTseq-" + mode + "] STAR output file not found -> " + in_file
    create_scripts(len(samples), commands, path_base, folder, output_dir)
    return submit_job_super("htseq-" + mode, 5, path_base + folder, wt, q, len(samples), timestamp)


def sam2sortbam(timestamp, path_base, folder, samples, wt, q):
    output = "results_sam2sortbam"
    secure_mkdir(path_base + folder, output)
    print "## SAM2SORTEDBAM"
    print "> Writing jobs for SAM2SORTEDBAM..."
    commands = list()
    ksamp = sortbysize(samples)
    proc_files = os.listdir(path_base + folder + "/results_star/")
    for sample in ksamp:
        in_file = path_base + folder + "/results_star/" + sample + "_Aligned.out.sam"
        if sample + "_Aligned.out.sam" in proc_files:
            out_file = path_base + folder + "/results_sam2sortbam/" + sample + ".sorted.bam"
            com = "java -jar " + config.path_picard + "/AddOrReplaceReadGroups.jar I=" + in_file + " O=" + out_file +" SO=coordinate RGID=id RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample 2> " + out_file + ".log"
            commands.append(com + sample_checker.replace("#FOLDER", path_base + folder + "/results_sam2sortbam").replace("#SAMPLE", sample))
        else:
            print "Warning: [SAM2SORTEDBAM] STAR output file not found -> " + in_file
    create_scripts(len(samples), commands, path_base, folder, output)
    return submit_job_super("sam2sortbam", 5, path_base + folder, wt, q, len(samples), timestamp)
