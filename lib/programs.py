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
slurm_partiton = 'aRNA-Seq'

def submit_job_super(pname, path, wt, nproc, q, ns, bsub_suffix, nstar, timestamp):
    folder = path.split("/")[-1]
    if pname.startswith("htseq"):
        job_name = "hts" + (pname.split("-")[1][0]).upper()
    elif pname.startswith("star-fusion"):
        job_name = "sFusion"
    else:
        job_name = pname[0:3]
    print "> Submitting " + pname + " job(s) to cluster..."
    job_ids = list()
    logs = list()
    for i in range(nstar):
        logname = path + "/logs/" + pname + "_" + str(i) + ".log"
        logs.append(logname)
        #if os.path.exists(logname):
        #    os.remove(logname)
        #job_ids.append(str(manager.submit_job(wt, str(nproc), slurm_partiton, logname, folder + "_" + job_name + "_" + str(i),
        job_ids.append(str(manager.submit_job(wt, "1", slurm_partiton, logname, folder + "_" + job_name + "_" + str(i),
                              path + "/results_" + pname + "/script_" + str(i) + ".sh", 1, path, bsub_suffix)))
    job_ids = "|".join(job_ids)
    logs = "|".join(logs)
    out = open(path + "/temp/pids.txt", 'a')
    #print "Path = " + path + "/temp/pids.txt"
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


def compute_mean_std(path_base, folder, samples, output, nproc, wt, q):
    ########################################################################
    ## Computes mean fragment length for kallisto single-end reads
    ########################################################################
    f  = list()
    s  = list()
    for i,j in sorted(samples.iteritems()):
        f.append("'"+j[0]+"'")
        s.append("'"+i+"'")
    s  = "samples = ["+",".join(s)+"]"
    f  = "files   = ["+",".join(f)+"]"
    o  = "output  = '"+output+"'"
    a1 = 'out = open(output, "w")\nprint >> out, "sample N M V"\nfor isamp in range(len(samples)):\n    f  = open(files[isamp],"r")\n    k  = 0'
    a2 = '    kt = 0\n    r  = list()\n    for i in f:\n        if (k % 4)==1:\n            l = float(len(i.rstrip()))\n            r.append(l)\n            kt += l\n        k += 1'
    a3 = '    f.close()\n    N = len(r)\n    M = float(kt)/N\n    t = 0\n    for i in r:\n        t += ((i-M)*(i-M))\n    S = (t/N)**0.5\n    if S < 0.001:\n        S = 0.001\n    print >> out, samples[isamp]+" "+str(int(N))+" "+str(round(M,3))+" "+str(round(S,3))'
    a4 = 'out.close()\n'
    out = open(path_base + folder + "/results_kallisto/script_stats.py",'w')
    print >> out, s
    print >> out, f
    print >> out, o
    print >> out, a1
    print >> out, a2
    print >> out, a3
    print >> out, a4
    out.close()
    print "> Submitting Kallisto stat reads job to cluster..."
    job_id = manager.submit_job(wt, str(min(int(nproc),len(samples))), q, path_base+folder+"/results_kallisto/log_cluster_stats.txt", "Jmean", "python "+path_base+folder+"/results_kallisto/script_stats.py", 0, path_base+folder, "")
    print "> Kallisto stats job ID: "+ job_id
    return job_id, path_base+folder+"/results_kallisto/log_cluster_stats.txt"


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
    #configuration, nchild, bsub_suffix = manager.get_bsub_arg(configuration, len(samples))
    bsub_suffix = manager.get_bsub_arg(configuration, len(samples))
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
    return submit_job_super("trimgalore", path_base + folder, wt, "1", q, len(samples), bsub_suffix, len(samples), timestamp)


def fastqc(timestamp, path_base, folder, samples, nproc, wt, q, tg):
    ########################################################################
    ## FastQC analysis
    ########################################################################
    print "## QC: FastQC"
    print "> Quality control with fastQC..."
    output = "results_fastqc"
    secure_mkdir(path_base + folder, "results_fastqc")
    output_folder = path_base + folder + "/results_fastqc"
    print "> Writing jobs for fastqc analysis..."
    bsub_suffix = manager.get_bsub_arg(nproc, len(samples))
    commands = list()
    ksamp = sortbysize(samples)
    for sample in ksamp:
        files = samples[sample]
        if not tg:
            if len(files) == 4:
                fnames = files[0] + " " + files[1]
            else:
                fnames = files[0]
        else:
            g = path_base + folder + "/results_trimgalore/"
            suf = ""
            if not files[0].split("/")[-1].endswith(".gz"):
                suf = ".gz"
            if len(files) == 4:
                fnames = g + files[0].split("/")[-1] + suf + " " + g + files[1].split("/")[-1] + suf
            else:
                fnames = g + files[0].split("/")[-1] + suf
        call = config.path_fastqc + " --extract -q -o " + output_folder + " " + fnames
        commands.append(call + sample_checker.replace("#FOLDER", output_folder).replace("#SAMPLE", sample))
    create_scripts(len(samples), commands, path_base, folder, output)
    return submit_job_super("fastqc", path_base + folder, wt, str(nproc), q, len(samples), bsub_suffix, len(samples), timestamp)


def kallisto(timestamp, path_base, folder, samples, path_index, bootstrap, nproc, wt, q, tg):
    output = "results_kallisto"
    secure_mkdir(path_base + folder, "results_kallisto")
    print "## RNAseq pseudoalignment with Kallisto"
    # Estimate counts in single-end datasss
    if len(samples[samples.keys()[0]]) == 2:
        print "> Estimating average and STD of fragment lengh required by Kalisto on single-read data..."
        outputT = path_base + folder + "/" + output + "/stats.txt"
        tid,log = compute_mean_std(path_base, folder, samples, outputT, "1", wt, q)
        vcrparser.job_wait(log, 10)
        f = open(outputT,'r')
        i = f.readline()
        stats = dict()
        for i in f:
            i = i.strip("\n").split(" ")
            stats[i[0]] = [i[2],i[3]]
        f.close()
    print "> Writing jobs for Kallisto pseudoalignment"
    bsub_suffix = manager.get_bsub_arg(nproc, len(samples))
    commands = list()
    ksamp = sortbysize(samples)
    for sample in ksamp:
        files = samples[sample]
        if not tg:
            if len(files) == 4:
                args = ""
                fnames = files[0]+" "+files[1]
            else:
                args = " --single -l mean -s var".replace("mean", stats[sample][0]).replace("var", stats[sample][1])
                fnames = files[0]
        else:
            g = path_base + folder + "/results_trimgalore/"
            suf = ""
            if not files[0].split("/")[-1].endswith(".gz"):
                suf = ".gz"
            if len(files) == 4:
                args = ""
                fnames = g + files[0].split("/")[-1] + suf + " " + g + files[1].split("/")[-1] + suf
            else:
                args = " --single -l mean -s var".replace("mean", stats[sample][0]).replace("var", stats[sample][1])
                fnames = g + files[0].split("/")[-1] + suf
        cmd = config.path_kallisto+" quant -b " + bootstrap + " -i " + path_index + " -o " + path_base+folder + "/results_kallisto/" + sample + args + " " + fnames
        commands.append(cmd + sample_checker.replace("#FOLDER", path_base + folder + "/" + output).replace("#SAMPLE", sample))
    create_scripts(len(samples), commands, path_base, folder, output)
    return submit_job_super("kallisto", path_base + folder, wt, str(nproc), q, len(samples), bsub_suffix, len(samples), timestamp)


def star(timestamp, path_base, folder, samples, nproc, wt, q, path_genome, star_params, tg):
    output = "results_star"
    secure_mkdir(path_base + folder, output)
    print "## RNAseq alignment with STAR..."
    print "> Writing jobs for STAR alignment..."
    bsub_suffix = manager.get_bsub_arg(nproc, len(samples))
    commands = list()
    ksamp = sortbysize(samples)
    for sample in ksamp:
        files = samples[sample]
        gg = ""
        fn1 = ""
        fn2 = ""
        if not tg:
            if len(files) == 2:
                fn1 = files[0]
            else:
                fn1 = files[0]
                fn2 = files[1]
            if files[0].endswith(".fastq.gz"):
                gg = " --readFilesCommand zcat"
        else:
            gg = " --readFilesCommand zcat"
            g = path_base + folder + "/results_trimgalore/"
            suf = ""
            if not files[0].split("/")[-1].endswith(".gz"):
                suf = ".gz"
            if len(files) == 2:
                fn1 = g + files[0].split("/")[-1] + suf
            else:
                fn1 = g + files[0].split("/")[-1] + suf 
                fn2 = g + files[1].split("/")[-1] + suf
        #command = "CT=$(lscpu --online --parse | wc -l)\n"
        #command = command + "CPUs=$(expr $CT - 6)\n"
        #command = command + config.path_star + " --quantMode GeneCounts --runThreadN $CPUs --genomeDir " + path_genome
        #command = command + " --readFilesIn " + fn + " --outFileNamePrefix " + path_base + folder + "/results_star/" + sample + "_" + gg
        #if len(star_params) > 0:
        #    command = command + star_params
        #commands.append(command + sample_checker.replace("#FOLDER", path_base + folder + "/results_star").replace("#SAMPLE", sample))
        with open('/share/code/lib/script_templates/star_script.sh','r') as script_template:
            command = script_template.read()
            script_template.close()
        command = command.replace("OUTPUT_DIR", path_base+folder+"/results_star/")
        command = command.replace("INPUT_FILES1", fn1)
        command = command.replace("INPUT_FILES2", fn2)
        command = command.replace("GENOME_DIR", path_genome)
        command = command.replace("SAMPLE_NAME", sample)
        command = command.replace("STAR_PARAMS", star_params)
        commands.append(command)
    create_scripts(len(samples), commands, path_base, folder, output)
    return submit_job_super("star", path_base + folder, wt, str(nproc), q, len(samples), bsub_suffix, len(samples), timestamp)


def htseq(timestamp, path_base, folder, samples, path_annotation, nproc, wt, q, mode, strand, countmode):
    output = "results_htseq-" + mode
    secure_mkdir(path_base + folder, output)
    print "## HTseq-count"
    print "> Writing jobs for HTseq-count " + mode + " analysis..."
    bsub_suffix = manager.get_bsub_arg(nproc, len(samples))
    commands = list()
    ksamp = sortbysize(samples)
    proc_files = os.listdir(path_base + folder + "/results_star/")
    for sample in ksamp:
        in_file = path_base + folder + "/results_star/" + sample + "_Aligned.out.bam"
        if sample + "_Aligned.out.bam" in proc_files:
            outputf= path_base + folder + "/results_htseq-" + mode + "/" + sample + ".tab"
            if mode == "gene":
                ld1 = config.path_htseq + strand + " -f bam -m " + countmode  + " -q " + in_file + " " + path_annotation
            else:
                ld1 = config.path_htseq + strand + " -f bam -m " + countmode  + " -i exon_id -q " + in_file + " " + path_annotation
            call = ld1 + " > " + outputf
            commands.append(call  + sample_checker.replace("#FOLDER", path_base + folder + "/" + output).replace("#SAMPLE", sample))
        else:
            print "Warning: [HTseq-" + mode + "] STAR output file not found -> " + in_file
    create_scripts(len(samples), commands, path_base, folder, output)
    return submit_job_super("htseq-" + mode, path_base + folder, wt, str(nproc), q, len(samples), bsub_suffix, len(samples), timestamp)


def sam2sortbam(timestamp, path_base, folder, samples, nproc, wt, q):
    output = "results_sam2sortbam"
    secure_mkdir(path_base + folder, output)
    print "## SAM2SORTEDBAM"
    print "> Writing jobs for SAM2SORTEDBAM..."
    bsub_suffix = manager.get_bsub_arg(nproc, len(samples))
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
    return submit_job_super("sam2sortbam", path_base + folder, wt, str(nproc), q, len(samples), bsub_suffix, len(samples), timestamp)
