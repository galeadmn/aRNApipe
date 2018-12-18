# -*- coding: utf-8 -*-
import random
import os
import time

def get_bsub_arg(pars, L):
    pars = pars.split("/")
    bsub_suffix = list()
    nproc  = int(pars[0])
    if pars[1] != "NA": # memory (Gb)
        mem = str(int(pars[1]) * 1024)
        bsub_suffix.append("-R rusage[mem="+mem+"]")
    if pars[2] != "NA":
        nchild = int(pars[2])
        bsub_suffix.append("-R span[hosts=1]")
    else:
        nchild = nproc
        nproc  = 1
    bsub_suffix = " ".join(bsub_suffix)
    if nchild > L:
        nchild = L
    return nproc, nchild, bsub_suffix


def job_status(path):
## RETURNS THE JOB STATUS GIVEN THE PATH TO THE LOG FILE THAT WILL BE GENERATED BY THE RESOURCE MANAGER
## Input:
##   path: Path to the log file that will be generated by the resource manager. In LSF-based systems this log file is specified using the 'bsub' argument '-o'
## Returns:
##   1:  The job has ended successfully based on the log file generated by the resource manager
##   0:  The job is currently running (log file has not been created)
##   -1: The job exited with error based on the log file generated by the resource manager
    if os.path.exists(path):
        status = ""
        f = open(path, 'r')
        for i in f:
            i = i.strip("\n")
            if i.startswith("Subject:") & i.endswith("Exited"):
                status = "error"
            elif i.startswith("Subject:") and i.endswith("Done"):
                status = "ok"
            else:
                if (status == "error") & i.startswith("TERM_REQUEUE_OWNER"):
                    status = "requeing"
        f.close()
        if status == "ok":
            return 1
        elif status == "error":
            return -1
        else:
            return 0
    else:
        return 0


def job_kill(uds):
## KILLS A JOB GIVEN ITS JOB ID
## Input:
##   uds: Job ID of the job that will be killed
## Returns:
##   Nothing
    os.system("bkill " + uds + " &>/dev/null")


# RETURNS DE PID OF A JOB ACCORDING TO THE OUTPUT OF BSUB SAVED INTO A FILE (PATH)
def get_pid(path, secs):
    time.sleep(secs)
    try:
        f   = open(path, 'r')
        uds = f.readline().rstrip().split(" ")[1].replace(">", "").replace("<", "")
        f.close()
        os.system("rm " + path)
        return uds
    except:
        # error opening or parsing the job submission ouptut file
        return "NA"


def submit_job(wt, n, q, output_log, jid, path2script, script, pt, bsub_suffix):
    ########################################################################
    ## Submits a job script to cluster
    ########################################################################
    ## wt:  cluster walltime
    ## n:   Number of procs  -R rusage[mem=32768]
    ## q:   Queue
    ## output_log:  Path to output log
    ## jid: Preffix
    # path2script: Path to script file
    # script: If 1 script, if 0 command
    r = str(random.randint(0, 10000))
    opts = "-R select[type==any] " + bsub_suffix
    if script == 1:
        command = 'bsub '+opts+' -W ' + wt + " -n " + n + " -q " + q + " -o " + output_log + " -J " + jid + " < " + path2script
    else:
        command = 'bsub '+opts+' -W ' + wt + " -n " + n + " -q " + q + " -o " + output_log + " -J " + jid + " '" + path2script + "'"
    os.system(command + " > " + pt + "/temp/temp" + jid + "_" + r + ".txt")
    uds = get_pid(pt + "/temp/temp" + jid + "_" + r + ".txt", 3)
    return uds