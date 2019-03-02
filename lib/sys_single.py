# -*- coding: utf-8 -*-
import os
import subprocess
import time

def get_bsub_arg(pars, L):
    pars = pars.split("/")
    nproc  = int(pars[0])
    return nproc, 1, ""


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
    return 1


# RETURNS DE PID OF A JOB ACCORDING TO THE OUTPUT OF BSUB SAVED INTO A FILE (PATH)
def get_pid(path, secs):
    return 0


def submit_job(wt, n, q, output_log, jid, path2script, script, pt):
    ########################################################################
    ## Submits a job script to CLUSTER
    ########################################################################
    ## wt:  Cluster walltime
    ## n:   Number of procs  -R rusage[mem=32768]
    ## q:   Queue
    ## output_log:  Path to output log
    ## jid: Preffix
    # path2script: Path to script file
    # script: If 1 script, if 0 command
    print path2script
    if script == 1:
        L = ["bash"]
        path2script = path2script.split(" ")
        for i in path2script:
            L.append(i)
        p = subprocess.Popen(L, shell = False)
    else:
        L = []
        path2script = path2script.split(" ")
        for i in path2script:
            L.append(i)
        p = subprocess.Popen(path2script)
    uds = str(p.pid)
    p = scan_process(p, 10, output_log)
    if p != 0:
        exit("Error executing: " + str(path2script))
    return uds


def scan_process(p, secs, log):
    dt = time.time()
    out = open(log, 'w')
    print >> out, "Started at " + time.strftime("%a %B %d %H:%M:%S %Y")
    out.close()
    vectMemo = [0.0, 0.0, 0.0]
    while p.poll() is None:
        time.sleep(secs)
        r = get_memo(p.pid)
        if r >= 0:
            vectMemo[0] +=1 # number of reads
            vectMemo[1] += r # summ
            if r > vectMemo[2]:
                vectMemo[2] = r # max
    out = open(log, 'a')
    print >> out, "Results reported on " + time.strftime("%a %B %d %H:%M:%S %Y")
    if p.returncode == 0:
        print >> out, "Subject: Done"
        print >> out, "Successfully completed."
    else:
        print >> out, "Subject: Exited"
        print >> out, "Error"
    print >> out, " CPU time : " + str(round(float(time.time() - dt),2)) + " sec."
    print >> out, " Max Memory : " + str(round(vectMemo[2],2)) + " MB"
    print >> out, " Average Memory : " + str(round(float(vectMemo[1])/vectMemo[0],2)) + " MB"
    print >> out, " Max Swap : 0 MB"
    out.close()
    return p.returncode


def get_memo(pid):
    try:
        pp = subprocess.Popen("ps aux | grep " + str(pid), shell = True, stdout=subprocess.PIPE)
        a = pp.communicate()
        a = a[0].split("\n")
        for i in a:
            if i.split()[1] == str(pid):
                return float(i.split()[5])/1000
        return -1
    except:
        return -1
