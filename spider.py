import os
import sys
import optparse
import lib.vcrparser as vcrparser
import lib.spider_stats as spider_stats
import lib.html_lib as html
import traceback

#########################################################################
# PARSER
#########################################################################
desc = "aRNApipe: SPIDER module"
parser = optparse.OptionParser(description = desc)
parser.add_option("-p", "--path", dest = "path", default = "", help = "Required: Path to the project folder")
(opt, args) = parser.parse_args()

#########################################################################
# INITIAL ARRANGEMENTS
#########################################################################
pathscript = os.path.dirname(sys.argv[0]) + "/R/" # PATH TO R SCRIPTS
project, path = html.check_project(opt.path) # CHECKS IF PROJECT FOLDER, SAMPLES FILE AND CONFIGURATION FILE EXIST
html.skeleton(path, os.path.dirname(sys.argv[0])) # CREATES THE SKELETON FOR THE HTML AND OUTPUT RESULTS
config = html.check_config(path) # PARSES THE CONFIGURATION FILE

try:
    samples = vcrparser.get_samples(path.replace(project, ""), project, path + "/samples.list") # PARSES THE SAMPLE FILE AND ASSOCIATED FASTQ FILES
except:
    samples = vcrparser.get_samples(path.replace(project, ""), project, path + "/samples.list", no_check=True)
f = open(path + "/samples.list", 'r')
h = f.readline()
samples_ordered = list()
for i in f:
    if len(i.split('\t')) > 0:
        samples_ordered.append(i.split("\t")[0])
f.close()
lmenu = html.get_menu(config, len(samples))

#########################################################################
# PROCESSES AND ARRANGES LOGS AND OUTPUT DATA FROM STAR, KALLISTO AND HTSEQ
#########################################################################
# TRIMGALORE
try:
    print "Calculating stats for trimgalore"
    spider_stats.stats_trimgalore(path)
except:
    print "Unexpected error generating stats of TRIMGALORE."

# BOWTIE2
btw = False
if config.has_key("bowtie2"):
    btw = True
    try:
        print "Calculating stats for bowtie2"
        spider_stats.stats_bowtie2(path)
    except Exception as ex:
        print traceback.format_exc()
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        print traceback.format_exc()
        print message

# STAR
try:
    spider_stats.summary_star(path + "/results_star/")
    spider_stats.stats_star(path + "/results_star/", samples_ordered, btw)  # GENERATES STATISTICS, COUNTS, RPKM AND ANNOTATION MATRICES
except Exception as ex:
    print traceback.format_exc()
    template = "An exception of type {0} occurred. Arguments:\n{1!r}"
    message = template.format(type(ex).__name__, ex.args)
    print traceback.format_exc()
    print message
# HTSEQ
try:
    spider_stats.stats_htseq(path + "/results_htseq-gene/", samples_ordered, "gene") # GENERATES STATISTICS, COUNTS, RPKM AND ANNOTATION MATRICES
except Exception as ex:
    print traceback.format_exc()
    template = "An exception of type {0} occurred. Arguments:\n{1!r}"
    message = template.format(type(ex).__name__, ex.args)
    print traceback.format_exc()
    print message

#########################################################################
# HTML SUMMARY WEBPAGE
#########################################################################
try:
    print "> Generating webpage with samples list and configuration..."
    print "  - " + path + "/HTML/summary.html"
    html_table = html.print_samples(path,config) # PROVIDES HTML TABLE WITH SAMPLES STATS
    html_table2 = html.config_file(path, "config.txt") # PROVIDES HTML TABLE WITH CONFIGURATION SETTINGS
    html.build_from_template("SUMMARY", project, "", html_table, html_table2, path+"/HTML/summary.html", os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_SUMMARY.html", lmenu)
except:
    print traceback.format_exc()
    print "  - Not ready"

#########################################################################
# HTML DOWNLOADS
#########################################################################
try:
    print "> Generating webpage with download links..."
    print "  - " + path + "/HTML/downloads.html"
    html.build_from_template("DOWNLOADS", project, "", "", "", path+"/HTML/downloads.html", os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_DOWNLOAD.html", lmenu)
except:
    print traceback.format_exc()
    print "  - Not ready"

#########################################################################
# HPC STATS
#########################################################################
try:
    print "> Generating webpage with HPC and LOG statistics..."
    # NOTE: The directory containing results must not be a subdirectory of the root ("/") directory
    #print "Rscript "+pathscript+"/log_stats.R " + path + "/outputs/log_stats.txt"
    #os.system("Rscript "+pathscript+"/log_stats.R " + path + "/outputs/log_stats.txt") # PLOT OF HPC USAGE
    #print "returned from R script"
    #html_table = html.print_table_default(path + "/outputs/log_stats.txt", 2, []) # PROVIDES HTML TABLE WITH HPC STATS
    fils = os.listdir(path + "/outputs/")
    gg = ""
    for i in fils:
        if i.endswith(".png") and i.startswith("log_stats_"):
            gg = gg + '<tr bgcolor="#A8A8A8"><td><center><b>Analysis run number '+i.split("_")[2].replace(".png","")+'</b></center></td></tr><tr bgcolor="#00CC66"><td><img src="../outputs/'+i+'" style="width:800px;"></td></tr>'
    html_table2 = "<table>" + gg + "</table>"
    html.build_from_template("HPC", project, "", html_table, html_table2, path+"/HTML/hpc.html", os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_HPC.html", lmenu)
except Exception as ex:
    print traceback.format_exc()
    template = "An exception of type {0} occurred. Arguments:\n{1!r}"
    message = template.format(type(ex).__name__, ex.args)
    print message

#########################################################################
# TRIM_GALORE
#########################################################################
try:
    if config.has_key("trimgalore"):
        print "> Generating webpage with TrimGalore/Cutadapt statistics..."
        print "  - " + path + "/HTML/trim.html"
        html_table = html.print_table_default(path + "/outputs/stats_trim.txt", -1, []) # PROVIDES HTML TABLE WITH HPC STATS
        data = html.bar_getdata(path + "/outputs/stats_trim_plot.txt",0,[],[])
        html.build_from_template("TrimGalore", project, data, html_table, "", path+"/HTML/trim.html", os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_TRIMG.html", lmenu)
except Exception as ex:
    print traceback.format_exc()
    template = "An exception of type {0} occurred. Arguments:\n{1!r}"
    message = template.format(type(ex).__name__, ex.args)
    print message

#########################################################################
# BOWTIE2
#########################################################################
try:
    if config.has_key("bowtie2"):
        print "> Generating webpage with Bowtie2 statistics..."
        print "  - " + path + "/HTML/bowtie2.html"
        html_table = html.print_table_default(path + "/outputs/stats_bowtie2.txt", -1, []) # PROVIDES HTML TABLE WITH HPC STATS
        data = html.bar_getdata(path + "/outputs/stats_bowtie2_plot.txt",0,[],[])
        html.build_from_template("Bowtie2", project, data, html_table, "", path+"/HTML/bowtie2.html", os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_BOWTIE2.html", lmenu)
except Exception as ex:
    print traceback.format_exc()
    template = "An exception of type {0} occurred. Arguments:\n{1!r}"
    message = template.format(type(ex).__name__, ex.args)
    print message

#########################################################################
# FASTQ
#########################################################################
try:
    if config.has_key("fastqc"):
        print "> Generating webpage with fastqc statistics..."
        print "  - " + path + "/HTML/fastqc.html"
        html_table = spider_stats.stats_fastq(path,samples,config) # PROVIDES HTML TABLE WITH FASTQ STATS
        html.build_from_template("FASTQC", project, "", html_table, "", path+"/HTML/fastqc.html", os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_FASTQC.html", lmenu)
except Exception as ex:
    print traceback.format_exc()
    template = "An exception of type {0} occurred. Arguments:\n{1!r}"
    message = template.format(type(ex).__name__, ex.args)
    print message

#########################################################################
# STAR_QC
#########################################################################
try:
    if config.has_key("star"):
        print "> Generating webpage with STAR statistics..."
        print "  - " + path + "/HTML/star.html"
        #html_table = html.print_table_default(path + "/outputs/star_unstranded_stats.txt", -1, []) # PROVIDES HTML TABLE WITH HPC STATS
        #data = html.bar_getdata (path + "/outputs/star_unstranded_stats.txt",0,range(1,6),[])
        #html.build_from_template("STAR", project, data, html_table, "", path+"/HTML/star.html", os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_STAR.html", lmenu)
        html_table = html.print_table_default(path + "/outputs/star_summary.txt", -1, []) # PROVIDES HTML TABLE WITH HPC STATS
        data = html.bar_getdata (path + "/outputs/star_summary.txt",0,range(1,8),[])
        html.build_from_template("STAR", project, data, html_table, "", path+"/HTML/star.html", os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_GALE_STAR.html", lmenu)
except Exception as ex:
    print traceback.format_exc()
    template = "An exception of type {0} occurred. Arguments:\n{1!r}"
    message = template.format(type(ex).__name__, ex.args)
    print message
    print "  - Not ready"

#########################################################################
# HTSEQ_QC
#########################################################################
try:
    for ij in ["htseq-gene", "htseq-exon"]:
        if config.has_key(ij):
            print "> Generating webpage with "+ij+" statistics..."
            print "  - " + path + "/HTML/"+ij+".html"
            html_table = html.print_table_default(path + "/outputs/"+ij+"_stats.txt", -1, []) # PROVIDES HTML TABLE WITH HPC STATS
            data = html.bar_getdata (path + "/outputs/"+ij+"_stats.txt",0,range(1,7),[])
            html.build_from_template(ij.upper(), project, data, html_table, "", path+"/HTML/" + ij + ".html", os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_HTSEQ.html", lmenu)
except Exception as ex:
    print traceback.format_exc()
    template = "An exception of type {0} occurred. Arguments:\n{1!r}"
    message = template.format(type(ex).__name__, ex.args)
    print message
    print "  - Not ready"

#########################################################################
# STAR & HTSEQ STATS ON COUNTS/RPKMS
#########################################################################
try:
    if len(samples) > 1:
        strandedness = config["programs"]["strandedness"].strip()
        if strandedness == "yes":
            n = {"star":["STAR","star_stranded"],"htseq-gene":["HTseq-count Gene", "htseq-gene"],"htseq-exon":["HTseq-count Exon", "htseq-exon"]}
        elif strandedness == "no":
            n = {"star":["STAR","star_unstranded"],"htseq-gene":["HTseq-count Gene", "htseq-gene"],"htseq-exon":["HTseq-count Exon", "htseq-exon"]}
        elif strandedness == "reverse":
            n = {"star":["STAR","star_stranded-reverse"],"htseq-gene":["HTseq-count Gene", "htseq-gene"],"htseq-exon":["HTseq-count Exon", "htseq-exon"]}
        for prog, pname in n.iteritems():
            if config.has_key(prog):
                #if os.path.isfile( path + "/outputs/" + pname[1]):
                os.system("Rscript "+pathscript+"/stats_algs.R " + path + "/outputs/ " + pname[1]) # PLOT OF HPC USAGE
                html_table = html.print_table_default(path + "/outputs/" + pname[1] + "_pca.txt", -1, [0, 1, 2, 3, 4, 6, 7, 8, 12, 13, 14, 15, 17, 18, 19])
                html.build_amcharts(os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_PCA.html", path + "/HTML/" + prog + "2.html", prog, pname, path, html_table, project, lmenu)
except Exception as ex:
    print traceback.format_exc()
    template = "An exception of type {0} occurred. Arguments:\n{1!r}"
    message = template.format(type(ex).__name__, ex.args)
    print message
    print "  - Not ready"
