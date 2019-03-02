import shutil
import os

modules = ["trimgalore", "bowtie2", "fastqc", "star", "htseq-gene", "htseq-exon"]
module_names = {"trimgalore":"",
                "bowtie2":"",
                "fastqc":"",
                "star":"",
                "htseq-gene":"",
                "htseq-exon":""}


def get_menu(config, ns):
    enabled = dict()
    for module in modules:
        if config.has_key(module):
            if int(config[module][0].split("/")[0]) > 0:
                enabled[module] = 1
    menu = list()
    menu.append('<h1><a #highlight="" href="./summary.html">Summary</a></h1>')
    menu.append('<h1><a #highlight="" href="./hpc.html">HPC statistics</a></h1>')
    if enabled.has_key("fastqc") or enabled.has_key("trimgalore"):
        menu.append('<h1>Raw QC:</h1>')
    if enabled.has_key("trimgalore"):
        menu.append('<h2><a #highlight="" href="./trim.html">- TrimGalore/Cutadapt</a></h2>')
    if enabled.has_key("bowtie2"):
        menu.append('<h2><a #highlight="" href="./bowtie2.html">- Bowtie2</a></h2>')
    if enabled.has_key("fastqc"):
        menu.append('<h2><a #highlight="" href="./fastqc.html">- FastQC</a></h2>')
    if enabled.has_key("star") or enabled.has_key("picard") or enabled.has_key("htseq-gene") or enabled.has_key("htseq-exon"):
        menu.append('<h1>Alignment QC:</h1>')
    if enabled.has_key("star"):
        menu.append('<h2><a #highlight="" href="./star.html">- STAR</a></h2>')
    if enabled.has_key("htseq-gene"):
        menu.append('<h2><a #highlight="" href="./htseq-gene.html">- HTseq-Gene</a></h2>')
    if enabled.has_key("htseq-exon"):
        menu.append('<h2><a #highlight="" href="./htseq-exon.html">- HTseq-Exon</a></h2>')
    if ns > 1:
        if enabled.has_key("star") or enabled.has_key("htseq-gene") or enabled.has_key("htseq-exon"):
            menu.append('<h1>Count statistics:</h1>')
            menu.append('<h2><a #highlight="" href="./downloads.html">- DOWNLOADS</a></h2>')
        if enabled.has_key("star"):
            menu.append('<h2><a #highlight="" href="./star2.html">- STAR</a></h2>')
        if enabled.has_key("htseq-gene"):
            menu.append('<h2><a  href="./htseq-gene2.html">- HTseq-Gene</a></h2>')
        if enabled.has_key("htseq-exon"):
            menu.append('<h2><a #highlight="" href="./htseq-exon2.html">- HTseq-Exon</a></h2>')
    menu = "\n".join(menu)
    return menu

def print_samples(path,config):
    analysis = ['trimgalore', 'bowtie2', 'fastqc', 'star', "htseq-gene", "htseq-exon"]
    sta= {"trimgalore":"TrimGalore", "bowtie2":"Bowtie2", "fastqc":"FastQC","star":"STAR", "htseq-gene":"HTseq-gene",
          "htseq-exon":"HTseq-exon"}
    # SAMPLES LIST
    samples = dict()
    f = open(path + "/samples.list",'r')
    hss = f.readline().strip("\n").split("\t")
    idx = []
    for i, ix in enumerate(hss):
        if ix.startswith('FASTQ'):
            idx.append(i)
    hs = [hss[0]] + [hss[i] for i in idx]
    for i in f:
        i = i.strip("\n").split("\t")
        if i[0] != "":
            samples[i[0]] = [i[j] for j in idx]
    f.close()
    # scan
    results  = dict()
    for i in analysis:
        if config.has_key(i):
            sok = dict()
            if os.path.exists(path + "/results_" + i + "/samples_ok.txt"):
                f = open(path + "/results_" + i + "/samples_ok.txt", 'r')
                for j in f:
                    sok[j.strip("\n")] = 1
                f.close()
            results[i] = dict()
            if i == "trimgalore":
                for x, y in sorted(samples.iteritems()):
                    res = []
                    if sok.has_key(x):
                        filess = y
                        for f in filess:
                            f = f.split("/")[-1]
                            link = "../results_trimgalore/" + f + "_trimming_report.txt"
                            link = '<a href="LINK" target="_blank">OK</a>'.replace("LINK", link)
                            res.append(link)
                    else:
                        res.append("FAIL")
                    results["trimgalore"][x] = " / ".join(res)
            if i == "bowtie2":
                for x, y in sorted(samples.iteritems()):
                    res = []
                    if sok.has_key(x):
                        filess = y
                        for f in filess:
                            f = f.split("/")[-1]
                            link = "../results_bowtie2/" + f + "_bowtie2_report.txt"
                            link = '<a href="LINK" target="_blank">OK</a>'.replace("LINK", link)
                            res.append(link)
                    else:
                        res.append("FAIL")
                    results["bowtie2"][x] = " / ".join(res)
            elif i=="fastqc":
                for x, y in sorted(samples.iteritems()):
                    res    = []
                    if sok.has_key(x):
                        filess = y
                        for f in filess:
                            f = f.split("/")[-1]
                            link = "../results_fastqc/"+f.replace(".fastq","").replace(".gz","")+"_fastqc/fastqc_report.html"
                            link = '<a href="LINK" target="_blank">OK</a>'.replace("LINK",link)
                            res.append(link)
                    else:
                        res.append("FAIL")
                    results["fastqc"][x] = " / ".join(res)
            elif i=="star":
                for x, y in sorted(samples.iteritems()):
                    res = []
                    if sok.has_key(x):
                        link = "../results_star/" + x + "_Aligned.out.sam"
                        link = '<a href="LINK" target="_blank">BAM-OK</a>'.replace("LINK", link)
                        res.append(link)
                        link = "../results_star/" + x + "_ReadsPerGene.out.tab"
                        link = '<a href="LINK" target="_blank">COUNTS-OK</a>'.replace("LINK", link)
                        res.append(link)
                        link = "../results_star/" + x + "_SJ.out.tab"
                        link = '<a href="LINK" target="_blank">SJ-OK</a>'.replace("LINK", link)
                        res.append(link)
                    else:
                        res.append("BAM-FAIL")
                        res.append("COUNTS-FAIL")
                        res.append("COUNTS-FAIL")
                    results["star"][x] = " / ".join(res)
            elif i=="htseq-gene":
                for x, y in sorted(samples.iteritems()):
                    res = []
                    if sok.has_key(x):
                        link = "../results_htseq-gene/" + x + ".tab"
                        link = '<a href="LINK" target="_blank">OK</a>'.replace("LINK", link)
                        res.append(link)
                    else:
                        res.append("FAIL")
                    results["htseq-gene"][x] = " / ".join(res)
            elif i=="htseq-exon":
                for x, y in sorted(samples.iteritems()):
                    res = []
                    if sok.has_key(x):
                        link = "../results_htseq-exon/" + x + ".tab"
                        link = '<a href="LINK" target="_blank">OK</a>'.replace("LINK", link)
                        res.append(link)
                    else:
                        res.append("FAIL")
                    results["htseq-exon"][x] = " / ".join(res)
    n = "<th bgcolor='#A8A8A8'>Sample</th>"
    for i in hs[1:]:
        n = n + "<th bgcolor='#A8A8A8'> Size "+i+" (Gb)</th>"
    for i in range(len(analysis)):
        if results.has_key(analysis[i]):
            n = n +"<th bgcolor='#A8A8A8'>"+sta[analysis[i]]+"</th>"
    thead = "<thead><tr>"+n+"</tr></thead>"
    tab = list()
    for i in sorted(samples.keys()):
        n = ["<td bgcolor='#A8A8A8'>"+i+"</td>"]
        for j in range(len(hs[1:])):
            try:
                n.append("<td bgcolor='#A8A8A8'>" + str(round(os.stat(samples[i][j]).st_size/1000000000.0, 2)) + "</td>")
            except:
                n.append("<td bgcolor='#A8A8A8'>NA</td>")
        for a in analysis:
            if results.has_key(a):
                cl = "#00CC66"
                if "FAIL" in results[a][i]:
                    cl = "#CC3300"
                n.append("<td bgcolor='"+cl+"'>"+results[a][i]+"</td>")
        tab.append("<tr>"+"".join(n)+"</tr>")
    return '<table id="DT" class="display">' + thead + "<tbody>" + "\n".join(tab) + "</tbody></table>"


def config_file(path, fname):
    f   = open(path + "/" + fname,'r')
    n   = list()
    for i in f:
        i = i.strip("\n").split("\t")
        if len(i) >= 2:
            g = "<tr><td bgcolor='#00CC66'>"+i[0]+"</td><td bgcolor='#00CC66'>"+i[1]+"</td></tr>"
        else:
            g = "<tr><td colspan='2' bgcolor='#C0C0C0'>" + i[0] + "</td></tr>"
        n.append(g)
    f.close()
    tab = "<table>"+"".join(n)+"</table>"
    return tab


def print_table_default(datafile, index, select):
    print datafile
    palette = ["#00FA9A", "#AFEEEE", "#D8BFD8", "#DEB887", "#D3D3D3", "#EEE8AA"]
    if not os.path.exists(datafile):
        print "returning empty string"
        return ""
    f = open(datafile, 'r')
    h = f.readline().strip("\n").split("\t")
    if len(select) == 0:
        select = range(len(h))
    n = ""
    for i in select:
        n = n + "<th align='center' bgcolor='#A8A8A8'>" + h[i] + "</th>"
    if index < 0:
        n = '<table id="DT" class="display"><thead><tr>' + n + '</tr></thead><tbody>'
    else:
        n = '<table><thead><tr>' + n + '</tr></thead><tbody>'
    M = dict()
    r = 0
    for i in f:
        i = i.strip("\n").split("\t")
        if len(i) > 0:
            temp = ""
            if index > -1:
                if not M.has_key(i[index]):
                    M[i[index]] = palette[r % len(palette)]
                    r += 1
            for k in select:
                j = i[k]
                if (j == "-1") or (j.startswith("NA ") or j.endswith(" NA") or j == "NA"):
                    temp = temp + "<td align='center' bgcolor='#CC3300'>" + j + "</td>"
                else:
                    if index < 0:
                        temp = temp + "<td align='center' bgcolor='#00CC66'>" + j + "</td>"
                    else:
                        temp = temp + "<td align='center' bgcolor='" + M[i[index]] + "'>" + j + "</td>"
            n = n + "<tr>" + temp + "</tr>"
    n += "</tbody></table>"
    return n


def check_project(path):
    print "> Checking project path..."
    if path == "":
        exit("Parameter -p is required.")
    if not os.path.exists(path):
        exit("Path to project folder not found.")
    if path.endswith("/"):
        path = path[0:(len(path)-1)]
    if not os.path.exists(path + "/config.txt"):
        exit("Project configuration file not found: " + path + "/config.txt")
    if not os.path.exists(path + "/samples.list"):
        exit("Project samples file not found: " + path + "/samples.list")
    project = path.split("/")
    project = project[len(project)-1] # project id
    return project, path


def skeleton(path, path2html):
    print "> Building HTML and OUTPUT folders skeletons..."
    print "  - Path: " + path
    print "  - Libs: " + path2html
    # Creates output directiories
    if os.path.exists(path + "/outputs"):
        os.system("rm -r " + path + "/outputs")
    os.mkdir(path + "/outputs")
    # Creates HTML directories
    n = os.listdir(path)
    if "HTML" in n:
        os.system("rm -r " + path + "/HTML")
    if not ("HTML" in os.listdir(path)):
        os.mkdir(path + "/HTML")
    os.mkdir(path + "/HTML/html")
    # Copy lightbox and jquery
    shutil.copy(path2html + "/html/style.css", path + "/HTML/html/style.css")
    shutil.copy(path2html + "/html/lytebox.js", path + "/HTML/html/lytebox.js")
    shutil.copy(path2html + "/html/lytebox.css", path + "/HTML/html/lytebox.css")
    shutil.copy(path2html + "/html/jquery-1.12.0.js", path + "/HTML/html/jquery-1.12.0.js")
    shutil.copy(path2html + "/html/jquery.dataTables.min.js", path + "/HTML/html/jquery.dataTables.min.js")
    shutil.copy(path2html + "/html/jquery.dataTables.min.css", path + "/HTML/html/jquery.dataTables.min.css")
    shutil.copy(path2html + "/html/dataTables.colReorder.min.js", path + "/HTML/html/dataTables.colReorder.min.js")
    shutil.copy(path2html + "/html/amcharts.js", path + "/HTML/html/amcharts.js")
    shutil.copy(path2html + "/html/serial.js", path + "/HTML/html/serial.js")
    shutil.copy(path2html + "/html/xy.js", path + "/HTML/html/xy.js")
    os.system("cp -r " + path2html + "/html/images " + path + "/HTML/html/")

def build_amcharts(input, output, prog, pname, path, html_table, project, lmenu):
    out = open(output, 'w')
    f   = open(input, 'r')
    it = ""
    for i in f:
        i = i.strip("\n")
        if i.startswith("#ITERATOR"):
            pattern = i.replace("#ITERATOR=","")
            if os.path.isfile(path + "/outputs/" + pname[1] + "_pca.txt"):
                f2 = open(path + "/outputs/" + pname[1] + "_pca.txt", 'r')
                h  = f2.readline()
                it = list()
                for j in f2:
                    np = pattern
                    j = j.strip("\n").split("\t")
                    for k in range(len(j)):
                        np = np.replace("#VAR"+str(len(j)-k-1),j[len(j)-k-1])
                    it.append(np)
                it = ", ".join(it)
                f2.close()
        else:
            if prog + "2.html" in i:
                i = i.replace("#HIGHLIGHT", 'style="color:#808080"')
            print >> out, i.replace("#LATMENU",lmenu).replace("#PROG", pname[0]).replace("#PROJECT", project).replace("#SITERATOR", it).replace("#HIGHTLIGHT", "").replace("#TABLE", html_table)
    f.close()
    out.close()


def check_config(path):
    # Parses the configuration file
    print "> Parsing configuration file..."
    try:
        z = ["trimgalore", "bowtie2", "fastqc", "star", "htseq-gene", "htseq-exon"]
        f = open(path + "/config.txt", 'r')
        analysis = dict()
        analysis["cluster"] = dict()
        analysis["programs"] = dict()
        for i in f:
            if not i.startswith("#"):
                i = i.strip("\n").split("\t")
                if len(i) > 1:
                    if i[0] in z:
                        if int(i[1].split("/")[0]) > 0:
                            analysis[i[0]] = [i[1], "results_" + i[0], dict()]
                    elif i[0] in ["wt", "q"]:
                        analysis["cluster"][i[0]] = i[1]
                    elif i[0] == "star_args_own":
                        if analysis["programs"]["star_args"] == "own":
                            analysis["programs"]["star_args"] = i[1]
                    elif i[0] == "genome_build":
                        analysis[i[0]] = i[1]
                    else:
                        analysis["programs"][i[0]] = i[1]
        f.close()
        return analysis
    except:
        exit("Error checking configuration file: " + path + "/config.txt")


# def print_config(config,path):
#     table = list()
#     table.append(["Analysis", "Processors", "Folder", "Timestamp",
#                   "TStart","TEnd","Success","CPU-Time","MaxMemo",
#                   "AveMemo","MaxSwap","Parameters"])
#     for i in modules:
#         if config.has_key(i):
#             n = check_log_cluster(path, config[i][1])
#             st = []
#             if len(config[i][2]) > 0:
#                 for v,w in config[i][2].iteritems():
#                     st.append(v+": "+w)
#             st = "<br>".join(st)
#             for j in range(len(n)):
#                 tt = [module_names[i]]+[config[i][0],"./"+config[i][1]]+n[j]+[st]
#                 table.append(tt)
#     n = ""
#     for i in table[0]:
#         n = n+"<th bgcolor='#A8A8A8'>"+i+"</th>"
#     n = ["<tr>"+n+"</tr>"]
#     for i in table[1:]:
#         temp = ""
#         for j in i:
#             if "NA" in j:
#                 temp = temp+"<td bgcolor='#CC3300'>"+j+"</td>"
#             else:
#                 temp = temp+"<td bgcolor='#00CC66'>"+j+"</td>"
#         n.append("<tr>"+temp+"</tr>")
#     return n


# def check_log_cluster(path, val):
#     t = os.listdir(path)
#     if not (val in t):
#         return ["NA","NA","NA","NA","NA","NA","NA"]
#     t  = os.listdir(path + "/" + val)
#     t2 = list()
#     for i in t:
#         if i.startswith("log_cluster_") and (("scheduler" in i) == False):
#             if len(i.split("_")) >= 4 :
#                 t2.append(i)
#     if len(t2) == 0:
#         return [["NA","NA","NA","NA","NA","NA","NA","NA"]]
#     n = list()
#     for jv in sorted(t2):
#         f = open(path + "/" + val + "/" + jv,'r')
#         ts = ""
#         te = ""
#         suc= "No"
#         cpu_time = ""
#         max_memo = ""
#         ave_memo = ""
#         max_swap = ""
#         pid      = "_".join(jv.split("_")[2:4]).replace(".txt","")
#         for i in f:
#             if i.startswith("Started at"):
#                 ts = i.rstrip().replace("Started at ","")
#             if i.startswith("Results reported on"):
#                 te = i.rstrip().replace("Results reported on ","")
#             if i.startswith("Successfully completed."):
#                 suc = "Yes"
#             if "CPU time :" in i:
#                 cpu_time = " ".join(i.rstrip().split()[3:5])
#             if "Max Memory :" in i:
#                 max_memo = " ".join(i.rstrip().split()[3:5])
#             if "Average Memory :" in i:
#                 ave_memo = " ".join(i.rstrip().split()[3:5])
#             if "Max Swap :" in i:
#                 max_swap = " ".join(i.rstrip().split()[3:5])
#         f.close()
#         n.append([pid,ts,te,suc,cpu_time,max_memo,ave_memo,max_swap])
#     return n


def bar_getdata (filename, head, cols_bar, cols_line):
    # LOAD DATA
    if not os.path.exists(filename):
        return ""
    f = open(filename, 'r')
    h = f.readline().strip("\n").split("\t")
    if (len(cols_bar)==0) and (len(cols_line)==0):
        cols_bar = range(1, len(h))
    D = list()
    for i in f:
        i = i.strip("\n").split("\t")
        n = list()
        for j in range(len(h)):
            try:
               if j == head:
                   n.append('"' + h[j] + '": "' + i[j] + '"')
               elif (j in cols_bar) or (j in cols_line):
                   if i[j] != "NA":
                       n.append('"' + h[j] + '": ' + i[j])
                   else:
                       n.append('"' + h[j] + '": 0')
            except IndexError:
               continue
        D.append("{" + ", ".join(n) + "}")
    D = "var chartData = [" + ", ".join(D) + "];"
    return D


def build_from_template(prog, project, data, html_table, html_table2, output, template, lmenu):
    out = open(output,'w')
    f = open(template, 'r')
    r = output.split("/")[-1]
    for i in f:
        i = i.strip("\n")
        if r in i:
            i = i.replace("#HIGHLIGHT", 'style="color:#808080"')
        print >> out, i.replace("#LATMENU",lmenu).replace("#PROG", prog).replace("#PROJECT", project).replace("#DATA", data).replace("#HIGHTLIGHT", "").replace("#TABLE2", html_table2).replace("#TABLE", html_table)
    f.close()
    out.close()
