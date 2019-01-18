#!/usr/bin/python

import cgi
import cgitb
cgitb.enable(format='text')
import os
from ete3 import Tree, TreeStyle, NodeStyle
import multiprocessing
import sys
import subprocess
from timeit import Timer
#import rpy2.robjects as r
import netifaces

ifaces_list = netifaces.interfaces()
ifaces = []
ifaces.append(ifaces_list[1])
ifaces.append(ifaces_list[2])
ifaces.append(ifaces_list[0])
ip = ''
for iface in ifaces:
    #prefer eth
    if str(iface).find("eth") > -1 and netifaces.ifaddresses(iface).has_key(2) == True:
        ip = netifaces.ifaddresses(iface)[2][0]["addr"]
        break
    if str(iface).find("wlan") > -1 and netifaces.ifaddresses(iface).has_key(2) == True:
        ip = netifaces.ifaddresses(iface)[2][0]["addr"]
        break
    if str(iface).find("lo") > -1:
        ip = netifaces.ifaddresses(iface)[2][0]["addr"]
        break


from pyvirtualdisplay import Display

log_dict = {}
tmp_log = open("input/tmp.log","r").readlines()
#lets put this in a dictionary
log_dict = {}
for line in tmp_log:
    line = str(line).strip()
    if str(line).find("=") > -1:
        line = str(line).split("=")
        key = line[0]
        value = line[1]
        log_dict[key] = value

def htmlTop():
    print("""Content-type:text/html\n\n
            <!DOCTYPE html>
            <html lang="en">
                <head>

                <style>
                    body {
                        background-color: #E8E8E8 ;
                    }
                    h2 {
                        color: black;
                        text-align: center;
                    }
                    p {
                        font-family: "Arial";
                        font-size: 20px;
                    }
                    </style>
                    <meta charset="utf-8"/>
                    <title>PhyloPi running</title>
                </head>
                <body>""")

def htmlTail():
          print("""</body>
            </html>""")

def render():
    result_dir = log_dict["result_dir"]
    project_name = log_dict["project_name"]
    newick_file_path = result_dir  + project_name + "_fastTree_tree.txt"
    t = Tree(open(newick_file_path).read())

    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = True
    ts.show_branch_support = True

    tree_out_pdf = newick_file_path = result_dir  + project_name + ".pdf"
    tree_out_png = newick_file_path = result_dir  + project_name + ".png"
    tree_out_svg = newick_file_path = result_dir  + project_name + ".svg"


    nbStyle = NodeStyle()
    nbStyle["fgcolor"] = "red"
    nbStyle["size"] = 10


    #find the primary samples, this can be found in the input fasta
    input_fasta = open(log_dict["fasta_input"],"r").readlines()
    fasta_input_list = []  #these are the primary samples
    for line in input_fasta:
        if str(line).find(">") > -1: #if fasta header is found
            fasta_input_list.append(str(line).replace(">","").replace("\n","").replace("\r","").replace(" ",""))

    #to make things easier, we want to mark samples in the current analysis
    for node in t.traverse():
        for nbSample in fasta_input_list:
            if  nbSample in node.name:
                node.set_style(nbStyle)


    t.render(tree_out_pdf,tree_style=ts)
    t.render(tree_out_png,tree_style=ts)
    t.render(tree_out_svg,tree_style=ts)


#def to create a R script to create a heatmap
def heatmap():
    result_dir = log_dict["result_dir"]
    project_name = log_dict["project_name"]
    aln_file_path = result_dir  + project_name
    fasta_input = log_dict["fasta_input"]
    RscriptFile = open("heatmap_Rscript.R", "w")


    rCode1 = '''
        library(seqinr)
        library(ape)
        library(d3heatmap)
        library(reshape2)
        library(RColorBrewer)
        library(htmlwidgets)
        library(jsonlite)
        library(yaml)
        library(plotly)
        '''

    rCode2 =  "\taln <- read.alignment('" + aln_file_path + "-mafft_trimal.cleaned.aln', format = 'fasta')"

    rCode3 = '''
        alnBin <- as.DNAbin(aln)
        alnDist <- dist.dna(alnBin)
        alnDist <- as.matrix(dist.dna(alnBin))
        alnDistLong <- melt(alnDist)
        alnDistLong$sameSample <- alnDistLong$Var1 == alnDistLong$Var2
        alnDistLong <- subset(alnDistLong, sameSample %in% FALSE)
        alnDistLong <- alnDistLong[,-4]
        colnames(alnDistLong) <- c('sample1', 'sample2', 'distance')
        alnDistLong$combined <- paste(alnDistLong$sample1, alnDistLong$sample2, sep = ",")
        alnDistLong$sorted <- apply(alnDistLong,
                                    1,
                                    function(x) {paste(sort(unlist(strsplit(x[4], split=","))), collapse = ",")})
        alnDistDedup <- alnDistLong[!duplicated(alnDistLong$sorted),]
        alnDistDedup <- alnDistDedup[, c(-4,-5)]
        alnDistDedup <- alnDistDedup[order(alnDistDedup$distance),]

        alnDist <- data.frame(alnDist,check.names = FALSE) #new
    '''
    rCode4 = '\tqSample <- names(read.fasta("' + fasta_input + '"))'

    rCode5 = '''
        rowNames <- rownames(alnDist)
        colNames <- colnames(alnDist)

        t <- list(
        family = "sans serif",
        size = 10,
        color = toRGB("black"))
        anno.df <- data.frame(q = character(),
                      rowNames = character())

        for (q in qSample){
            tmp <- cbind(q, rowNames)
            anno.df <- rbind(anno.df, tmp)
        }

        focus.df <- alnDistLong
        focus.df$prim <- focus.df$sample1 %in% qSample

        focus.df <- focus.df %>% 
          filter(prim == TRUE | sample1 == sample2)

        p <- plot_ly(x = focus.df$sample2,
                    y = focus.df$sample1,
                    z = focus.df$distance,
                    type = "heatmap", colors = brewer.pal(11, "RdYlGn"), 
                    zmin = 0.0, zmax = 0.03,  xgap = 2, ygap = 1) %>% 
            layout(margin = list(l = 100, r = 10, b = 100, t = 10, pad = 4), 
                             yaxis = list(tickfont = list(size = 10), showspikes = TRUE),
                             xaxis = list(tickfont = list(size = 10), showspikes = TRUE))
'''




    R_table_out = 'write.csv(alnDist,"/var/www/cgi-bin/phylophile/' + aln_file_path + '_heatmap.csv")'
    R_long_table_out = 'write.csv(alnDistDedup,"/var/www/cgi-bin/phylophile/' + aln_file_path + '_dists.csv", row.names = FALSE)'

    R_save_cmd = "htmlwidgets::saveWidget(as_widget(p), '/var/www/cgi-bin/phylophile/"\
                 + aln_file_path + "_heatmap.html')"

    RscriptFile.write(rCode1)
    RscriptFile.write("\n")
    RscriptFile.write(rCode2)
    RscriptFile.write("\n")
    RscriptFile.write(rCode3)
    RscriptFile.write("\n")
    RscriptFile.write(rCode4)
    RscriptFile.write("\n")
    RscriptFile.write(rCode5)
    RscriptFile.write("\n")
    RscriptFile.write(R_table_out)
    RscriptFile.write("\n")
    RscriptFile.write(R_long_table_out)
    RscriptFile.write("\n")
    RscriptFile.write(R_save_cmd)
    RscriptFile.close()
    #Rscript_cmd = "Rscript heatmap_Rscript.R" 
    os.popen("Rscript heatmap_Rscript.R")

#main program
if __name__ == "__main__":
        try:
            htmlTop()

            #print information out to the webpage
            print '''
            <h3>  &#8658; Update blastn database</h3>

            <h3>  &#8658; Search blastn database for the specified number of best matches</h3>

            <h3>  &#8658; Multiple alignment of submitted and retrieved sequences using MAFFT</h3>

            <h3>  &#8658; Cleaning of alignment with trimal</h3>

            <h3>  &#8658; Calculating and rendering the distance matrix using R</h3>

            <h3>  &#8658; Phylogenetic inference using fasttree</h3>

            <h3>  &#8658; Rendering of tree using ete3</h3>

            '''



            display = Display(visible=0, size=(1280, 948))
            display.start()
        #time heatmap
            #t1 = Timer(heatmap)
            #heatmapTime = t1.timeit(1)
            #projectName = log_dict["project_name"]
            #line = projectName +  ",heatmapTime," + str(heatmapTime) + "\n"
            #timeFile = open("timeFile.csv","a")
            #timeFile.write(line)
            #timeFile.close()

            heatmap()

            #time render
            #t2 = Timer(render)
            #renderTime = t2.timeit(1)

            #projectName = log_dict["project_name"]

            #line = projectName +  ",renderTime," + str(renderTime) + "\n"
            #timeFile = open("timeFile.csv","a")
            #timeFile.write(line)
            #timeFile.close()

            render()
            display.stop()

            target_url = '<meta http-equiv="refresh" content="0;url=http://' + ip + '/cgi-bin/phylophile/zip_and_archive.py" /> '
            print target_url
            htmlTail()

        except:
            cgi.print_exception()
