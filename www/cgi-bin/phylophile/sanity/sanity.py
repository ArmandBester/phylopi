#!/usr/bin/env python

#   Phylopi: Purpose built phylogenetic pipeline for the HIV drug
#   resistance testing facility.

#   Copyright (C) 2018 Armand Bester, University of the Free State

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.


import os
from operator import itemgetter

import cgi
import shutil
import cgitb; cgitb.enable()

from pyvirtualdisplay import Display

import rpy2.robjects as r

fn = ''
form = cgi.FieldStorage()
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
                    <title>Sanity</title>
                </head>
                <body>""")
def htmlTail():
          print("""</body>
            </html>""")
       
          
def save_uploaded_file():
    # Get filename here.
    fileitem = form['filename']    
    global fn

    # Test if the file was uploaded
    if fileitem.filename:
        # strip leading path from file name to avoid 
        # directory traversal attacks
        fn = os.path.basename(fileitem.filename)
        fn = str(fn).replace(" ","_")
        open('/var/www/cgi-bin/phylophile/sanity/' + str(fn), 'wb').write(fileitem.file.read())
        message = 'The file ' + fn + ' was uploaded successfully'
    else:
        message = 'No file was uploaded'
        print message

def makeBlastdb():
    makeblastdb_cmd = "../makeblastdb -in HIV1HXB2.fasta -dbtype nucl -title sanity -out sanity"
    os.popen(makeblastdb_cmd)

def doBlastn():
    blastn_cmd = "../blastn -query " + fn  + " -db sanity -evalue 0.1 -outfmt 10 -out sanity.csv -num_threads 4"
    os.popen(blastn_cmd)
    rm_cmd = "rm " + fn
    os.popen(rm_cmd)
    
def sanity():
    r.r('''
        library(ggplot2)
                ''')
    
    r.r('''
        blastTable <- read.csv("sanity.csv", header = FALSE, stringsAsFactors = TRUE)
        colnames(blastTable) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                      "qstart", "qend", "sstart", "send", "evalue", "bitscore")
        
        yRange <- seq(3, nrow(blastTable)+2,1)  
        blastTable <- cbind(blastTable,yRange)


        POL <- data.frame(c(2253),c(2550))
        colnames(POL) <- c("sstart","send")
        
        RT <- data.frame(c(2550),c(3870))
        colnames(RT) <- c("sstart","send")
        
        IN <- data.frame(c(4230),c(5096))
        colnames(IN) <- c("sstart","send")
        
        
        
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        
        
        blastTable$qseqid <- as.factor(blastTable$qseqid)
        blastTable$sstart <- as.integer(blastTable$sstart)
        blastTable$send <- as.integer(blastTable$send)
        
        if (nrow(blastTable) > 7){
            adjSize = 50.0/nrow(blastTable)
          } else{
            adjSize = 6
          }
        
        p = ggplot() +
          geom_segment(data = blastTable, aes(x = sstart,y = yRange, 
                                              xend = send,yend = yRange,
                                              color = blastTable$qseqid)) + 
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                axis.title.y=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                legend.position = "none") +
          labs(x = "Position relative to HXB2", y = NULL)+
          geom_text(data = blastTable, aes(x= (sstart +send)/2, 
                                           y = yRange, label = qseqid), size = adjSize)+
          geom_rect(data = POL, aes(xmin=sstart, xmax=send, ymin = 1, ymax = 1.25), color = "green", fill="green")+
          geom_rect(data = RT, aes(xmin=sstart, xmax=send, ymin = 1, ymax = 1.25), color = "blue", fill = "blue")+
          geom_rect(data = IN, aes(xmin=sstart, xmax=send, ymin = 1, ymax = 1.25), color = "red", fill = "red")
        
        ggsave("sanity.svg", p)
        
                ''')
    
    
    

#main program          
if __name__ == "__main__":
        try:
            
            
            htmlTop()
            save_uploaded_file()
            makeBlastdb()
            doBlastn()
            sanity()
            
            
            png_path = "sanity.svg"
            
            print "<img src=" + png_path + ' height="800" width="800"/>'
            
            print "<p>"
            
            target_url = '''<input type=button onClick="parent.location='http://''' + ip + '''/index.html'" value='Return to input'>'''
            print target_url
            
            htmlTail()
        except:
            cgi.print_exception()


