#!/usr/bin/python

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


import cgi
import cgitb; cgitb.enable()
import multiprocessing
import os
from timeit import Timer
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


def trimal():
    result_dir = log_dict["result_dir"]
    project_name = log_dict["project_name"]
    mafft_aln = log_dict["mafft_aln"]
    trimal_file_path = mafft_aln + "_trimal"
    trimal_cmd = "./trimal -in " + mafft_aln + " -out " + trimal_file_path + " -htmlout " + trimal_file_path + ".html -terminalonly -gt 1"
    
    os.popen(trimal_cmd)
    
    trimal_cleaned_file_path = str(trimal_file_path).replace(".aln","") + ".cleaned.aln"
    trimal_file = open(trimal_file_path,"r").readlines()
    clean_output = open(trimal_cleaned_file_path,"w")
    for line in trimal_file:
        if str(line).find(">") > -1:
            line = line
        else:
            line = line.replace(" ","_")
        clean_output.write(line)
    clean_output.close()
    

#main program
if __name__ == "__main__":
        try:
            htmlTop()
            
            #print information out to the webpage
            print '''            
            <h3>  &#8658; Update blastn database</h3>
            
            <h3>  &#8658; Search blastn database for the specified number of best matches</h3>
            
            <h3>  &#8658; Multiple alignment of submitted and retrieved sequences using MAFFT</h3>
            
            <h3>  &#8658;Cleaning of alignment with trimal</h3>
            
            <h3>Calculating and rendering the distance matrix using R</h3>
            
            <h3>Phylogenetic inference using fasttree</h3>
            
            <h3>Rendering of tree using ete3</h3>
            
            '''
            #time trimal function
            #t = Timer(trimal)
            #trimalTime = t.timeit(1)
            #projectName = log_dict["project_name"]
            #line = projectName +  ",trimalTime," + str(trimalTime) + "\n"
            #timeFile = open("timeFile.csv","a")
            #timeFile.write(line)
            #timeFile.close()
            
            trimal()

            #now redirect to fasttree
            target_url =  '<meta http-equiv="refresh" content="0;url=http://' + ip + '/cgi-bin/phylophile/fasttree.py" /> '
            print target_url
            
            htmlTail()
        except:
            cgi.print_exception()
            print '''
            Consider doing the sanity check on the first page
            <form>
            '''
            target_url = '''<input type=button onClick="parent.location='http://''' + ip + '''/index.html'" value='Return to input'>'''
            print target_url
            print '''
            </form>
            '''
