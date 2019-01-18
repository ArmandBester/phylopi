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
import shutil
import os
import time
import cgitb; cgitb.enable()
from collections import OrderedDict
import subprocess
import sys
import multiprocessing
import netifaces
from timeit import Timer

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


form = cgi.FieldStorage()

log_dict = {}

cpus_avail = str(multiprocessing.cpu_count()-1)  #for blast
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
#start by parsing some info from the tmp_log file
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

#get values from the conf file in bastnDB folder
conf = open("blastnDB/config.conf","r").readlines()
blastnDB = ''
original_fasta = ''
for line in conf:
    if str(line).find("db_name") > -1:
        blastnDB = str(line).split("=")[1]
    if str(line).find("fasta") > -1:
        original_fasta = str(line).split("=")[1]

def blast():        
    ##blast analysis
    fasta_input = log_dict["fasta_input"]
    reps = log_dict["reps"]
    result_dir = log_dict["result_dir"]
    project_name = log_dict["project_name"]
    blast_database = log_dict["blast_database"]
    blast_cmd = "./blastn -query " + fasta_input + " -task blastn -db blastnDB/" + blast_database + " -word_size 7 -outfmt 10 -num_alignments " + str(reps) + " -num_threads " + cpus_avail + " -out " + result_dir + project_name + ".blastn"
    os.system(blast_cmd)
       
    ######################################################
    ########## parse the blastn results ##################
    blastn_result_file = result_dir + project_name + ".blastn"
    blastn_result = open(blastn_result_file).readlines()    
    blastn_query_list = []
    blastn_subject_list = []
    for line in blastn_result:
        line = str(line).replace("\n","").split(",")
        blastn_query = line[0]
        blastn_subject = line[1]
        blastn_query_list.append(blastn_query)
        blastn_subject_list.append(blastn_subject)
    #the subject list might contain duplicates, which we will remove
    #just to make sure that nothing slips throught, do the same with the query list
    blastn_query_list = sorted(set(blastn_query_list),key=blastn_query_list.index)
    blastn_subject_list = sorted(set(blastn_subject_list),key=blastn_subject_list.index)
    #combine the two lists without creating duplicates
    blastn_combined_list = blastn_query_list + blastn_subject_list
    blastn_combined_list = sorted(set(blastn_combined_list),key=blastn_combined_list.index)
    #create a temporary object containing both the fasta used to create the blastdb as well as the input fasta
    blast_db_fasta = "blastnDB/" + blast_database + ".fasta"
    tmp_combined_fasta = open(original_fasta).read() + "\n" + open(fasta_input,"r").read()
    blastn_query_list = []
    blastn_subject_list = []
    for line in blastn_result:
        line = str(line).replace("\n","").split(",")
        blastn_query = line[0]
        blastn_subject = line[1]
        blastn_query_list.append(blastn_query)
        blastn_subject_list.append(blastn_subject)
    #the subject list might contain duplicates, which we will remove
    #just to make sure that nothing slips throught, do the same with the query list
    blastn_query_list = sorted(set(blastn_query_list),key=blastn_query_list.index)
    blastn_subject_list = sorted(set(blastn_subject_list),key=blastn_subject_list.index)
    #combine the two lists without creating duplicates
    blastn_combined_list = blastn_query_list + blastn_subject_list
    blastn_combined_list = sorted(set(blastn_combined_list),key=blastn_combined_list.index)
    #create a temporary object containing both the fasta used to create the blastdb as well as the input fasta
    blast_db_fasta = original_fasta
    tmp_combined_fasta = open(blast_db_fasta,"r").read() + "\n" + open(fasta_input,"r").read()
    #create a fasta dictionary
    fasta_chunks = str(tmp_combined_fasta).split(">")
    sequence_dict = {}
    fasta_header_list = []
    for chunk in fasta_chunks:
        if len(chunk) > 1: #for some reason the first item in the fasta_chunks is empty
            chunk_split = str(chunk).split("\n",1) #only split using the first occurance of \n
            header = str(chunk_split[0]).replace("\r","").replace("\n","") #the first item in the list is the header
            #remove any spaces from the header and replace with _
            sequence = chunk_split[1] #the second item in the list is the sequence, might contain \n.  We only did one split for chunk split
            #remove any \n from the seqeunce
            sequence = str(sequence).replace("\n","")
            sequence_dict[header] = sequence
    #now we have a sequence dictionary containing all past and present sequences,
    #this might contain duplicates, but our lists does not so no matter
    #create a fasta file needed for the phylogenetics
    fasta_for_phylo_path = result_dir + project_name + "-blast_hits_added.fasta"
    fasta_for_phylo = open(fasta_for_phylo_path,"w")
    for entry in blastn_combined_list:
        fasta_header = ">" + entry + "\n"
        sequence = sequence_dict[entry] + "\n"
        fasta_for_phylo.write(fasta_header)
        fasta_for_phylo.write(sequence)
    fasta_for_phylo.close()



#main program          
if __name__ == "__main__":
          #blast
        try:
            htmlTop()
            
            #print information out to the webpage
            print '''            
            <h3>  &#8658; Update blastn database</h3>
            
            <h3>  &#8658; Search blastn database for the specified number of best matches</h3>
            
            <h3>Multiple alignment of submitted and retrieved sequences using MAFFT</h3>
            
            <h3>Cleaning of alignment with trimal</h3>
            
            <h3>Calculating and rendering the distance matrix using R</h3>
            
            <h3>Phylogenetic inference using fasttree</h3>
            
            <h3>Rendering of tree using ete3</h3>
            
            '''

            #time the blast function, thus blast, retrieval and duplicate removal
            #t = Timer(blast)
            #blastTime = t.timeit(1)
            #projectName = log_dict["project_name"]
            #line = projectName +  ",blast," + str(blastTime) + "\n"
            #timeFile = open("timeFile.csv","a")
            #timeFile.write(line)
            #timeFile.close()
            
            blast()
            
            
            #now redirect to blast mafft            
            target_url = '<meta http-equiv="refresh" content="0;url=http://' + ip + '/cgi-bin/phylophile/mafft.py" /> '
            print target_url
            
            htmlTail()
              
        except:
            cgi.print_exception()
          
          

              
