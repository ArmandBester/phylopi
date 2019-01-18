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
import cgitb
cgitb.enable(format='text')
import os
import zipfile
from os.path import basename
import shutil
import sqlite3
from datetime import date, datetime
import string
import random


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

result_dir = log_dict["result_dir"]
project_name = log_dict["project_name"]
project_file_path = result_dir  + project_name
     
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
                    <title>Results</title>
                    
                </head>
                <body>""")
def htmlTail():
          print("""</body>
            </html>""")

def random_postfix(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))
postfix = random_postfix()

#def zip_results():
    
    
def zip_results():
    zip_archive = "output/" + log_dict["project_name"] + "_"+ postfix +  ".zip"
    zip_input = "output/" + log_dict["project_name"]
    zipper(zip_input, zip_archive)
    #copy the output zip file to the archive dir
    src = zip_archive
    dst = str(zip_archive).replace("output/","archive/")
    shutil.copyfile(src,dst)

def zipper(dir, zip_file):
    zip = zipfile.ZipFile(zip_file, 'w', compression=zipfile.ZIP_DEFLATED)
    root_len = len(os.path.abspath(dir))
    for root, dirs, files in os.walk(dir):
        archive_root = os.path.abspath(root)[root_len:]
        for f in files:
            fullpath = os.path.join(root, f)
            archive_name = os.path.join(archive_root, f)
            #print f
            zip.write(fullpath, archive_name, zipfile.ZIP_DEFLATED)
    zip.close()
    return zip_file
    
    
    '''
    z = zipfile.ZipFile(zip_archive,"w")
    folder_contents  = os.listdir("output/" + log_dict["project_name"])
    for file_item in folder_contents:
        file_path = "output/" + log_dict["project_name"] + "/" + file_item
        z.write(file_path,basename(file_path))
    z.close()
    '''
    
    

def add_to_db():
    #if the db does not excist create one ant the table, else pass
    if os.path.isfile("db/phyloStore.db"):
        pass
    else:
        db = sqlite3.connect("db/phyloStore.db")
        conn = db.cursor()
        conn.execute('''CREATE TABLE results
                     (id INTEGER PRIMARY KEY not null,
                     name TEXT,
                     pr_samples TEXT,
                     sec_samples TEXT,
                     reps TEXT,
                     date TEXT,
                     file_loc TEXT);''')
        
        
    #add result entries
    #find the last id in the DB so that we can increment for the next id
    conn = sqlite3.connect("db/phyloStore.db")
    db_entry_name = log_dict["project_name"]    
    db_date = date.today()   
    result_file_location = "archive/" + log_dict["project_name"] + "_" + postfix + ".zip"    
    reps = log_dict["reps"]
    
    #find the primary samples, this can be found in the input fasta
    input_fasta = open(log_dict["fasta_input"],"r").readlines()
    fasta_input_list = []  #these are the primary samples
    for line in input_fasta:
        if str(line).find(">") > -1: #if fasta header is found
            fasta_input_list.append(str(line).replace(">","").replace("\n","").replace("\r","").replace(" ",""))
    print "<p>"
    print "Input: "
    print len(fasta_input_list)
    print " "
    print fasta_input_list
    
    #find the seconday or blast hit samples, this can be found by parsing the blastn file in results
    blast_hits = open(log_dict["result_dir"] + log_dict["project_name"] + ".blastn","r").readlines()
    blast_hit_list = []
    for line in blast_hits:
        line = str(line).split(",")
        hit = line[1]
        blast_hit_list.append(hit)
    print "<p>"    
    blast_hit_list = list(set(blast_hit_list))
    print "BLAST Hits: "
    print len(blast_hit_list)
    print " "
    print blast_hit_list
    blast_hits_minus_Q = set(blast_hit_list) - set(fasta_input_list)
    blast_hits_minus_Q = list(blast_hits_minus_Q)
    print "<p>"
    print "BLAST hits without query sequences: "
    print len(blast_hits_minus_Q)
    print " "
    print blast_hits_minus_Q
    
    db_entry = (db_entry_name, str(fasta_input_list).replace("[","").replace("]","").replace("'",""),\
                str(blast_hits_minus_Q).replace("[","").replace("]","").replace("'",""), reps, str(db_date), result_file_location)
    
    conn.execute('INSERT INTO results(name, pr_samples, sec_samples,reps, date, file_loc) VALUES (?,?,?,?,?,?)', db_entry)
    conn.commit()
    conn.close()
    
#main program
if __name__ == "__main__":
        try:
            htmlTop()

            zip_results()
            print "<h2>Phylogenetic Results</h2>"
            # print warning fir the mapping sanity check
            q_mapping = open("q_mapping", "r").readlines()
            q_mapping = q_mapping[1:]
            if len(q_mapping) > 0:
                print "<b>Warning: the following sequences seems wrong for this analysis</b> <p>"
                for line in q_mapping:
                    line = str(line).split(",")[1]
                    line = '<font color="red">' + line + '</font>'
                    print line
                    print "<p>"
                print "Check your alignment and your queries mapped against the HXB2 reference, <p>"
                print "using the links below. <p>"
            print "Click on the link below to download your results"
            print "<p>"
            print('''<a href="output/''' + log_dict["project_name"] + "_" + postfix +  ".zip" + '''">Download</a>''')
            print "<p>"
            target_url = '''<input type=button onClick="parent.location='http://''' + ip + '''/index.html'" value='Return to input'>'''
            print target_url
            print "<p>"
            #create links for the alignment and dictance matrix
            dist_matrix = '<a href="http://' + ip + '/cgi-bin/phylophile/' + project_file_path + '_heatmap.html"' + ' target="_blank">View Distance Matrix</a> '
            html_aln = '<a href="http://' + ip + '/cgi-bin/phylophile/' + project_file_path + '-mafft.aln_trimal.html"' + ' target="_blank">View MAFFT alingment</a> '
            q_map_link = '<a href="http://' + ip + '/cgi-bin/phylophile/' + result_dir + '/q_mapping.svg"' + ' target="_blank">View the mapped queries</a> '
            print dist_matrix
            print "   "
            print html_aln
            print "   "
            print q_map_link
            print "<p>"

            


            png_path = "output/" + log_dict["project_name"]  + "/" + log_dict["project_name"] + ".svg"
            
            print "<img src=" + png_path +"/>"
            
            add_to_db()
            
            #print stuff about what has and has not been added to the blast database
            print "<p>"
            print "Samples added to blast database: "
            print log_dict["Samples added "]
            print "<p>"
            print "Samples NOT added to blast database: "
            print log_dict["Samples not added "]    
            htmlTail()
        except:
            cgi.print_exception()




