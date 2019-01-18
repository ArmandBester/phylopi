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
import time
import sys
import multiprocessing
import shutil


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
    
    
    

os.environ['PYTHON_EGG_CACHE'] = '/home/pi/python_eggs/'


#create a switch for adding or not adding unique new uploads to the blast database
add_new_unique_seqs = 1

fn = ''


form = cgi.FieldStorage()

blast_database = ""
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



def save_uploaded_file():
    # Get filename here.
    fileitem = form['filename']
    
    global fn
  

    # Test if the file was uploaded
    if fileitem.filename:
        # strip leading path from file name to avoid 
        # directory traversal attacks
        fn = os.path.basename(fileitem.filename)
        open('/var/www/cgi-bin/phylophile/input/' + str(fn).replace(" ","_"), 'wb').write(fileitem.file.read()) #use replace to remove white spaces

        message = 'The file ' + fn + ' was uploaded successfully'
   
    else:
        message = 'No file was uploaded'
        print message

def cln_output():
    #this is called to remove old file from these directories
    os.popen("rm -r output/")
    
def add_new_to_db_fasta():
    #create a header list of the main fasta file
    main_header_list = []
    
    conf = open("blastnDB/config.conf","r").readlines()
    for line in conf:
        if str(line).find("db_name=") > -1:
            db_name = str(line).split("=")[1].replace("\n","")
        if str(line).find("fasta=") > -1:
            main_fasta = str(line).split("=")[1].replace("\n","")
        
    fasta_for_db = open(main_fasta).readlines()
    
    for line in fasta_for_db:
        if str(line).find(">") > -1:
            main_header_list.append(str(line).replace(">","").replace("\n",""))
            
    #create a header list of the new to be added fasta file
    new_header_list = []
    new_fasta_dict = {}
    
    input_file = "input/" + str(fn).replace(" ","_")
    input_fasta = open(input_file, "r").read()
    input_fasta = str(input_fasta).replace("\t", "\n").replace("\r", "\n").replace("\n\n", "\n")
    input_fasta_chunks = str(input_fasta).split(">")[1:]
    for chunck in input_fasta_chunks:
        chunck = str(chunck).split("\n",1)
        header = chunck[0]
        sequence = chunck[1]
        new_header_list.append(header)
        new_fasta_dict[header] = sequence
    unique_list = []
    non_unique_list = []
    for fasta_header in new_header_list:
        if fasta_header in  main_header_list:
            non_unique_list.append(fasta_header)
        else:
            unique_list.append(fasta_header)
    
    #append the unique items to the main fasta file
    fasta_updating_handle = open(main_fasta, "a")
    fasta_updating_handle.write("\n")
    for item in unique_list:
        new_header = ">" + item + "\n"
        new_sequence = new_fasta_dict[item] + "\n"
        fasta_updating_handle.write(new_header)
        fasta_updating_handle.write(new_sequence)
    fasta_updating_handle.close()
    
    tmp_log_file = open("input/tmp.log","a") #apend to the log file
    tmp_log_file.write("Samples added = ")
    tmp_log_file.write(str(unique_list).replace("[","").replace("]",""))
    tmp_log_file.write("\n")
    tmp_log_file.write("Samples not added = ")
    tmp_log_file.write(str(non_unique_list).replace("[","").replace("]",""))
    tmp_log_file.write("\n")
    
    #REFORMAT BLASTDB
    #only required if new entries have been added
    if len(unique_list) > 0:
        #print "Reformating BLASTn databse " + db_name
        
        makeblastdb_cmd = "./makeblastdb -dbtype nucl -in " + main_fasta + " -input_type fasta -title " + db_name + " -out blastnDB/" + db_name
        #print makeblastdb_cmd
        tmp_log_file.write(makeblastdb_cmd)
        tmp_log_file.write("\n")
        os.popen(makeblastdb_cmd)
        
    
    
    
    
    
    
        


def db_checks():

    #print "   \nDOING SOME CHECKS <p>"
    #do some checks
    #see if config.conf exists in blastnDB
    if os.path.isfile("blastnDB/config.conf") == True:
        pass
        #print "\nFound the config file in directory, blastnDB.<p>"
    else:
        print "\nERROR: Could not find config.conf in directory, blastnDB<p>"
        print "Did you remember to setup your blast database first?<p>"
        print '<form>'
        target_url_return = '''<input type=button onClick="parent.location='http://''' + ip + '''/index.html'" value='Return to input'>'''
        print target_url_return
        print '''<p>
                </form>'''
        quit()

    ###check if the the name of the DB in the config file exists
    global blast_database
    
    conf = open("blastnDB/config.conf","r").readlines()
    for line in conf:
        if str(line).find("db_name") > -1:
            blast_database = str(line).split("=")[1]
            blast_db = blast_database
    nhr = "blastnDB/" + blast_database + ".nhr"
    nin = "blastnDB/" + blast_database + ".nin"
    nsq = "blastnDB/" + blast_database + ".nsq"
        



#main program          
if __name__ == "__main__":
        try:
            
            htmlTop()
            
            #print information out to the webpage
            print '''
             
            <h3> &#8658; Update blastn database </h3>
            
            <h3>Search blastn database for the specified number of best matches </h3>
            
            <h3>Multiple alignment of submitted and retrieved sequences using MAFFT</h3>
            
            <h3>Cleaning of alignment with trimal</h3>
            
            <h3>Calculating and rendering the distance matrix using R</h3>
            
            <h3>Phylogenetic inference using fasttree</h3>
            
            <h3>Rendering of tree using ete3</h3>
            
            '''
            
            
            
            
            
            cln_output()  
            #before we process a new input fasta file, we first need to make sure the dir is empty
            fasta_rm_cmd = "rm /var/www/cgi-bin/phylophile/input/*"
            os.popen(fasta_rm_cmd)
            #create a log file to keep track
            tmp_log = open("/var/www/cgi-bin/phylophile/input/tmp.log","a")
              
            #now we can add the new fasta file to the input dir for processing
            save_uploaded_file()
            
            
            add_new_to_db_fasta()  
            #lets do some blast db checks, just in case              
            db_checks()              
              


            #input fasta check
            input_file_list = os.listdir("input")
            fasta_list = []
            for entry in input_file_list:
                if str(entry).find(".fasta") > -1 or str(entry).find(".fas") > -1:
                    fasta_list.append(entry)

            project_name = str(fasta_list[0]).replace(".fasta","").replace(".fas","")
            fasta_input = "input/" + str(fasta_list[0])
            tmp_log.write("fasta_input=")
            tmp_log.write(fasta_input+"\n")
            tmp_log.write("project_name=")
            tmp_log.write(project_name+"\n")
            result_dir = "output/" + project_name + "/"
            tmp_log.write("result_dir=")
            tmp_log.write(result_dir+"\n")
            tmp_log.write("blast_database=")
            tmp_log.write(blast_database+"\n")
            #print "Your project name is: " + project_name + "<p>"
            out_path = "output/" + project_name + "/"
            
            out_path = "output/" + project_name + "/"
            if os.path.isdir(out_path) == False:
                os.makedirs(out_path)

            
            
            #look for the radio value on how many blast reps are wanted
            reps = form.getvalue("hits")
            #print "Number of best blast matches to include: "
            #print reps
            tmp_log.write("reps=")
            tmp_log.write(reps+"\n")
            tmp_log.close()
            
            
                       
            #target_url = '<meta http-equiv="refresh" content="0;url=http://' + ip + '/cgi-bin/phylophile/blast.py" /> '
            target_url = '<meta http-equiv="refresh" content="0;url=http://' + ip + '/cgi-bin/phylophile/sanity.py" /> '
            print target_url


            htmlTail()
        except:
            cgi.print_exception()





              
