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
from collections import OrderedDict, Counter
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


form = cgi.FieldStorage()
db_name = form.getvalue("Name")


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
                    <title>Making DB</title>
                </head>
                <body>""")
def htmlTail():
          print("""</body>
            </html>""")



def save_uploaded_file():
    #remove all old blastnDB files
    os.popen("rm blastnDB/*")
    # Get filename here.
    fileitem = form['filename']
    # Test if the file was uploaded
    if fileitem.filename:
        # strip leading path from file name to avoid 
        # directory traversal attacks
        fn = os.path.basename(fileitem.filename)
        open('/var/www/cgi-bin/phylophile/blastnDB/' + str(fn).replace(" ","_"), 'wb').write(fileitem.file.read()) #use replace to remove white spaces

        message = 'The file ' + fn + ' was uploaded successfully'
   
    else:
        message = 'No file was uploaded'
    print message



def parse_fasta_makeblastdb():
    
    fileitem = form['filename']
    fn = "blastnDB/" + os.path.basename(fileitem.filename)
    fasta_for_db = open(fn,"r").read()

    #clean the file, from tabs, \r and others
    fasta_for_db = str(fasta_for_db).replace("\t","\n").replace("\r","\n").replace("\n\n","\n")
    #split the fasta file into chunks
    fasta_for_db_chunks_list = str(fasta_for_db).split(">")
    
    #the first entry in this list is empty, remove it
    fasta_for_db_chunks_list = fasta_for_db_chunks_list[1:]
    #parse the fasta chunk list and split it into heades and sequences
    db_header_list = []
    db_fasta_dict = {}
    for item in fasta_for_db_chunks_list:
        item = str(item).split("\n",1)
        db_header_list.append(item[0])
        db_fasta_dict[item[0]] = item[1]
    
    duplicate_header_list = [item for item, count in Counter(db_header_list).items() if count > 1]
    
    if len(duplicate_header_list) == 0:
        print "Found no duplicates fasta header names in your input file.<p>"
        print "Your blastn database was created<p>"
    if len(duplicate_header_list) > 0:
        print "<p>Found duplicate fasta header names, if these were repeats, please include a date or number with the repeat<p>"
        print "ERROR:  Your blast database was not formatted <p>"
        print "Below are the culprits:<p>"
        for item in duplicate_header_list:
            print item
            print "<p>"
        exit()
    
    #everything is fine and we can continue
    tmp_fasta = "tmp_fasta.fasta"
    blast_input_fasta = open(tmp_fasta,"w") #this will be overwritten if it exists
    
    #now we can itterate trhough the fasta header list and use the iterations to look up the sequence in the dictionary
    for header in db_header_list:
        sequence = db_fasta_dict[header]
        blast_input_fasta.write(">")
        blast_input_fasta.write(header)
        blast_input_fasta.write("\n") #a new line
        blast_input_fasta.write(sequence) 
        blast_input_fasta.write("\n") #a new line
    
    blast_input_fasta.close()
    #rm *.nin *.nsq *.nhr
    makeblastdb_cmd = "./makeblastdb -dbtype nucl -in " + tmp_fasta + " -input_type fasta -title " + db_name + " -out blastnDB/" + db_name
        
    os.popen(makeblastdb_cmd)
    
    conf = open("blastnDB/config.conf","w")
    conf.write("db_name=")
    conf.write(db_name)
    conf.write("\n")
    conf.write("fasta=")
    conf.write(fn)
    conf.close()
    #overrwite the original faste with the tmp_faste which has been unified
    shutil.move(tmp_fasta,fn)




#main program          
if __name__ == "__main__":
        try:
                       
            htmlTop()
            
            save_uploaded_file()

            parse_fasta_makeblastdb()
            
            print '<form>'
            target_url = '''<input type=button onClick="parent.location='http://''' + ip + '''/index.html'" value='Return to input'>'''
            print target_url       
                    
            print  '</form>'
            
            htmlTail()
            
        except:
            cgi.print_exception()





              
