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
                    
                    </style>
                    
                    <meta charset="utf-8"/>
                    <title>Results</title>
                    <h2>Find sequences in DB</h2>
                </head>
                <body>""")
def htmlTail():
          print("""</body>
            </html>""")
          
#get values from the conf file in bastnDB folder
conf = open("blastnDB/config.conf","r").readlines()
blastnDB = ''
original_fasta = ''
for line in conf:
    if str(line).find("db_name") > -1:
        blastnDB = str(line).split("=")[1]
    if str(line).find("fasta") > -1:
        original_fasta = str(line).split("=")[1]
        
        

fasta_name = "blastnDB/" + str(open("blastnDB/config.conf").read()).replace("\n","") + ".fasta"
fasta_dict = {}
fasta_list = []
def html_parse_db_fasta():
    main_fasta_file = open(original_fasta,"r").readlines()
    fasta_header = ""
    fasta_seq = ""
    for line in main_fasta_file:
        if str(line).find(">") > -1:
            fasta_header = str(line).replace(" ", "_").replace("\n","").replace("\r","").replace(">","") #remove spaces
            fasta_seq = "" #reset the fasta seq to empty, new header is detected, thus new squence, see line below
            fasta_list.append(fasta_header)
            fasta_header
        else:
            fasta_seq = fasta_seq + line #if the fasta sequence goes across multiple lines conacanate it until > is found again
            fasta_dict[fasta_header] =  fasta_seq
            
            


query_list = []
contains_list = []
def find_fasta():
    form = cgi.FieldStorage()
    query = form.getvalue("samples")
    contains = form.getvalue("contains")
    query = str(query).replace(";","").replace(" ","").replace("'","")
    query_list = str(query).split(",")
    if contains == None:
        for entry in query_list:
            if fasta_dict.has_key(entry):
                print ">" + entry  + "<br />"
                print fasta_dict[entry]  + "<br />"
            else:
                print  ">" + entry  + "<br />"
                print "not found<br />"
                
    if contains == "Contains":
        result_list = []
        for entry in query_list:
            for sample in fasta_list:
                if str(sample).find(entry) > -1:
                    result_list.append(sample)
        for result in result_list:
            print ">" + result + "<br />"
            print fasta_dict[result] +"<br />"

    print "<br />"

#main program
if __name__ == "__main__":
        try:
            htmlTop()
            
            print '''<form>
                    <input type=button onClick="parent.location='http://%s/index.html'" value='Return to input'>
                    <input type=button onClick="parent.location='http://%s/find_samples.html'" value='Return to Find'>
                    <p>
                    </form>''' %(ip, ip)
            
            html_parse_db_fasta()
            find_fasta()

            
            
              
            htmlTail()
        except:
            cgi.print_exception()
