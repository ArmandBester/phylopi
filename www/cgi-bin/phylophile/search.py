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
import sqlite3

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
                    <h2>Previous Results</h2>
                </head>
                <body>""")
def htmlTail():
          print("""</body>
            </html>""")


def search():
    form = cgi.FieldStorage()
    name = form.getvalue("Name")
    pr_sample = form.getvalue("PriSample")
    sec_sample = form.getvalue("SecSample")
    #return name, pr_sample, sec_sample

    db = sqlite3.connect("db/phyloStore.db")
    conn = db.cursor()
    entry_list = [] #create a list of row tuples
    #if all fields are empty return everything
    query = "select * from results where pr_samples like '%"+ pr_sample + "%' and name like '%" + name + "%' and sec_samples like '%" + sec_sample + "%'"

    for row in conn.execute(query): #("select * from results"):
        entry_list.append(row)
    for entry in entry_list:
        ID = entry[0]
        NAME = entry[1]
        REPS = entry[4]
        DATE = entry[5]
        RESULT = entry[6]
        
        print "Result for: "
        print "<b>"        
        print NAME
        print "</b>"
        print "<br>"
        print "# Hits &emsp; Date"
        print "<br>"
        print REPS
        print "&emsp;&emsp;&emsp;"
        print DATE
        print "&emsp;"
        print ('''<a href="http://''' + ip + '''/cgi-bin/phylophile/''' + RESULT + '''">Download</a>''') 
        print "<p>"
        print "<p>"
    print "<p>"
    print "<p>"



#main program
if __name__ == "__main__":
        try:
            htmlTop()
            
            print '<form>'
            target_url1 = '''<input type=button onClick="parent.location='http://''' + ip + '''/index.html'" value='Return to input'>'''
            target_url2 = '''<input type=button onClick="parent.location='http://''' + ip + '''/search.html'" value='Return to search'>'''
            print target_url1
            print target_url2
            print '</form>'
            
              
            search()
            
            
              
            htmlTail()
        except:
            cgi.print_exception()
