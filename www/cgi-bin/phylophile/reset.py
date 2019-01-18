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
                    <title>Initialize</title>
                </head>
                <body>""")
def htmlTail():
          print("""</body>
            </html>""")




#main program          
if __name__ == "__main__":
        try:
            
            htmlTop()
            
            reset_cmd = "rm archive/* db/* input/* output/* blastnDB/*"
            os.popen(reset_cmd)
            #print information out to the webpage
            print '''
             
            <h3>Your server data has been deleted. </h3>          
            '''
            target_url = '''<input type=button onClick="parent.location='http://''' + ip + '''/index.html'" value='Return to input'>'''
            print target_url
            
            
            
            htmlTail()
        except:
            cgi.print_exception()





              
