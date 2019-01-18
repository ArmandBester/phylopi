# PhyloPi:  

PhyloPi is an automated phylogenetics pipeline designed for the drug resistance testing facility.  It uses BLAST to update a rolling database with each submission as well as retrieving most similar sequences to include in the analysis.

Our manuscript is available at https://www.biorxiv.org/content/early/2018/07/13/367946

                   
##Installation instructions

###On a Raspberry Pi 3

The image for PhyloPi can be downloaded from 

http://scholar.ufs.ac.za:8080/xmlui/handle/11660/7638


We recommend installing/flashing this image on a SD-card with at least 16GB of space. You will need a SD card reader and of course a Raspberry Pi 3. The easiest and most generic method for flashing the image to the SD card is to use etcher, available for Windows, MacOS and Linux.

Download the PhyloPi_release_ddmmyyyy.zip file from the provided link. Uncompress this zip archive resulting image is a bzip2 compressed archive and conveniently etcher can use this file as input without the need to first decompress it.

After flashing, insert the SD card into the Raspberry Pi 3 and power it up using a suitable power supply. You don’t have to connect a screen to the Pi, our pipeline was designed to be headless. After about a minute you should see a WiFi hotspot broadcast as PhyloPi, connect to this with the default password (phylophile). Next you need to expand the file system of the operating system. On Linux or MacOS use any terminal of your choice, for Windows we suggest using Putty.

For more information on connecting to the Pi via ssh, please see:

Windows:
https://www.raspberrypi.org/documentation/remote-access/ssh/windows.md

Linux and MacOS
https://www.raspberrypi.org/documentation/remote-access/ssh/unix.md

Remember the IP address of the Pi is 192.168.1.1, you will need that for the connection.

Username: pi

Password: phylophile

After establishing a ssh connection, type the command

sudo raspi-config

Use the arrow keys to navigate to ‘Advanced Options’ and select ‘Expand File System’. Select OK and ‘Finish’. At the prompt, ‘Do you want to reboot now?’ select ‘Yes’. Close your SSH session.

The phylogenetic pipeline is now ready to be set up. Open a browser (we suggest Firefox, Google Chrome or Safari) and navigate to http://192.168.1.1. If this is the first time you use PhyloPi you will first need to setup a blast database. Click on the button ‘BLASTdb conf’, browse to the fasta file you want to use as a BLAST database. The fasta headers should not contain any commas, semi-colons or spaces. Good practice is to use hyphens, underscores or periods as seperators in fasta headers. Give your database a name. Click on Format DB. You are now ready to do some analysis.

If for any reason you want to reset your PhyloPi to ‘factory state’ follow the link:

http://192.168.1.1/reset.html

and follow the instructions. This will erase all data!

Links:
Raspberry Pi3 info: https://www.raspberrypi.org/products/raspberry-pi-3-model-b/

etcher: https://etcher.io/ A multiplatform aplication to flash the PhyloPi image onto a SD card.

Putty: https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html

###Installing on another platform

If you want to install our pipeline on another platform, let say a virtual machine which you can run on your Windows PC or Mac, please contact us for an image.

Alternatively, you can roll your own.  The source code is available at

https://github.com/ArmandBester/phylopi

Be sure to have binaries compile for your architecture, the ones from github is for ARM. 


If you need any help, have any suggestions or would like to contribute, contact us on:

besterpa@ufs.ac.za
or
besterpa.sci@gmail.com


<br>
<br>
<br>
<br>
<br>
<br>


