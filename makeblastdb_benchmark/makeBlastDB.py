#!/usr/bin/env python3

import random
import os
import time
from timeit import Timer
import subprocess

inputFile = open("cleaned_hiv-db.fasta", "r").read()
resultsFile = open("makeBlastTimes.txt", "a")

inputFile = inputFile.split(">")

#kList = (10, 100, 500, 1000, 5000, 10000)
kList1 =  random.sample(range(100, 500), 50)
kList2 =  random.sample(range(500, 1000), 50)
kList3 =  random.sample(range(1000, 5000), 50)
kList4 =  random.sample(range(5000, 8000), 50)
kList5 =  random.sample(range(8000, 10000), 50)
kList = kList1 + kList2 + kList3 + kList4 + kList5

def makeDb():
    samples = random.sample(inputFile, k)
    tmpFasta = open("tmp.Fasta", "w")
    for sample in samples:
        tmpFasta.write(">{}".format(sample))
    tmpFasta.close()
    makeDbCmd = "./makeblastdb -in tmp.Fasta -dbtype nucl -title tmp -out tmp"
    subprocess.call(makeDbCmd, shell = True)


for k in kList:
    #makeDb()

    t = Timer(makeDb)
    makeDbTime = t.timeit(1)
    resultsFile.write("{}, {}\n".format(k, makeDbTime))
    
resultsFile.close()

    

    



