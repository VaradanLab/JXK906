"""
Created on Thu Aug 10 21:01:43 2017

@author: jaspreet
"""
import os
import sys
import subprocess as sub
#import pandas as pd
#import numpy as np
import time
import re
import datetime as dt
from os.path import join, getsize

# function for base recalibration after add and replace read groups
 
path = "/mnt/pan/Data16/jxk906/PROJECTS/TNBC/AA/"
files_in_directory= os.listdir(path)
print("____Number of files after Add_replace_grps____")
print([".bam files",len([f for f in files_in_directory if(f.endswith(".bam"))]),
       ".bai files",len([f for f in files_in_directory if(f.endswith(".bai"))])])
print("____Proceeding with Base recalibration____") ## need to add a bit of time
time.sleep(1)
print("Slurm options to be default y/n:")
slurm_options= input()
if slurm_options!="Y":
    print("hours")
    hr=input()
    print("minutes")
    mins=input()
    mints=str(mins).zfill(2)
    print("seconds")
    sec=input()
    secs=str(sec).zfill(2)
    t=":".join([str(hr), str(mints), str(secs)])
    print("memory required")
    mem=input()
    print("nodes required")
    nodes=input()
    print("no of cores")
    cpus=input()
out_dir = "/mnt/pan/Data16/jxk906/PROJECTS/TNBC/BASERECALIBRATED_FILES"
time.sleep(1)

for f in files_in_directory:
    # CHUNK 1 getting the size of files initial start
   if (f.endswith(".bam")):
     i=re.split("\.",f)[0]
     print(i)
     size=((os.path.getsize(os.path.join(path+f)))/1000000) ## giving size in bytes 
     print("____AddOrReplacedGrps file size____")
     time.sleep(1)
     print("%s"%f+"____"+str(size)+"MB")
     time.sleep(1)
    # CHUNK 2 making SLURM SCRIPTS
     print("____Building SLURM SCRIPTS now____")
     time.sleep(1)
     print("____Writing Slurm for %s_____"%i)
     time.sleep(1)
     start_time=time.time()
     slurm=open("/mnt/pan/Data16/jxk906/PROJECTS/TNBC/BR_scripts/%s.slurm"%i,"w")
     if slurm_options=="Y":
        options=["#!/bin/bash","\n","#SBATCH --mail-user=jxk906@case.edu","\n",
             "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%i,"\n","#SBATCH --time=80:00:00",
             "\n","#SBATCH --nodes=1","\n","#SBATCH --cpus-per-task=1","\n","#SBATCH --output=%s.out"%i,
             "\n","#SBATCH --mem=8gb","\n","#SBATCH -A lxh259","\n", "module load gatk","\n"]
     else:
         options=["#!/bin/bash","\n","#SBATCH --mail-user=jxk906@case.edu","\n",
             "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%i,"\n","#SBATCH --time=%s"%t,
             "\n","#SBATCH --nodes=%s"%nodes,"\n","#SBATCH --cpus-per-task=%s"%cpus,"\n","#SBATCH --output=%s.out"%i,
             "\n","#SBATCH --mem=%sgb"%mem,"\n","#SBATCH -A lxh259","\n", "module load gatk","\n"]
     calib_grp = i+".recal.grp"
     grp_file="/".join([out_dir,calib_grp]) # out_dir
     for j in range(len(options)):
         #print(options[j])
         slurm.write(options[j])
     slurm.write("\n")
     slurm.write("/mnt/pan/Data16/Tools/jdk1.8.0_102/bin/java -Xmx8g -jar $GATK -T BaseRecalibrator -R /mnt/pan/Data16/jxk906/REFERENCE/hg19.fa -I %s -knownSites /mnt/pan/Data16/jxk906/split_vcf/Merged_dbsnp.vcf -o %s"%(os.path.join(path+f),grp_file))
     slurm.write("\n")
     slurm.write("wait")
     time.sleep(1)
   else:
      continue
end_time=time.time()
print("___Elapsed time %isec ____"%(end_time-start_time))
slurm.close()
# for execution of slurm scripts and looking at output of slurm scripts and the size of actual files
#1. path where slurm scripts are
#2. path where my output od scripts are (which is same as my scrpit path)
#3. path where my sample files are being written

def knowJobStatus(uid):
    '''
    Name:           Run slurm scripts
    Description:        uses sbatch
    Keyword argument:
    uid -- user id 
    '''
    p=sub.Popen(['squeue','-u','%s'%uid], stdout=sub.PIPE)
    soutput,sinput=p.communicate()
    return soutput

script_dir="/mnt/pan/Data16/jxk906/PROJECTS/TNBC/BR_scripts"
slurm_list=[f for f in os.listdir(script_dir) if f.endswith(".slurm")]
print("Total scripts %s"%len(slurm_list))
sample_dir= "/mnt/pan/Data16/jxk906/PROJECTS/TNBC/BASERECALIBRATED_FILES/" ## output directory
for i in slurm_list:
    s1= "/".join([script_dir,i])
    sub.call("sbatch %s"%s1, shell=True)
    print(knowJobStatus("jxk906"))
  
## NEED TO LOOK INTO HOW TO OPEN THE OUTPUT SCRIPTS and the actual filesize at times
## for this i want a max limit in which my code will open all files one by one
## andshould print last lines of .out files after every 15 mins of current time
## but once i see all files output as successfull loop should break
## if there is an error in any file stop the loop

slurm_out_list=[f for f in os.listdir(script_dir) if f.endswith(".out")]
d=(dt.timedelta(hours=24).days)*24 ## subject of input or alteration
current_time=time.time()
expected_time=current_time+(d*3600)
print("Total_out_scripts = %s"%len(slurm_out_list)) ## to know how many have been 
#print("no of files you want to open at a time:")
#a_files=input()
def GetLastNLines(n, fileName):
    '''
    Name:           Get LastNLines
    Description:        Gets last n lines using Unix tail
    Output:         returns last n lines of a file
    Keyword argument:
    n -- number of last lines to return
    filename -- Name of the file you need to tail into
    '''
    p=sub.Popen(['tail','-n',str(n),fileName], stdout=sub.PIPE)
    soutput,sinput=p.communicate()
    return soutput
def forOutputPrint(slurm_out_list, sample_dir):
        for i in slurm_out_list:
            s2= "/".join([script_dir,i])
            split = re.split("\.",i)[0]+".recal.grp"
            grp_file_size=((os.path.getsize(os.path.join(sample_dir+split)))/1000000)
            print(GetLastNLines(10,"%s"%s2))
            print("\n")
            print("%s :: %sMB"%(re.split("\.", i)[0]+".recal.grp", grp_file_size))
            time.sleep(5)

from threading import Timer, Thread

def call_at_interval(time, callback, args):
    while True:
        timer = Timer(time, callback, args=args)
        timer.start()
        print("continue (Y/N)")
        inp=input()
        if inp!="Y":
            timer.cancel()
        else:
            timer.join()

def setInterval(time, callback, *args):
    Thread(target=call_at_interval, args=(time, callback, args)).start()

setInterval(60,forOutputPrint,slurm_out_list,sample_dir)

            
    
