#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 01:26:21 2017

@author: jxk906
"""

import os
import sys
import subprocess as sub
import pandas as pd
import numpy as np
import time
import re
import datetime as dt
from os.path import join, getsize
global PYTHON_PATH
print("ENTER : Location of python executable files")
PYTHON_PATH=input()
sys.path.append(os.path.abspath(PYTHON_PATH))
import FUNCTIONS_FOR_PIPELINE as fp

''' CHUNK 00 GETTING THE REQUIRED DIRECTORIES and LOCATION OF TOOLS REQUIRED'''

def AskLocations():
    '''
    Name: Ask Locations
    Description > This function is required for the pipeline (must be called irrespective of any function called)
                > Asks the locations of PROJECT,REF_GENOME,vcf, JAVA, b_mod, BWA, g_mod, GATK, PICARD, FASTQ 
    '''
    global PROJECT, REF_GENOME, vcf, JAVA, b_mod, BWA, g_mod, GATK, PICARD, FASTQ, MAIL, ACCOUNT
    # global variables to be written inside the function
    print("> Input locations for required files and directories [AVOID using '/' in the end of all paths]")
    time.sleep(2)
    print("> ENTER PROJECT LOCATION [all folders for processing will be created at this location]")
    PROJECT = input()
    print(">ENTER LOCATION OF REFERENCE GENOME (hg19/Grch37)")
    REF_GENOME= input()
    print("> ENTER PATH OF dbSNP_VCF FILE:")
    vcf=input()
    print(">ENTER LOCATION TO Java (version-1.8.x)")
    JAVA= input()
    print(">ENTER LOCATION TO PICARD.jar (version-2.x.x)")
    PICARD=input()
    print(">Load Slurm BWA module: (y/n)")
    b_in=input()
    if b_in != "y":
        print(">ENTER LOCATION TO BWA (version-0.7.x)")
        b_mod= "\t"
        BWA=input()
    else:
        print("Slurm Module for BWA (version-0.7.12) Is Loaded")
        b_mod = "module load bwa"
        BWA= "bwa"
    print(">Load GATK module in SLURM: (y/n)")
    gatk_in = input()
    if gatk_in != "y":
        print("> ENTER LOCATION TO GATK.jar (version-3.6 or higher)")
        g_mod = "\t"
        GATK=input()
    else:
        print("> Slurm Module for GATK (version-3.6) Is Loaded")
        g_mod="module load gatk"
        GATK ="$GATK"
    print(">ENTER LOCATION OF FASTQ FILES")
    FASTQ = input()
    print(">ENTER USER MAIL FOR SLURM SCRIPT")
    MAIL = input()
    print(">ENTER ACCOUNT WITH CPU CORES ASSIGNED FOR SLURM JOBS")
    ACCOUNT = input()
    return PROJECT, REF_GENOME, vcf, JAVA, b_mod, BWA, g_mod, GATK, PICARD, FASTQ, MAIL, ACCOUNT

#PROJECT, REF_GENOME, vcf, JAVA, b_mod, BWA, g_mod, GATK, PICARD, FASTQ = AskLocations()

''' CHUNK 01 ALIGNMNET '''
#c+=1 # for process 1
def alignment(PROJECT, FASTQ):
        PROC="Alignment"
        INIT_FASTQ_FILES= FASTQ
        OUT_DIR_ALIGN="/".join([PROJECT,PROC])
        if not os.path.exists(OUT_DIR_ALIGN):
            os.mkdir(OUT_DIR_ALIGN)
        files_in_directory= os.listdir(INIT_FASTQ_FILES)
        slurm_scripts = "/".join([OUT_DIR_ALIGN, "alignment"])
        if not os.path.exists(slurm_scripts):
            os.mkdir(slurm_scripts)
        print("> Number of FASTQ files :  %s"%len([f for f in files_in_directory if(f.endswith(".fastq"))]))
        time.sleep(1)
        print(">Proceeding with ALIGNMENT OF PAIRED READS WITH REFERENCE:")
        time.sleep(1)
        print("> SLURM OPTIONS DEFAULT y/n:")
        slurm_options= input()
        if slurm_options!="y":
            print("hours")
            hr=input()
            print("minutes")
            mins=input()
            mints=str(mins).zfill(2)
            print("seconds")
            sec=input()
            secs=str(sec).zfill(2)
            t=":".join([str(hr), str(mints), str(secs)])
            print("slurm memory required")
            mem=input()
            print("nodes required")
            nodes=input()
            print("no of cores (>=1 or <=8)")
            cpus=input()
            
        time.sleep(1)
        # CHUNK 01 : MAKING SLURM SCRIPTS
        print("> Require file name seperator:")
        sep = input()
        time.sleep(1)
        files_for_match = files_in_directory[1:]
        pairs = {}
        for f in files_in_directory:
            for fm in files_for_match:
                if f.endswith(".fastq") and fm.endswith(".fastq"):
                    p1 = f.split(".")[0].split(sep) ## input seperator in filename (# here in my case it is _)
                    p2 = fm.split(".")[0].split(sep) ## input the seperator in filename
                    if "R1" in p1 and "R2" in p2:
                        i= p1.index("R1"); j = p2.index("R2")
                        if "".join(p1[:i]) == "".join(p2[:j]):
                            key = "_".join(p2[:j])
                            print("Matched pairs %s and %s"%(f,fm))
                            pairs[key]=[f,fm]
                            pairs[key].append([f,fm])
                        #else:
                            #print("pairs not found")
                    elif "R2" in p1 and "R1" in p2:
                        i= p1.index("R2") ; j = p2.index("R1")
                        if "".join(p1[:i]) == "".join(p2[:j]):
                            key = "_".join(p2[:j])
                            print("Matched pairs %s and %s"%(f,fm))
                            pairs[key]=[f,fm]
                            pairs[key].append([f,fm])
                        #else:
                         #   print("pairs not found")
                    else:
                        continue
        for k,v in pairs.items():
             print("> Writing Slurm for %s"%k)
             time.sleep(1)
             slurm=open("%s/%s.slurm"%(slurm_scripts,k),"w")
             if slurm_options=="y":
                options=["#!/bin/bash","\n","#SBATCH --mail-user=%s"%MAIL,"\n",
                     "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%i,"\n","#SBATCH --time=15:00:00",
                     "\n","#SBATCH --nodes=1","\n","#SBATCH --cpus-per-task=1","\n","#SBATCH --output=%s.out"%i,
                     "\n","#SBATCH --mem=8gb","\n","#SBATCH -A %s"%ACCOUNT,"\n","#SBATCH --workdir=%s"%slurm_scripts,"\n"]
             else:
                 options=["#!/bin/bash","\n","#SBATCH --mail-user=%s"%MAIL,"\n",
                     "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%i,"\n","#SBATCH --time=%s"%t,
                     "\n","#SBATCH --nodes=%s"%nodes,"\n","#SBATCH --cpus-per-task=%s"%cpus,"\n","#SBATCH --output=%s.out"%i,
                     "\n","#SBATCH --mem=%sgb"%mem,"\n","#SBATCH -A %s"%ACCOUNT,"\n","#SBATCH --workdir=%s"%slurm_scripts,"\n"]
             out_sam = i+".sam" ### need the extension option to be inputted.
             sam_file="/".join([OUT_DIR_ALIGN,out_sam]) # out_dir for baserecalibrator and -BQSR for Print Reads
             for j in range(len(options)):
                 #print(options[j])
                 slurm.write(options[j])
             slurm.write("\n")
             slurm.write("%s mem -t %s %s %s %s > %s"%(BWA,cpus,REF_GENOME,"/".join([INIT_FASTQ_FILES,v[0]]),"/".join([INIT_FASTQ_FILES,v[1]]),sam_file))
             slurm.write("\n")
             slurm.write("wait")
             time.sleep(1)
        
        slurm.close()
        print(">EXECUTION OF SCRIPTS")
        time.sleep(1)
        # CHUNK 02: EXECUTION OF SLURM SCRIPTS
        slurm_list=[f for f in os.listdir(slurm_scripts) if f.endswith(".slurm")]
        print("Total scripts %s"%len(slurm_list))  ## output directory same as out_dir
        
        for i in slurm_list:
            s1= "/".join([slurm_scripts,i]) 
            sub.call("sbatch %s"%s1, shell=True)
        
        slurm_out_list=[f for f in os.listdir(slurm_scripts) if f.endswith(".out")]
        print("Total_out_scripts = %s"%len(slurm_out_list))
        
        # CHUNK 03:CODES AFTER COMPLETION OF JOBS
        fp.GetErrors(slurm_scripts)
        slurm_script_align = slurm_scripts
        #slurm_out_list_align = slurm_out_list
        return  [slurm_script_align , OUT_DIR_ALIGN, INIT_FASTQ_FILES]

    
''' CHUNK 02 SAMSORT '''
def samsort(PROJECT,OUT_DIR_ALIGN):
            print(">STEP2: SORTING SAM USING SortSam")
            PROC2="Sorting"
            INIT_DIR_SORTING =  OUT_DIR_ALIGN ## output directory of aligned files 
            OUT_DIR_SORTING = "/".join([PROJECT,PROC2]) ## ooutput directory for sorted files will be the directory made using project location and name of the process
            if not os.path.exists(OUT_DIR_SORTING):
               os.mkdir(OUT_DIR_SORTING)
            files_in_directory= os.listdir(INIT_DIR_SORTING) ## getting the list of aligned files 
            slurm_scripts = "/".join([OUT_DIR_SORTING, "sort_scripts"])
            if not os.path.exists(slurm_scripts):
                os.mkdir(slurm_scripts)
            print("> Number of SAM files after ALIGNMENT:")
            print([".sam files",len([f for f in files_in_directory if(f.endswith(".sam") or f.endswith(".bam"))])])
            time.sleep(1)
            print("> Proceeding with SORTING STEP") ## need to add a bit of time
            time.sleep(1)
            print("> ENTER PICARD TOOL SHOWN: SortSam") ## need to see whether to add it or not cz unlike gatk toolsonly one is required
            tool=input()
            print("java memory required")
            jav = input()
            print("> SLURM OPTIONS DEFAULT y/n:")
            slurm_options= input()
            if slurm_options!="y":
                print("hours")
                hr=input()
                print("minutes")
                mins=input()
                mints=str(mins).zfill(2)
                print("seconds")
                sec=input()
                secs=str(sec).zfill(2)
                t=":".join([str(hr), str(mints), str(secs)])
                print("slurm memory required")
                mem=input()
                print("nodes required")
                nodes=input()
                print("no of cores (>=1 or <=8)")
                cpus=input()
                
            time.sleep(1)
            # CHUNK 01 : MAKING SLURM SCRIPTS
            print("> Building SLURM SCRIPTS")
            time.sleep(1)
            
            for f in files_in_directory:
               if f.endswith(".sam") or f.endswith(".bam"):
                 i=re.split("\.",f)[0]
                 print("> Writing Slurm for %s"%i)
                 time.sleep(1)
                 slurm=open("%s/%s.slurm"%(slurm_scripts,i),"w")
                 if slurm_options=="y":
                    options=["#!/bin/bash","\n","#SBATCH --mail-user=%s"%MAIL,"\n",
                         "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%i,"\n","#SBATCH --time=15:00:00",
                         "\n","#SBATCH --nodes=1","\n","#SBATCH --cpus-per-task=1","\n","#SBATCH --output=%s.out"%i,
                         "\n","#SBATCH --mem=8gb","\n","#SBATCH -A %s"%ACCOUNT,"\n","#SBATCH --workdir=%s"%slurm_scripts,"\n"]
                 else:
                     options=["#!/bin/bash","\n","#SBATCH --mail-user=%s"%MAIL,"\n",
                         "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%i,"\n","#SBATCH --time=%s"%t,
                         "\n","#SBATCH --nodes=%s"%nodes,"\n","#SBATCH --cpus-per-task=%s"%cpus,"\n","#SBATCH --output=%s.out"%i,
                         "\n","#SBATCH --mem=%sgb"%mem,"\n","#SBATCH -A %s"%ACCOUNT,"\n","#SBATCH --workdir=%s"%slurm_scripts,"\n"]
                 out_sam = i+".bam" ### need the extension option to be inputted.
                 sam_file="/".join([OUT_DIR_SORTING,out_sam]) # out_dir for baserecalibrator and -BQSR for Print Reads
                 for j in range(len(options)):
                     #print(options[j])
                     slurm.write(options[j])
                 slurm.write("\n")
                 slurm.write("%s -Xmx%sg -jar %s %s INPUT=%s OUTPUT=%s SORT_ORDER=coordinate CREATE_INDEX=TRUE VALIDATION_STRINGENCY=LENIENT"%(JAVA,jav,PICARD,tool,"/".join([INIT_DIR_SORTING,f]),sam_file))
                 slurm.write("\n")
                 slurm.write("wait")
                 time.sleep(1)
               else:
                  continue
            
            slurm.close()
            print(">EXECUTION OF SCRIPTS")
            time.sleep(1)
            # CHUNK 02: EXECUTION OF SLURM SCRIPTS
            slurm_list=[f for f in os.listdir(slurm_scripts) if f.endswith(".slurm")]
            print("Total scripts %s"%len(slurm_list))  ## output directory same as out_dir
            
            for i in slurm_list:
                s1= "/".join([slurm_scripts,i]) 
                sub.call("sbatch %s"%s1, shell=True)
            
            slurm_out_list=[f for f in os.listdir(slurm_scripts) if f.endswith(".out")]
            print("Total_out_scripts = %s"%len(slurm_out_list))
            
            # CHUNK 03:CODES AFTER COMPLETION OF JOBS
            fp.GetErrors(slurm_scripts)
            slurm_scripts_samsort = slurm_scripts
            #slurm_out_list_samsort = slurm_out_list()
            return [slurm_scripts_samsort, OUT_DIR_SORTING, INIT_DIR_SORTING]


''' CHUNK 03 ADDREPLACE GRPS '''
def  addreplacegrp(PROJECT, OUT_DIR_SORTING):
            print("> STEP 3: ADDING READ GROUPS")
            PROC4="Addorreplacereadgrps"
            #init_time=time.time()
            INIT_DIR_ADDREPLACEGRPS=  OUT_DIR_SORTING
            OUT_DIR_ADDREPLACEGRPS = "/".join([PROJECT,PROC4])
            if not os.path.exists(OUT_DIR_ADDREPLACEGRPS):
                os.mkdir(OUT_DIR_ADDREPLACEGRPS)
            files_in_directory= os.listdir(INIT_DIR_ADDREPLACEGRPS)
            slurm_scripts = "/".join([OUT_DIR_ADDREPLACEGRPS, "addreplace_scripts"])
            if not os.path.exists(slurm_scripts):
                os.mkdir(slurm_scripts)
            print("> Number of BAM files after SORTING:")
            print([".bam files",len([f for f in files_in_directory if(f.endswith(".bam"))]),
                   ".bai files",len([f for f in files_in_directory if(f.endswith(".bai"))])])
            print("> Proceeding with ADD_REPLACE_REAG_GROUP STEP") ## need to add a bit of time
            time.sleep(1)
            print("> ENTER PICARD TOOL SHOWN HERE: AddOrReplaceReadGroups ")
            tool=input()
            print("java memory required")
            jav = input()
            print("> SLURM OPTIONS DEFAULT y/n:")
            slurm_options= input()
            if slurm_options!="y":
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
                print("no of cores (>=1 or <=8)")
                cpus=input()
                
            time.sleep(1)
            # CHUNK 01 : MAKING SLURM SCRIPTS
            print("Enter seperator in filename")
            sep = input()
            time.sleep(1)
            
            for f in files_in_directory:
               if (f.endswith(".bam")):
                 i=f.split(".")[0]
                 if "N" in i.split(sep):
                     RGSM="NORMAL"
                 else:
                     RGSM="TUMOR"
                 print("> Writing Slurm for %s"%i)
                 time.sleep(1)
                # start_time=time.time()
                 slurm=open("%s/%s.slurm"%(slurm_scripts,i),"w")
                 if slurm_options=="y":
                    options=["#!/bin/bash","\n","#SBATCH --mail-user=%s"%MAIL,"\n",
                         "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%i,"\n","#SBATCH --time=15:00:00",
                         "\n","#SBATCH --nodes=1","\n","#SBATCH --cpus-per-task=1","\n","#SBATCH --output=%s.out"%i,
                         "\n","#SBATCH --mem=8gb","\n","#SBATCH -A %s"%ACCOUNT,"\n","#SBATCH --workdir=%s"%slurm_scripts,"\n"]
                 else:
                     options=["#!/bin/bash","\n","#SBATCH --mail-user=%s"%MAIL,"\n",
                         "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%i,"\n","#SBATCH --time=%s"%t,
                         "\n","#SBATCH --nodes=%s"%nodes,"\n","#SBATCH --cpus-per-task=%s"%cpus,"\n","#SBATCH --output=%s.out"%i,
                         "\n","#SBATCH --mem=%sgb"%mem,"\n","#SBATCH -A %s"%ACCOUNT,"\n","#SBATCH --workdir=%s"%slurm_scripts,"\n"]
                 out_bam = i+".bam" ### need the extension option to be inputted.
                 bam_file="/".join([OUT_DIR_ADDREPLACEGRPS,out_bam]) # out_dir for baserecalibrator and -BQSR for Print Reads
                 for j in range(len(options)):
                     #print(options[j])
                     slurm.write(options[j])
                 slurm.write("\n")
                 slurm.write("%s -Xmx%sg -jar %s %s INPUT=%s OUTPUT=%s RGID=id RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=%s CREATE_INDEX=TRUE VALIDATION_STRINGENCY=LENIENT"%(JAVA,jav,PICARD,tool,"/".join([INIT_DIR_ADDREPLACEGRPS,f]),bam_file,RGSM))
                 slurm.write("\n")
                 slurm.write("wait")
                 time.sleep(1)
               else:
                  continue
              
            slurm.close()
            print("> EXECUTION OF ADDING READ GROUP SCRIPTS")
            # CHUNK 02: EXECUTION OF SLURM SCRIPTS
            slurm_list=[f for f in os.listdir(slurm_scripts) if f.endswith(".slurm")]
            print("Total scripts %s"%len(slurm_list))  ## output directory same as out_dir
            
            for i in slurm_list:
                s1= "/".join([slurm_scripts,i]) 
                sub.call("sbatch %s"%s1, shell=True)
            
            slurm_out_list=[f for f in os.listdir(slurm_scripts) if f.endswith(".out")]
            print("Total_out_scripts = %s"%len(slurm_out_list))
            
            # CHUNK 03:CODES AFTER COMPLETION OF JOBS
            fp.GetErrors(slurm_scripts)
            slurm_script_addrep = slurm_scripts
            #slurm_out_list_addrep = slurm_out_list
            return [slurm_script_addrep , OUT_DIR_ADDREPLACEGRPS,INIT_DIR_ADDREPLACEGRPS]


''' CHUNK 04 MARK DUPLICATE '''
def markduplicates(PROJECT, OUT_DIR_ADDREPLACEGRPS):
        print(">STEP4: REMOVING DUPLICATES BY: MarkDuplicates")
        PROC3="Duplicates"
        #init_time=time.time()
        INIT_DIR_DUPLICATE =  OUT_DIR_ADDREPLACEGRPS 
        OUT_DIR_DUPLICATE = "/".join([PROJECT,PROC3])
        os.mkdir(OUT_DIR_DUPLICATE)
        files_in_directory= os.listdir(INIT_DIR_DUPLICATE)
        slurm_scripts = "/".join([OUT_DIR_DUPLICATE, "duprem_scripts"])
        os.mkdir(slurm_scripts)
        print("> Number of BAM files after ADDING READ GROUPS:")
        print([".bam files",len([f for f in files_in_directory if(f.endswith(".bam"))]),
               ".bai files",len([f for f in files_in_directory if(f.endswith(".bai"))])])
        time.sleep(2)
        print("> Proceeding with DUPLICATE REMOVAL STEP") ## need to add a bit of time
        time.sleep(1)
        print("> ENTER PICARD TOOL SHOWN HERE: MarkDuplicates")
        tool=input()
        print("java memory required")
        jav = input()
        print("> SLURM OPTIONS DEFAULT y/n:")
        slurm_options= input()
        if slurm_options!="y":
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
            print("no of cores (>=1 or <=8)")
            cpus=input()
            
        time.sleep(1)
        # CHUNK 01 : MAKING SLURM SCRIPTS
        print("> Building SLURM SCRIPTS")
        time.sleep(1)
        
        for f in files_in_directory:
           if (f.endswith(".bam")):
             i=re.split("\.",f)[0]
             print("> Writing Slurm for %s"%i)
             time.sleep(1)
             #start_time=time.time()
             slurm=open("%s/%s.slurm"%(slurm_scripts,i),"w")
             if slurm_options=="y":
                options=["#!/bin/bash","\n","#SBATCH --mail-user=%s"%MAIL,"\n",
                     "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%i,"\n","#SBATCH --time=15:00:00",
                     "\n","#SBATCH --nodes=1","\n","#SBATCH --cpus-per-task=1","\n","#SBATCH --output=%s.out"%i,
                     "\n","#SBATCH --mem=8gb","\n","#SBATCH -A %s"%ACCOUNT,"\n","#SBATCH --workdir=%s"%slurm_scripts,"\n"] # IF PICARD TO BE LOADED OR NOT
             else:
                 options=["#!/bin/bash","\n","#SBATCH --mail-user=%s"%MAIL,"\n",
                     "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%i,"\n","#SBATCH --time=%s"%t,
                     "\n","#SBATCH --nodes=%s"%nodes,"\n","#SBATCH --cpus-per-task=%s"%cpus,"\n","#SBATCH --output=%s.out"%i,
                     "\n","#SBATCH --mem=%sgb"%mem,"\n","#SBATCH -A %s"%ACCOUNT,"\n","#SBATCH --workdir=%s"%slurm_scripts,"\n"]
             out_bam = i+".bam" ### need the extension option to be inputted.
             bam_file="/".join([OUT_DIR_DUPLICATE,out_bam]) # out_dir for baserecalibrator and -BQSR for Print Reads
             metric  = i+".metrics.txt"
             met_file= "/".join([OUT_DIR_DUPLICATE, metric])
             for j in range(len(options)):
                 #print(options[j])
                 slurm.write(options[j])
             slurm.write("\n")
             slurm.write("%s -Xmx%sg -jar %s %s REMOVE_DUPLICATES=True CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT INPUT=%s OUTPUT=%s %s"%(JAVA,jav,PICARD,tool,"/".join([INIT_DIR_DUPLICATE,f]),bam_file, met_file))
             slurm.write("\n")
             slurm.write("wait")
             time.sleep(1)
           else:
              continue
          
        slurm.close()
        print("> EXCETION OF DUPILCATE REMOVAL SCRIPTS")
        # CHUNK 02: EXECUTION OF SLURM SCRIPTS
        slurm_list=[f for f in os.listdir(slurm_scripts) if f.endswith(".slurm")]
        print("Total scripts %s"%len(slurm_list))  ## output directory same as out_dir
        
        for i in slurm_list:
            s1= "/".join([slurm_scripts,i]) 
            sub.call("sbatch %s"%s1, shell=True)
        
        slurm_out_list=[f for f in os.listdir(slurm_scripts) if f.endswith(".out")]
        print("Total_out_scripts = %s"%len(slurm_out_list))
        
        # CHUNK 03:CODES AFTER COMPLETION OF JOBS
        fp.GetErrors(slurm_scripts)
        slurm_scripts_markduplicates = slurm_scripts
        #slurm_out_list_markduplicates =  slurm_out_list
        return [slurm_scripts_markduplicates, OUT_DIR_DUPLICATE,INIT_DIR_DUPLICATE]


''' CHUNK 05 BASE RECALIBRATION '''
def baserecalibration(PROJECT, OUT_DIR_DUPLICATE):
            print("> step 5: PROCEEDING WITH BASE RECALIBRATION ")
            time.sleep(2)
            ''' CHUNK 05A BASE RECALIBRATE - CREATING RECAL.GRP FILE'''
            PROC5="Baserecalibration"
            #init_time=time.time()
            INIT_DIR_BASERECALIBRATION=  OUT_DIR_DUPLICATE
            OUT_DIR_BASERECALIBRATION = "/".join([PROJECT,PROC5])
            if not os.path.exists(OUT_DIR_BASERECALIBRATION):
                os.mkdir(OUT_DIR_BASERECALIBRATION)
            files_in_directory= os.listdir(INIT_DIR_BASERECALIBRATION)
            slurm_scripts_grp = "/".join([OUT_DIR_BASERECALIBRATION, "baserecalib_recalgrp_scripts"])
            slurm_scripts_bam = "/".join([OUT_DIR_BASERECALIBRATION, "printreads_scripts"])
            if not os.path.exists(slurm_scripts_grp):
                os.mkdir(slurm_scripts_grp)
            if not os.path.exists(slurm_scripts_bam):
                os.mkdir(slurm_scripts_bam)
            print("> Number of BAM files after DUPICATES REMOVED:")
            print([".bam files",len([f for f in files_in_directory if(f.endswith(".bam"))]),
                   ".bai files",len([f for f in files_in_directory if(f.endswith(".bai"))])])
            time.sleep(2)
            print("java memory required")
            jav = input()
            print("> SLURM OPTIONS DEFAULT y/n:")
            slurm_options= input()
            if slurm_options!="y":
                print(">SLURM OPTIONS REQUIRED(int)")
                time.sleep(1)
                print("Hours")
                hr=input()
                print("Minutes")
                mins=input()
                mints=str(mins).zfill(2)
                print("Seconds")
                sec=input()
                secs=str(sec).zfill(2)
                t=":".join([str(hr), str(mints), str(secs)])
                print("Slurm Memory required >=16")
                mem=input()
                print("Nodes required")
                nodes=input()
            print("no of cpus per thread required for parallelization (-nct option used)")
            print ("2<=(-nct)<8-Will distribute the job among cpus and will decrease the time taken")
            cpus=input()
            time.sleep(1)
            # CHUNK 01 : MAKING SLURM SCRIPTS
            print("> Building SLURM SCRIPTS FOR RECAL_GRP FILES AND FOLLOWING BAM FILES")
            time.sleep(1)
            
            for f in files_in_directory:
               if (f.endswith(".bam")):
                 i=re.split("\.",f)[0]
                 print("> Writing Slurm for %s"%i)
                 time.sleep(1)
                 #start_time=time.time()
                 slurm=open("%s/%s.slurm"%(slurm_scripts_grp,i),"w")
                 slurm1=open("%s/%s.slurm"%(slurm_scripts_bam,i),"w")
                 if slurm_options=="y":
                    options=["#!/bin/bash","\n","#SBATCH --mail-user=%s"%MAIL,"\n",
                         "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%i,"\n","#SBATCH --time=30:00:00",
                         "\n","#SBATCH --nodes=1","\n","#SBATCH --cpus-per-task=%s"%cpus,"\n","#SBATCH --output=%s.out"%i,
                         "\n","#SBATCH --mem=16gb","\n","#SBATCH -A %s"%ACCOUNT,"\n", "%s"%g_mod,"\n"]
                 else:
                     options=["#!/bin/bash","\n","#SBATCH --mail-user=%s"%MAIL,"\n",
                         "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%i,"\n","#SBATCH --time=%s"%t,
                         "\n","#SBATCH --nodes=%s"%nodes,"\n","#SBATCH --cpus-per-task=%s"%cpus,"\n","#SBATCH --output=%s.out"%i,
                         "\n","#SBATCH --mem=%sgb"%mem,"\n","#SBATCH -A %s"%ACCOUNT,"\n","%s"%g_mod,"\n"]
                 calib_grp= i+".recal.grp" ### need the extension option to be inputted.
                 grp_file="/".join([OUT_DIR_BASERECALIBRATION,calib_grp]) # out_dir for baserecalibrator and -BQSR for Print Reads
                 calib_bam=i+".bam"
                 bam_file="/".join([OUT_DIR_BASERECALIBRATION,calib_bam])
                 for j in range(len(options)):
                     #print(options[j])
                     slurm.write(options[j])
                     slurm1.write(options[j])
                 slurm.write("#SBATCH --workdir=%s"%slurm_scripts_grp)
                 slurm.write("\n")
                 slurm.write("\n")
                 slurm1.write("#SBATCH --workdir=%s"%slurm_scripts_bam)
                 slurm1.write("\n")
                 slurm1.write("\n")
                 slurm.write("%s -Xmx%sg -jar %s -T BaseRecalibrator -nct %s -R %s -I %s -knownSites %s -o %s "%(JAVA,jav,GATK,cpus,REF_GENOME,"/".join([INIT_DIR_BASERECALIBRATION,f]),vcf,grp_file))
                 slurm1.write("%s -Xmx%sg -jar %s -T PrintReads -nct %s -R %s -I %s -BQSR %s -o %s "%(JAVA,jav,GATK,cpus,REF_GENOME,"/".join([INIT_DIR_BASERECALIBRATION,f]),grp_file,bam_file))
                 slurm.write("\n")
                 slurm1.write("\n")
                 slurm.write("wait")
                 slurm1.write("wait")
                 time.sleep(1)
               else:
                  continue
              
            slurm.close()
            print("> EXECUTION OF BASE_RECALIBRATOR SCRIPTS")
            # CHUNK 02: EXECUTION OF SLURM SCRIPTS
            slurm_list_grp=[f for f in os.listdir(slurm_scripts_grp) if f.endswith(".slurm")]
            print("Total scripts %s"%len(slurm_list_grp))  ## output directory same as out_dir
            time.sleep(1)
            for i in slurm_list_grp:
                s1= "/".join([slurm_scripts_grp,i]) 
                sub.call("sbatch %s"%s1, shell=True)
            time.sleep(2)
            slurm_out_list_grp=[f for f in os.listdir(slurm_scripts_grp) if f.endswith(".out")]
            print("Total_out_scripts = %s"%len(slurm_out_list_grp))
            
            ## CHUNK 02A
            fp.GetErrors(slurm_scripts_grp)
            #time.sleep(5)
            #success_out=fp.Success(slurm_scripts_grp,slurm_out_list_grp)
            #fp.forOutputLog(OUT_DIR_BASERECALIBRATION,slurm_scripts_grp,INIT_DIR_BASERECALIBRATION,success_out, "recalgrp", ".recal.grp")
            #time.sleep(2)
            print("> DONE WITH RECAL.GRP FILES")
            time.sleep()
            print("> EXECUTION OF PRINT READS SCRIPTS")
            # CHUNK 03: EXECUTION OF SLURM SCRIPTS
            slurm_list_bam=[f for f in os.listdir(slurm_scripts_bam) if f.endswith(".slurm")]
            print("Total scripts %s"%len(slurm_list_bam))  ## output directory same as out_dir
            time.sleep(1)
            for i in slurm_list_bam:
                s1= "/".join([slurm_scripts_bam,i]) 
                sub.call("sbatch %s"%s1, shell=True)
            
            time.sleep(2)
            slurm_out_list_bam=[f for f in os.listdir(slurm_scripts_bam) if f.endswith(".out")]
            print("Total_out_scripts = %s"%len(slurm_out_list_bam))
            
            # CHUNK 03A:CODES AFTER COMPLETION OF JOBS
            fp.GetErrors(slurm_scripts_bam)
            slurm_scripts_baserecalib =  slurm_scripts_bam
            #slurm_out_list_baserecalib = slurm_out_list_bam
            return [slurm_scripts_baserecalib, OUT_DIR_BASERECALIBRATION, INIT_DIR_BASERECALIBRATION]
        


''' CHUNK 06 INDEL REALIGNMENT '''
def indelrealignment(PROJECT, OUT_DIR_BASERECALIBRATION):
            print("> step 6: PROCEEDING WITH INDEL REALIGNMENT ")
            time.sleep(2)
            PROC6="Indelrealigner"
            #init_time=time.time()
            INIT_DIR_INDEL=  OUT_DIR_BASERECALIBRATION
            OUT_DIR_INDEL = "/".join([PROJECT,PROC6])
            if not os.path.exists(OUT_DIR_INDEL):
                os.mkdir(OUT_DIR_INDEL)
            files_in_directory= os.listdir(INIT_DIR_INDEL)
            slurm_scripts_interval = "/".join([OUT_DIR_INDEL, "targetinterval_scripts"])
            slurm_scripts_bam = "/".join([OUT_DIR_INDEL, "indelrealign_scripts"])
            if not os.path.exists(slurm_scripts_interval):
                os.mkdir(slurm_scripts_interval)
            if not os.path.exists(slurm_scripts_bam):
                os.mkdir(slurm_scripts_bam)
            print("> Number of BAM files after BASE RECALIBRATION")
            print([".bam files",len([f for f in files_in_directory if(f.endswith(".bam"))]),
                   ".bai files",len([f for f in files_in_directory if(f.endswith(".bai"))])])
            time.sleep(2)
            print("java memory required")
            jav = input()
            print("> SLURM OPTIONS DEFAULT y/n:")
            slurm_options= input()
            if slurm_options!="y":
                print(">SLURM OPTIONS REQUIRED(int)")
                time.sleep(1)
                print("Hours")
                hr=input()
                print("Minutes")
                mins=input()
                mints=str(mins).zfill(2)
                print("Seconds")
                sec=input()
                secs=str(sec).zfill(2)
                t=":".join([str(hr), str(mints), str(secs)])
                print("Slurm Memory required >=16")
                mem=input()
                print("Nodes required")
                nodes=input()
                print("no of cores (>=1 or <=8)")
                cpus=input()
                
            time.sleep(1)
            # CHUNK 01 : MAKING SLURM SCRIPTS
            print("> Building SLURM SCRIPTS FOR TARGET INTERVAL FILES & FOLLOWING BAM FILES")
            time.sleep(1)
            
            for f in files_in_directory:
               if (f.endswith(".bam")):
                 i=re.split("\.",f)[0]
                 print("> Writing Slurm for %s"%i)
                 time.sleep(1)
                 #start_time=time.time()
                 slurm=open("%s/%s.slurm"%(slurm_scripts_interval,i),"w")
                 slurm1=open("%s/%s.slurm"%(slurm_scripts_bam,i),"w")
                 if slurm_options=="y":
                    options=["#!/bin/bash","\n","#SBATCH --mail-user=%s"%MAIL,"\n",
                         "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%i,"\n","#SBATCH --time=15:00:00",
                         "\n","#SBATCH --nodes=1","\n","#SBATCH --cpus-per-task=1","\n","#SBATCH --output=%s.out"%i,
                         "\n","#SBATCH --mem=16gb","\n","#SBATCH -A %s"%ACCOUNT,"\n","%s"%g_mod,"\n"]
                 else:
                     options=["#!/bin/bash","\n","#SBATCH --mail-user=%s"%MAIL,"\n",
                         "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%i,"\n","#SBATCH --time=%s"%t,
                         "\n","#SBATCH --nodes=%s"%nodes,"\n","#SBATCH --cpus-per-task=%s"%cpus,"\n","#SBATCH --output=%s.out"%i,
                         "\n","#SBATCH --mem=%sgb"%mem,"\n","#SBATCH -A %s"%ACCOUNT,"\n","%s"%g_mod,"\n"]
                 target_interval= i+".realigner.intervals" ### need the extension option to be inputted.
                 target_file="/".join([OUT_DIR_INDEL,target_interval]) # out_dir for baserecalibrator and -BQSR for Print Reads
                 indel_bam=i+".bam"
                 bam_file="/".join([OUT_DIR_INDEL,indel_bam])
                 for j in range(len(options)):
                     #print(options[j])
                     slurm.write(options[j])
                     slurm1.write(options[j])
                 slurm.write("#SBATCH --workdir=%s"%slurm_scripts_interval)
                 slurm.write("\n")
                 slurm1.write("\n")
                 slurm1.write("#SBATCH --workdir=%s"%slurm_scripts_bam)
                 slurm1.write("\n")
                 slurm1.write("\n")
                 slurm.write("%s -Xmx%sg -jar %s -T RealignerTargetCreator -R %s -I %s -o %s "%(JAVA,jav,GATK,REF_GENOME,"/".join([INIT_DIR_INDEL,f]),target_file))
                 slurm1.write("%s -Xmx%sg -jar %s -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s "%(JAVA,jav,GATK,REF_GENOME,"/".join([INIT_DIR_INDEL,f]),target_file,bam_file))
                 slurm.write("\n")
                 slurm1.write("\n")
                 slurm.write("wait")
                 slurm1.write("wait")
                 time.sleep(1)
               else:
                  continue
              
            slurm.close()
            print("> EXECUTION OF TARGET INTERVAL SCRIPTS")
            # CHUNK 02: EXECUTION OF SLURM SCRIPTS
            slurm_list_interval=[f for f in os.listdir(slurm_scripts_interval) if f.endswith(".slurm")]
            print("Total scripts %s"%len(slurm_list_interval))  ## output directory same as out_dir
            time.sleep(1)
            for i in slurm_list_interval:
                s1= "/".join([slurm_scripts_interval,i]) 
                sub.call("sbatch %s"%s1, shell=True)
            time.sleep(2)
            slurm_out_list_interval=[f for f in os.listdir(slurm_scripts_interval) if f.endswith(".out")]
            print("Total_out_scripts = %s"%len(slurm_out_list_interval))
            
            ## CHUNK 02A
            fp.GetErrors(slurm_scripts_interval)
            print("> DONE WITH TARGET CREATOR WITHOUT ERRORS")
            time.sleep(1)
            print("> EXECUTION OF INDEL REALIGN SCRIPTS")
            # CHUNK 03: EXECUTION OF SLURM SCRIPTS
            slurm_list_bam=[f for f in os.listdir(slurm_scripts_bam) if f.endswith(".slurm")]
            print("Total scripts %s"%len(slurm_list_bam))  ## output directory same as out_dir
            time.sleep(1)
            for i in slurm_list_bam:
                s1= "/".join([slurm_scripts_bam,i]) 
                sub.call("sbatch %s"%s1, shell=True)
            
            time.sleep(2)
            slurm_out_list_bam=[f for f in os.listdir(slurm_scripts_bam) if f.endswith(".out")]
            print("Total_out_scripts = %s"%len(slurm_out_list_bam))
            
            # CHUNK 03A:CODES AFTER COMPLETION OF JOBS
            fp.GetErrors(slurm_scripts_bam)
            slurm_scripts_indel = slurm_scripts_bam
            #slurm_out_list_indel = slurm_out_list_bam
            return [slurm_scripts_indel,OUT_DIR_INDEL, INIT_DIR_INDEL]

def hsmetrics(PROJECT,INIT):
            print(">STEP: RUNNING HSMETRICS")
            PROC="HSMETRICS"
            #init_time=time.time()
            #INIT_DIR_METRICS =  OUT_DIR_ALIGN ## output directory of aligned files 
            OUT_DIR_HSMETRICS = "/".join([PROJECT,PROC]) ## ooutput directory for sorted files will be the directory made using project location and name of the process
            if not os.path.exists(OUT_DIR_HSMETRICS):
                os.mkdir(OUT_DIR_HSMETRICS)
            files_in_directory= os.listdir(INIT) ## getting the list of aligned files 
            slurm_scripts = "/".join([OUT_DIR_HSMETRICS, "HSMETRICS_scripts"])
            if not os.path.exists(slurm_scripts):
                os.mkdir(slurm_scripts)
            print("> Number of BAM FILES:")
            print([".bam files",len([f for f in files_in_directory if(f.endswith(".bam"))])])
            time.sleep(1)
            print("> Proceeding with forming  metrics") ## need to add a bit of time
            time.sleep(1)
            print("> ENTER PICARD TOOL SHOWN:CalculateHsMetrics ") ## need to see whether to add it or not cz unlike gatk toolsonly one is required
            tool=input()
            print("> ENTER LOCATION OF BAIT's FILE")
            bait = input()
            print(">ENTER LOCATION OF TARGET FILE")
            target = input()
            print("> SLURM OPTIONS DEFAULT y/n:")
            slurm_options= input()
            if slurm_options!="y":
                print("hours")
                hr=input()
                print("minutes")
                mins=input()
                mints=str(mins).zfill(2)
                print("seconds")
                sec=input()
                secs=str(sec).zfill(2)
                t=":".join([str(hr), str(mints), str(secs)])
                print("slurm memory required")
                mem=input()
                print("nodes required")
                nodes=input()
                print("no of cores (>=1 or <=8)")
                cpus=input()
                
            time.sleep(1)
            # CHUNK 01 : MAKING SLURM SCRIPTS
            print("> Building SLURM SCRIPTS")
            time.sleep(1)
            
            for f in files_in_directory:
               if (f.endswith(".bam")):
                 i=re.split("\.",f)[0]
                 print("> Writing Slurm for %s"%i)
                 time.sleep(1)
                 #start_time=time.time()
                 slurm=open("%s/%s.slurm"%(slurm_scripts,i),"w")
                 if slurm_options=="y":
                    options=["#!/bin/bash","\n","#SBATCH --mail-user=%s"%MAIL,"\n",
                         "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%i,"\n","#SBATCH --time=80:00:00",
                         "\n","#SBATCH --nodes=1","\n","#SBATCH --cpus-per-task=1","\n","#SBATCH --output=%s.out"%i,
                         "\n","#SBATCH --mem=8gb","\n","#SBATCH -A %s"%ACCOUNT,"\n", "#SBATCH --workdir=%s"%slurm_scripts, "\n"]
                 else:
                     options=["#!/bin/bash","\n","#SBATCH --mail-user=%s"%MAIL,"\n",
                         "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%i,"\n","#SBATCH --time=%s"%t,
                         "\n","#SBATCH --nodes=%s"%nodes,"\n","#SBATCH --cpus-per-task=%s"%cpus,"\n","#SBATCH --output=%s.out"%i,
                         "\n","#SBATCH --mem=%sgb"%mem,"\n","#SBATCH -A %s"%ACCOUNT,"\n","#SBATCH --workdir=%s"%slurm_scripts, "\n"]
                 out_txt = i+".txt" ### need the extension option to be inputted.
                 txt_file="/".join([OUT_DIR_HSMETRICS,out_txt]) # out_dir for baserecalibrator and -BQSR for Print Reads
                 for j in range(len(options)):
                     #print(options[j])
                     slurm.write(options[j])
                 slurm.write("\n")
                 slurm.write("%s -Xmx8g -jar %s %s I=%s O=%s  R=%s BAIT_INTERVALS=%s TARGET_INTERVALS=%s"%(JAVA, PICARD,tool,"/".join([INIT,f]),txt_file, REF_GENOME, bait,target))
                 slurm.write("\n")
                 slurm.write("wait")
                 time.sleep(1)
               else:
                  continue
            
            slurm.close()
            print(">EXECUTION OF SCRIPTS")
            time.sleep(1)
            # CHUNK 02: EXECUTION OF SLURM SCRIPTS
            slurm_list=[f for f in os.listdir(slurm_scripts) if f.endswith(".slurm")]
            print("Total scripts %s"%len(slurm_list))  ## output directory same as out_dir
            for i in slurm_list:
                s1= "/".join([slurm_scripts,i]) 
                sub.call("sbatch %s"%s1, shell=True)
            time.sleep(5)
            slurm_out_list=[f for f in os.listdir(slurm_scripts) if f.endswith(".out")]
            print("Total_out_scripts = %s"%len(slurm_out_list))
            
            # CHUNK 03:CODES AFTER COMPLETION OF JOBS
            fp.GetErrors(slurm_scripts)
            slurm_scripts_hsmetrics = slurm_scripts
            #slurm_out_list_samsort = slurm_out_list()
            return [slurm_scripts_hsmetrics, OUT_DIR_HSMETRICS, INIT]        

##############################################################################################################################
#print("\t*********** SUCCESSFULLY COMPLETED DATA PROCESSING PIPELINE ****************\t")
