#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 19:31:22 2017

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
print("location of python files")
PYTHON_PATH = input()
sys.path.append(os.path.abspath(PYTHON_PATH))
import FUNCTIONS_FOR_PIPELINE as fp




print("> Proceeding with MUTATION CALLING after DATA PROCESSING")
time.sleep(1)
'''
The mutation caller used are varscan and mutect lod-10 this pipeline has three functions:
1. Variation caller using Mutect version 1.1.7
2. Mpileup function using samtools for generating mpileup file as a varscan Input
3. Somatic tool for varscan variant calling
def AskLocations():
    
    Name: Ask Locations
    Description > This function is required for the pipeline (must be called irrespective of any function called)
                > Asks the locations of PROJECT,REF_GENOME,vcf, JAVA, b_mod, BWA, g_mod, GATK, PICARD, FASTQ 
    
    global PROJECT, REF_GENOME, vcf, JAVA, PICARD, g_mod, GATK,sam_mod, samtools, target_intervals,INIT, MAIL, ACCOUNT
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
    print(">Load Slurm SAMTOOLS module: (y/n)")
    b_in=input()
    if b_in != "y":
        print(">ENTER LOCATION TO SAMTOOLS (version-1.2 or higher)")
        sam_mod= "\t"
        samtools=input()
    else:
        print("Slurm Module for SAMTOOLS (version-1.5) Is Loaded")
        sam_mod = "module load samtools"
        samtools= "samtools"
    print(">ENTER LOCATION OF PER CHROMOSOME TARGET INTERVAL FILE")
    target_intervals = input()
    print(">ENTER LOCATION OF INIT FILES")
    INIT = input()
    print(">ENTER USER MAIL FOR SLURM SCRIPT")
    MAIL = input()
    print(">ENTER ACCOUNT WITH CPU CORES ASSIGNED FOR SLURM JOBS")
    ACCOUNT = input()
    return PROJECT, REF_GENOME, vcf, JAVA,PICARD, g_mod, GATK, sam_mod, samtools, target_intervals, INIT, MAIL, ACCOUNT

'''
global PROJECT, REF_GENOME, vcf, JAVA, g_mod, GATK, sam_mod, samtools, target_intervals,INIT, MAIL, ACCOUNT 
PROJECT = "/mnt/pan/Data16/jxk906/PROJECTS/TNBC"
REF_GENOME="/mnt/pan/Data16/jxk906/REFERENCE/hg19.fa"
vcf="/mnt/pan/Data16/jxk906/split_vcf/Merged_dbsnp.vcf"
JAVA="/mnt/pan/Data16/Tools/jdk1.8.0_91/bin/java"
g_mod="module load gatk"
GATK = "$GATK"
sam_mod = "module load samtools"
samtools = "samtools"
target_intervals = "/mnt/pan/Data16/jxk906/target_file"
INIT = "/mnt/pan/Data16/jxk906/PROJECTS/TNBC/INDELREALIGNER"
MAIL = "jxk906@case.edu"
ACCOUNT = "lxh259"
## for getting the T/N pairs without having files :
def get_pairs():
        '''
        Name : get_pairs
        Description: Returns the dictionary of tumor and normal pairs
        '''
        global INIT_FILES_PATH, sep
        print("> Locate Tumor and Normal pairs Enter Location of .bam Files ")
        INIT_FILES_PATH=input()
        print("> Require seperator of file name")
        sep= input()
        files_in_directory = os.listdir(INIT_FILES_PATH)
        files_for_match = files_in_directory[1:]
        pairs = {}
        for f in files_in_directory:
            for fm in files_for_match:
                if f.endswith(".bam") and fm.endswith(".bam"):
                    p1 = f.split(".")[0].split(sep) ## input seperator in filename (# here in my case it is _)
                    p2 = fm.split(".")[0].split(sep) ## input the seperator in filename
                    if "T" in p1 and "N" in p2:
                        i= p1.index("T"); j = p2.index("N")
                        if p1[:i] == p2[:j]:
                            print("Matched Pairs %s and %s"%(f,fm))
                            pairs[f]=fm 
                    else:
                        continue
        return pairs, sep
        print("NUMBER OF PAIRS %s"%len(pairs))

## NOTE -- the dictionary has KEYS = TUMOR VALUES=NORMALS
''' The mutect or mpileup require bam files in pairs which is why getpairs is the required funtion to run prior to mutectl10 or mpileup,
PAIRS= getpairs(), PAIRS is a python dictionary with KEYS=TUMOR samples and VALUES=NORMAL samples '''
### READ THE TUMOR AND NORMAL FILE with pair

def get_pairs_from_file():
        print("Location of Tab delimited file with Tumor and matched normal pairs")
        FILE_PATH = input()
        print("input file separator")
        sep1 = input()
        df_pairs = pd.read_table(FILE_PATH, sep = sep1, index_col=False, dtype=None)
        pairs = {}
        cols = df_pairs.columns
        vals = df_pairs.values
        if cols[0]=="Tumor" and cols[1]== "Normal":
            for i in vals:
                pairs[i[0]]=i[1]
        elif cols[0]=="Normal" and cols[1]=="Tumor":
            for i in vals:
                pairs[i[1]]=i[0]
        else:
            "None"
        return pairs
        print("NUMBER OF PAIRS %s"%len(pairs))
        
def mutect_l_10(PROJECT,INIT, pairs):
    PROC = "mutectl10" #
    INIT_DIR_MUT= INIT_FILES_PATH# initial files need to be indel realigner
    OUT_DIR_MUT="/".join([PROJECT,PROC])#
    if not os.path.exists(OUT_DIR_MUT):
        os.mkdir(OUT_DIR_MUT)#
    slurm_scripts = "/".join([OUT_DIR_MUT, "mutectl10_scripts"])#
    if not os.path.exists(slurm_scripts):
        os.mkdir(slurm_scripts)#
    print("> Number of BAM files :") ## has to be indel realigned bam files
    print([".bam files",len([f for f in os.listdir(INIT_DIR_MUT) if(f.endswith(".bam"))])])
    time.sleep(1)
    print(">Proceeding with MUTECT:")
    time.sleep(1)
    print(">ENTER full jar file location of MUTECT:") 
    MUT = input()
    print(">ENTER location of java (version - 1.7.0_80):")
    jav = input()
    print("> ENTER TOOL shown here: MuTect") ## need to see whether to add it or not cz unlike gatk toolsonly one is required
    tool=input()
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
    for k,v in pairs.items():
        print("> Writing Slurm for %s"%k)
        time.sleep(1)
        #start_time=time.time()
        slurm=open("%s/%s.slurm"%(slurm_scripts,k.split(sep)[0]),"w")
        if slurm_options=="y":
            options=["#!/bin/bash","\n","#SBATCH --mail-user=%s"%MAIL,"\n",
             "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%k.split(sep)[0],"\n","#SBATCH --time=80:00:00",
             "\n","#SBATCH --nodes=1","\n","#SBATCH --cpus-per-task=%s"%cpus,"\n","#SBATCH --output=%s.out"%k.split(sep)[0],
             "\n","#SBATCH --mem=16gb","\n","#SBATCH -A %s"%ACCOUNT,"\n","#SBATCH --workdir=%s"%slurm_scripts,"\n"]
        else:
            options=["#!/bin/bash","\n","#SBATCH --mail-user=%s"%MAIL,"\n",
             "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%k.split(sep)[0],"\n","#SBATCH --time=%s"%t,
             "\n","#SBATCH --nodes=%s"%nodes,"\n","#SBATCH --cpus-per-task=%s"%cpus,"\n","#SBATCH --output=%s.out"%k.split(sep)[0],
             "\n","#SBATCH --mem=%sgb"%mem,"\n","#SBATCH -A %s"%ACCOUNT,"\n","#SBATCH --workdir=%s"%slurm_scripts,"\n"]
        out_text = k.split("_")[0]+".stats.txt" ### need the extension option to be inputted.
        text_file="/".join([OUT_DIR_MUT,out_text]) # out_dir for baserecalibrator and -BQSR for Print Reads
        out_vcf = k.split("_")[0] + ".vcf"
        vcf_file = "/".join([OUT_DIR_MUT, out_vcf])
        for j in range(len(options)):
         #print(options[j])
            slurm.write(options[j])
        slurm.write("\n")
        slurm.write("%s -jar %s -T %s -R %s -I:normal %s -I:tumor %s --dbsnp %s --tumor_lod 10 --fraction_contamination 0.01 -o %s -vcf %s"%(jav,MUT,tool,REF_GENOME,"/".join([INIT_DIR_MUT,v]), "/".join([INIT_DIR_MUT,k]),vcf,text_file,vcf_file))
        slurm.write("\n")
        slurm.write("wait")
        time.sleep(1)

    slurm.close()
    print(">EXECUTION OF MUTECT SCRIPTS")
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
    fp.GetErrors(slurm_scripts) ## NEED TO CHECK GETERRORS AND SUCCESS FUNCTION.
    time.sleep(180)
    #Done, no_msg =fp.Success(slurm_scripts)
    #fp.forOutputLog(OUT_DIR_ALIGN, slurm_scripts,INIT_FASTQ_FILES, success_out, "align", ".sam")
   # THIS HAS TO BE MODIFIED 
    #no_m=[k for k, v in no_msg.items()]
    #don=[k for k, v in Done.items()]
    #for i in no_m:
     #   if not i in don:
      #      don.append(i)
       # else:
        #    continue
    out_log= open("/".join([slurm_scripts,"mutect-L10.log"]), "w") ## making the combined file of last lines
    list_of_init_files = os.listdir(INIT_FILES_PATH)
    init_bam_files=[f for f in list_of_init_files if f.endswith(".bam")] ## initial bam files will be the same
    c=0
    out = [f for f in os.listdir(slurm_scripts) if f.endswith(".out")]
    #total_read_list = []
    #reads_left=[]
    final= {}
    ## Mutect stat_out, vcf_out and slurm_out are in the format TN05.stats.txt, TN05.vcf and TN05.out 
    ## BUT there will be 2 indel realigned files for every sample
    ## Therefore this creating a final dictionary with KEY=out_file(TN05.out) and VALUE = [TN05_T.bam, TN05_N.bam]
    ## easy to write the size and for reassurance of path
    for k in out:
        list_indel_files = []
        for i in init_bam_files:
            if k.split('.')[0]== i.split('.')[0].split(sep)[0]:
                list_indel_files.append(i)
            else:
                final[k]=list_indel_files
    for k,v in final.items():
            #for i in init_bam_files:
            #  if i.split(".")[0]==k.split(".")[0]:
                  c+=1
                  init_path_1= "/".join([INIT_FILES_PATH,v[0]])  # joining path to script CANNOT TAKE INTIAL PATH CZ THE THERE WILL BE TWO FILES
                  init_bam_size_1 =((os.path.getsize(init_path_1))/1000000000)
                  
                  init_path_2= "/".join([INIT_FILES_PATH,v[1]])  # joining path to script CANNOT TAKE INTIAL PATH CZ THE THERE WILL BE TWO FILES
                  init_bam_size_2=((os.path.getsize(init_path_2))/1000000000)
                  
                  recent_output_vcf= "/".join([OUT_DIR_MUT,k.split('.')[0]+".vcf"])
                  out_vcf_size=((os.path.getsize(recent_output_vcf))/1000000000)
                  
                  recent_output_stat= "/".join([OUT_DIR_MUT,k.split('.')[0]+".stats.txt"])
                  out_stat_size=((os.path.getsize(recent_output_stat))/1000000000)
                 
                  out_log.write("%s. %s\n"%(c,k))
                  #out_log.write("\n")
                  out_log.write("INPUT_1: %s - %s gb\n"%(init_path_1,init_bam_size_1))
                  out_log.write("INPUT_2: %s - %s gb\n"%(init_path_2,init_bam_size_2))
                 # out_log.write("\n")
                  out_log.write("OUTPUT_VCF: %s - %s gb\n"%(recent_output_vcf,out_vcf_size))
                  out_log.write("OUTPUT_STAT.TXT: %s - %s gb\n"%(recent_output_stat,out_stat_size))
                  #out_log.write("\n")
                  slurm_out_file= "/".join([slurm_scripts,k])
                  p3= sub.Popen(['grep','-i','Program Args',slurm_out_file],stdout=sub.PIPE) 
                  soutput,sinput=p3.communicate()
                  out_log.write("%s\n"%(" ".join(((soutput.decode("utf-8")).split())[4:])))## writing the command used
                  lines=fp.GetLastNLines(20,slurm_out_file)
                  for i in lines[-8:-5]:
                      out_log.write("%s"%(" ".join([j for j in i])))
                      out_log.write("\n")
                  #out_log.write("%s"%lines)
                  #for j in lines[-6]:
                  #print([c,i,lines[-6][-4:-1]])
                  #reads_left.append([c,k,int(lines[-6][-4])-int(lines[-6][4])])
                  #total_read_list.append([c,k,lines[-6][-4]])
                  #out_log.write("%s,%s,%s"%(c,i,lines[-6][-4:-1][0]))
                     #out_log.write("%s"%i[-4:-1])
                     #out_log.write("\n")
                  out_log.write("\n")
             # else:
              #    continue
    out_log.close()
    #df_total_reads = pd.DataFrame.from_records(total_read_list)
    #df_reads_left = pd.DataFrame.from_records(reads_left)
    #return df_reads_left, df_total_reads
 #   end_time=time.time()
  #  print("> ELAPSED  TIME FOR ALIGNMENT %s hrs"%((init_time-end_time)/3600))
   # time.sleep(6)'''

def mutect2(PROJECT,INIT, pairs):
    PROC = "Mutect2" #
    INIT_DIR_MUT= INIT_FILES_PATH# initial files need to be indel realigner
    OUT_DIR_MUT="/".join([PROJECT,PROC])#
    if not os.path.exists(OUT_DIR_MUT):
        os.mkdir(OUT_DIR_MUT)#
    slurm_scripts = "/".join([OUT_DIR_MUT, "Mutect2_scripts"])#
    if not os.path.exists(slurm_scripts):
        os.mkdir(slurm_scripts)#
    print("> Number of BAM files :") ## has to be indel realigned bam files
    print([".bam files",len([f for f in os.listdir(INIT_DIR_MUT) if(f.endswith(".bam"))])])
    time.sleep(1)
    print(">Proceeding with MUTECT2:")
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
    for k,v in pairs.items():
        print("> Writing Slurm for %s and %s"%(k,v))
        time.sleep(1)
        #start_time=time.time()
        name = "/".join([k.split(".")[0], v.split(".")[0]])
        slurm=open("%s/%s.slurm"%(slurm_scripts,name),"w")
        if slurm_options=="y":
            options=["#!/bin/bash","\n","#SBATCH --mail-user=%s"%MAIL,"\n",
             "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%name,"\n","#SBATCH --time=80:00:00",
             "\n","#SBATCH --nodes=1","\n","#SBATCH --cpus-per-task=%s"%cpus,"\n","#SBATCH --output=%s.out"%name,
             "\n","#SBATCH --mem=16gb","\n","#SBATCH -A %s"%ACCOUNT,"\n","#SBATCH --workdir=%s"%slurm_scripts,"\n"]
        else:
            options=["#!/bin/bash","\n","#SBATCH --mail-user=%s"%MAIL,"\n",
             "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%name,"\n","#SBATCH --time=%s"%t,
             "\n","#SBATCH --nodes=%s"%nodes,"\n","#SBATCH --cpus-per-task=%s"%cpus,"\n","#SBATCH --output=%s.out"%name,
             "\n","#SBATCH --mem=%sgb"%mem,"\n","#SBATCH -A %s"%ACCOUNT,"\n","#SBATCH --workdir=%s"%slurm_scripts,"\n"]
        out_vcf = name + ".vcf"
        vcf_file = "/".join([OUT_DIR_MUT, out_vcf])
        for j in range(len(options)):
             #print(options[j])
            slurm.write(options[j])
        slurm.write("\n")
        slurm.write("%s -Xmx16g -jar %s -T MuTect2 -R %s -I:normal %s -I:tumor %s --dbsnp %s --output_mode EMIT_VARIANTS_ONLY -o %s"%(JAVA,GATK,tool,REF_GENOME,"/".join([INIT_DIR_MUT,v]), "/".join([INIT_DIR_MUT,k]),vcf,vcf_file))
        slurm.write("\n")
        slurm.write("wait")
        time.sleep(1)
    
        slurm.close()
        print(">EXECUTION OF MUTECT2 SCRIPTS")
        time.sleep(1)
    # CHUNK 02: EXECUTION OF SLURM SCRIPTS
        slurm_list=[f for f in os.listdir(slurm_scripts) if f.endswith(".slurm")]
        print("Total scripts %s"%len(slurm_list))  ## output directory same as out_dir
        
        for i in slurm_list:
            s1= "/".join([slurm_scripts,i]) 
            sub.call("sbatch %s"%s1, shell=True)
        
        slurm_out_list=[f for f in os.listdir(slurm_scripts) if f.endswith("%s.out"%i.split(".")[0])]
        print("Total_out_scripts = %s"%len(slurm_out_list))
        
        # CHUNK 03:CODES AFTER COMPLETION OF JOBS
        fp.GetErrors(slurm_scripts) ## NEED TO CHECK GETERRORS AND SUCCESS FUNCTION.
        time.sleep(180)

    out_log= open("/".join([slurm_scripts,"Mutect2.log"]), "w") ## making the combined file of last lines
    c=0
    out = [f for f in os.listdir(slurm_scripts) if f.endswith(".out")]
    #total_read_list = []
    #reads_left=[]
    final= {}
    ## Mutect stat_out, vcf_out and slurm_out are in the format TN05.stats.txt, TN05.vcf and TN05.out 
    ## BUT there will be 2 indel realigned files for every sample
    ## Therefore this creating a final dictionary with KEY=out_file(TN05.out) and VALUE = [TN05_T.bam, TN05_N.bam]
    ## easy to write the size and for reassurance of path
    for k in out:
        f = k.split("/")
        f1 = f[0]+".bam"
        f2 = f[1]+".bam"
        final[k]=[f1,f2]
                
    for k,v in final.items():
            #for i in init_bam_files:
            #  if i.split(".")[0]==k.split(".")[0]:
                  c+=1
                  init_path_1= "/".join([INIT_FILES_PATH,v[0]])  # joining path to script CANNOT TAKE INTIAL PATH CZ THE THERE WILL BE TWO FILES
                  init_bam_size_1 =((os.path.getsize(init_path_1))/1000000000)
                  
                  init_path_2= "/".join([INIT_FILES_PATH,v[1]])  # joining path to script CANNOT TAKE INTIAL PATH CZ THE THERE WILL BE TWO FILES
                  init_bam_size_2=((os.path.getsize(init_path_2))/1000000000)
                  
                  recent_output_vcf= "/".join([OUT_DIR_MUT,k.split('.')[0]+".vcf"])
                  out_vcf_size=((os.path.getsize(recent_output_vcf))/1000000000)
                  
                  #recent_output_stat= "/".join([OUT_DIR_MUT,k.split('.')[0]+".stats.txt"])
                  #out_stat_size=((os.path.getsize(recent_output_stat))/1000000000)
                 
                  out_log.write("%s. %s\n"%(c,k))
                  #out_log.write("\n")
                  out_log.write("INPUT_1: %s - %s gb\n"%(init_path_1,init_bam_size_1))
                  out_log.write("INPUT_2: %s - %s gb\n"%(init_path_2,init_bam_size_2))
                 # out_log.write("\n")
                  out_log.write("OUTPUT_VCF: %s - %s gb\n"%(recent_output_vcf,out_vcf_size))
                  #out_log.write("OUTPUT_STAT.TXT: %s - %s gb\n"%(recent_output_stat,out_stat_size))
                  #out_log.write("\n")
                  slurm_out_file= "/".join([slurm_scripts,k])
                  p3= sub.Popen(['grep','-i','Program Args',slurm_out_file],stdout=sub.PIPE) 
                  soutput,sinput=p3.communicate()
                  out_log.write("%s\n"%(" ".join(((soutput.decode("utf-8")).split())[4:])))## writing the command used
                  lines=fp.GetLastNLines(20,slurm_out_file)
                  for i in lines[-8:-5]:
                      out_log.write("%s"%(" ".join([j for j in i])))
                      out_log.write("\n")

                  reads_left.append([k,v,int(lines[-6][-4])-int(lines[-6][4])])
                  total_read_list.append([k,v,lines[-6][-4]])
                  out_log.write("\n")
             # else:
              #    continue
    out_log.close()
    df_total_reads = pd.DataFrame.from_records(total_read_list, columns=["Sample1", "Sample2","Total_Reads"])
    df_reads_left = pd.DataFrame.from_records(reads_left,columns = ["Sample1", "Sample2","Reads_Remaining"])
    merged = pd.merge(df_total_reads,df_reads_left,how="inner", on="Sample1")
    merged.to_csv("/".join([slurm_scripts,"Reads_per_Sample_after_mutect2.csv"], sep="\t", index=False))

def mutect2_with_intervals(PROJECT,INIT, pairs):
    PROC = "Mutect2" #
    INIT_DIR_MUT= INIT_FILES_PATH# initial files need to be indel realigner
    OUT_DIR_MUT="/".join([PROJECT,PROC])#
    if not os.path.exists(OUT_DIR_MUT):
        os.mkdir(OUT_DIR_MUT)#
    slurm_scripts = "/".join([OUT_DIR_MUT, "Mutect2_scripts"])#
    if not os.path.exists(slurm_scripts):
        os.mkdir(slurm_scripts)#
    print("> Number of BAM files :") ## has to be indel realigned bam files
    print([".bam files",len([f for f in os.listdir(INIT_DIR_MUT) if(f.endswith(".bam"))])])
    time.sleep(1)
    print(">Proceeding with MUTECT2:")
    print(">ENTER LOCATION OF PER CHROMOSOME TARGET INTERVAL FILE")
    target_intervals = input()
    time.sleep(1)
    print("> ENTER TOOL shown here: MuTect2") ## need to see whether to add it or not cz unlike gatk toolsonly one is required
    tool=input()
    time.sleep(1)
    intv = sorted(os.listdir(target_intervals))
    for i in sorted(intv):
        print("> Writing Slurm for %s"%i.split(".")[0])
        time.sleep(1)
        #start_time=time.time()
        #slurm=open("%s/%s.slurm"%(slurm_scripts,k.split(sep)[0]),"w")
        for k,v in sorted(pairs.items()):
            intv_sam = ".".join([k.split(sep)[0],i.split(".")[0]])
            slurm=open("%s/%s.slurm"%(slurm_scripts,intv_sam),"w")
            options=["#!/bin/bash","\n","#SBATCH --mail-user=%s"%MAIL,"\n",
                 "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%intv_sam,"\n","#SBATCH --time=30:00:00",
                 "\n","#SBATCH --nodes=1","\n","#SBATCH --cpus-per-task=1","\n","#SBATCH --output=%s.out"%intv_sam,
                 "\n","#SBATCH --mem=8gb","\n","#SBATCH -A %s"%ACCOUNT,"\n","#SBATCH --workdir=%s"%slurm_scripts,"\n", "%s"%g_mod, "\n"]
            out_vcf = intv_sam + ".vcf"
            vcf_file = "/".join([OUT_DIR_MUT, out_vcf])
            for j in range(len(options)):
             #print(options[j])
                slurm.write(options[j])
            slurm.write("\n")
            slurm.write("%s -jar %s -T %s -R %s -I:normal %s -I:tumor %s --dbsnp %s -L %s --output_mode EMIT_VARIANTS_ONLY --tumor_lod 10 -o %s"%(JAVA,GATK,tool,REF_GENOME,"/".join([INIT_DIR_MUT,v]), "/".join([INIT_DIR_MUT,k]),vcf,"/".join([target_intervals,i]),vcf_file))
            slurm.write("\n")
            slurm.write("wait")
            time.sleep(1)
    
        slurm.close()
    print(">EXECUTION OF MUTECT SCRIPTS")
    time.sleep(1)
    # CHUNK 02: EXECUTION OF SLURM SCRIPTS
    slurm_list=sorted([f for f in os.listdir(slurm_scripts) if f.endswith(".slurm")])
    print("Total scripts %s"%len(slurm_list))  ## output directory same as out_dir
    print("Slurm_list ending .slurm")
    print(slurm_list)
    list_of_samples = []
    if len(sorted(slurm_list))>80:
        remain = len(slurm_list)-(len(slurm_list)%80)
        divide = int(len(slurm_list)/80)
        for i in range(divide):
            lis = []
            first = i*80
            end = (i+1)*80
            list_1 = slurm_list[first:end]
            for j in list_1:
                lis.append( ".".join(j.split(".")[0:2])+".out")
                s1= "/".join([slurm_scripts,j]) 
                sub.call("sbatch %s"%s1, shell=True)
                
            #slurm_out_list=[f for f in os.listdir(slurm_scripts) if f.endswith(".out")]
            list_of_samples.append(list_1)   
            print("Total_out_scripts = %s"%len(lis))
            fp.GetErrors_mutect2_interval(slurm_scripts, list_1) ## NEED TO CHECK GETERRORS AND SUCCESS FUNCTION.
            time.sleep(180)
        list_2 = slurm_list[remain:]
        for j in list_2:
                s1= "/".join([slurm_scripts,j]) 
                sub.call("sbatch %s"%s1, shell=True)
        
        slurm_out_list=[f for f in os.listdir(slurm_scripts) if f.endswith(".out")]
        print("Total_out_scripts = %s"%len(slurm_out_list))
        fp.GetErrors_mutect2_interval(slurm_scripts, list_2) ## NEED TO CHECK GETERRORS AND SUCCESS FUNCTION.
        time.sleep(180)

    out_log= open("/".join([slurm_scripts,"Mutect2.log"]), "w") ## making the combined file of last lines
    list_of_init_files = sorted(os.listdir(INIT_FILES_PATH))
    init_bam_files=sorted([f for f in list_of_init_files if f.endswith(".bam")]) ## initial bam files will be the same
    c=0
    out = [f for f in os.listdir(slurm_scripts) if f.endswith(".out")]
    #total_read_list = []
    #reads_left=[]
    final= {}
    ## Mutect stat_out, vcf_out and slurm_out are in the format TN05.stats.txt, TN05.vcf and TN05.out 
    ## BUT there will be 2 indel realigned files for every sample
    ## Therefore this creating a final dictionary with KEY=out_file(TN05.out) and VALUE = [TN05_T.bam, TN05_N.bam]
    ## easy to write the size and for reassurance of path
    for k in sorted(out):
        chr_int = ".".join(k.split('.')[0:2])
        list_indel_files = []
        for i in init_bam_files:
            if k.split('.')[0]== i.split('.')[0].split(sep)[0]:
                list_indel_files.append(i)
            else:
                final[chr_int]=list_indel_files
    print("END FINAL dictionary ")
    print(final)
    for k,v in final.items():
            #for i in init_bam_files:
            #  if i.split(".")[0]==k.split(".")[0]:
                  c+=1
                  init_path_1= "/".join([INIT_FILES_PATH,v[0]])  # joining path to script CANNOT TAKE INTIAL PATH CZ THE THERE WILL BE TWO FILES
                  init_bam_size_1 =((os.path.getsize(init_path_1))/1000000000)
                  
                  init_path_2= "/".join([INIT_FILES_PATH,v[1]])  # joining path to script CANNOT TAKE INTIAL PATH CZ THE THERE WILL BE TWO FILES
                  init_bam_size_2=((os.path.getsize(init_path_2))/1000000000)
                  
                  recent_output_vcf= "/".join([OUT_DIR_MUT,k+".vcf"])
                  out_vcf_size=((os.path.getsize(recent_output_vcf))/1000000000)
                  
                  #recent_output_stat= "/".join([OUT_DIR_MUT,k.split('.')[0]+".stats.txt"])
                  #out_stat_size=((os.path.getsize(recent_output_stat))/1000000000)
                 
                  out_log.write("%s. %s\n"%(c,k))
                  #out_log.write("\n")
                  out_log.write("INPUT_1: %s - %s gb\n"%(init_path_1,init_bam_size_1))
                  out_log.write("INPUT_2: %s - %s gb\n"%(init_path_2,init_bam_size_2))
                 # out_log.write("\n")
                  out_log.write("OUTPUT_VCF: %s - %s gb\n"%(recent_output_vcf,out_vcf_size))
                  #out_log.write("OUTPUT_STAT.TXT: %s - %s gb\n"%(recent_output_stat,out_stat_size))
                  #out_log.write("\n")
                  slurm_out_file= "/".join([slurm_scripts,k])
                  p3= sub.Popen(['grep','-i','Program Args',slurm_out_file],stdout=sub.PIPE) 
                  soutput,sinput=p3.communicate()
                  out_log.write("%s\n"%(" ".join(((soutput.decode("utf-8")).split())[4:])))## writing the command used
                  lines=fp.GetLastNLines(20,slurm_out_file)
                  for i in lines[-8:-5]:
                      out_log.write("%s"%(" ".join([j for j in i])))
                      out_log.write("\n")
                  #out_log.write("%s"%lines)
                  #for j in lines[-6]:
                  #print([c,i,lines[-6][-4:-1]])
                  #reads_left.append([c,k,int(lines[-6][-4])-int(lines[-6][4])])
                  #total_read_list.append([c,k,lines[-6][-4]])
                  #out_log.write("%s,%s,%s"%(c,i,lines[-6][-4:-1][0]))
                     #out_log.write("%s"%i[-4:-1])
                     #out_log.write("\n")
                  out_log.write("\n")
             # else:
              #    continue
    out_log.close()
    #df_total_reads = pd.DataFrame.from_records(total_read_list, columns=["Sample1", "Sample2","Total_Reads"])
    #df_reads_left = pd.DataFrame.from_records(reads_left,columns = ["Sample1", "Sample2","Reads_Remaining"])
    #merged = pd.merge(df_total_reads,df_reads_left,how="inner", on="Sample1")
    #merged.to_csv("/".join([slurm_scripts,"Reads_per_Sample_after_mutect2.csv"], sep="\t", index=False))
def mutect2_with_intervals_1(tumor,normal):
    PROC = "Mutect2" #
    sample = tumor.split("_")[0]
    INIT_DIR_MUT= INIT# initial files need to be indel realigner
    OUT_DIR_MUT="/".join([PROJECT,PROC])#
    out_sample = "/".join([OUT_DIR_MUT,tumor.split("_")[0]])
    if not os.path.exists(OUT_DIR_MUT):
        os.mkdir(OUT_DIR_MUT)#
    if not os.path.exists(out_sample):
        os.mkdir(out_sample)
    slurm_scripts = "/".join([out_sample, "Mutect2_scripts"])#
    if not os.path.exists(slurm_scripts):
        os.mkdir(slurm_scripts)#
    #print("> Number of BAM files :") ## has to be indel realigned bam files
    print("Tumor Sample : %s"%tumor)
    print("Normal Sample: %s"%normal)
    time.sleep(1)
    print(">Proceeding with MUTECT2:")
    time.sleep(1)
    #print("> ENTER TOOL shown here: MuTect2") ## need to see whether to add it or not cz unlike gatk toolsonly one is required
    #tool=input()
    #time.sleep(1)
    intv = sorted(os.listdir(target_intervals))
    for i in sorted(intv):
        if i.startswith(".panfs"):
            continue
        else:
            print("> Writing Slurm for %s"%i.split(".")[0])
            time.sleep(1)
            intv_sam = ".".join([tumor.split("_")[0],i.split(".")[0]])
            slurm=open("%s/%s.slurm"%(slurm_scripts,intv_sam),"w")
            options=["#!/bin/bash","\n","#SBATCH --mail-user=%s"%MAIL,"\n",
                 "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%intv_sam,"\n","#SBATCH --time=48:00:00",
                 "\n","#SBATCH --nodes=1","\n","#SBATCH --cpus-per-task=1","\n","#SBATCH --output=%s.out"%intv_sam,
                 "\n","#SBATCH --mem=16gb","\n","#SBATCH -A %s"%ACCOUNT,"\n","#SBATCH --workdir=%s"%slurm_scripts,"\n", "%s"%g_mod, "\n"]
            out_vcf = intv_sam + ".vcf"
            vcf_file = "/".join([out_sample, out_vcf])
            for optn in range(len(options)):
             #print(options[j])
                slurm.write(options[optn])
            slurm.write("\n")
            slurm.write("%s -Xmx8g -jar %s -T MuTect2 -R %s -I:normal %s -I:tumor %s --dbsnp %s -L %s --output_mode EMIT_VARIANTS_ONLY --tumor_lod 10 -o %s"%(JAVA,GATK,REF_GENOME,"/".join([INIT_DIR_MUT,normal]),"/".join([INIT_DIR_MUT,tumor]),vcf,"/".join([target_intervals,i]),vcf_file))
            slurm.write("\n")
            slurm.write("wait")
            time.sleep(1)
    
    slurm.close()
    print(">EXECUTION OF MUTECT SCRIPTS")
    time.sleep(1)
    # CHUNK 02: EXECUTION OF SLURM SCRIPTS
    slurm_list=sorted([f for f in os.listdir(slurm_scripts) if f.endswith(".slurm")])
    print("Total scripts %s"%len(slurm_list))  ## output directory same as out_dir
    print("Slurm_list ending .slurm")
    print(slurm_list)
    #list_of_samples = []
    for run in slurm_list:
        s1= "/".join([slurm_scripts,run]) 
        sub.call("sbatch %s"%s1, shell=True)
    fp.GetErrors_with_intervals(slurm_scripts,sample)
    time.sleep(180)

    #out_log= open("/".join([slurm_scripts,"Mutect2.log"]), "w") ## making the combined file of last lines
    #out_Script = sorted([f for f in os.listdir(slurm_scripts) if f.endswith(".out")])
    out_sample1 = sorted([f for f in os.listdir(out_sample) if f.endswith(".vcf")])
    m= ".".join(out_sample1[0].split(".")[0:2])+".txt"
    sub.call("%s -Xmx4g -jar /mnt/pan/Data16/Tools/snpEff/SnpSift.jar extractFields %s CHROM POS ID REF ALT QUAL FILTER GEN[0].GT GEN[0].AD[0] GEN[0].AD[1] GEN[0].AF GEN[1].GT GEN[1].AD[0] GEN[1].AD[1] GEN[1].AF > %s"%(JAVA, "/".join([out_sample,out_sample1[0]]), "/".join([out_sample,m])),shell=True)
    df = pd.read_table("/".join([out_sample,m]), dtype=None, sep="\t",index_col=False)
    for z in sorted(out_sample1[1:]):
        m= ".".join(z.split(".")[0:2])+".txt"
        sub.call("%s -Xmx4g -jar /mnt/pan/Data16/Tools/snpEff/SnpSift.jar extractFields %s CHROM POS ID REF ALT QUAL FILTER GEN[0].GT GEN[0].AD[0] GEN[0].AD[1] GEN[0].AF GEN[1].GT GEN[1].AD[0] GEN[1].AD[1] GEN[1].AF > %s"%(JAVA, "/".join([out_sample,z]), "/".join([out_sample,m])),shell=True)
        df1 = pd.read_table("/".join([out_sample,m]), dtype=None, sep="\t",index_col=False)
        df = pd.concat([df,df1])
    val = df.values
    df1 = pd.DataFrame.from_records(val, columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "TUMOR_GT", "TUMOR_AD_REF","TUMOR_AD_ALT", "TUMOR_AF", "NORMAL_GT","NORMAL_AD_REF","NORMAL_AD_ALT", "NORMAL_AF"])
    df1.to_csv("/".join([out_sample,tumor.split("_")[0]+"_allvariants.txt"]), sep="\t", index=False)

    
def mpileup(PROJECT,INIT, pairs):
    PROC = "Mpileup" #
    INIT_DIR_MPILEUP= INIT_FILES_PATH # BASE RECALIBRATION
    OUT_DIR_MPILEUP="/".join([PROJECT,PROC])# 
    if not os.path.exists(OUT_DIR_MPILEUP): 
        os.mkdir(OUT_DIR_MPILEUP)# 
    #files_in_directory= os.listdir(INIT_DIR_MUT)
    slurm_scripts = "/".join([OUT_DIR_MPILEUP, "MPILEUP_SCRIPTS"])#
    if not os.path.exists(slurm_scripts):
        os.mkdir(slurm_scripts)#
    print("> Number of BAM files :")
    print([".bam files",len([f for f in os.listdir(INIT_DIR_MPILEUP) if(f.endswith(".bam"))])])
    time.sleep(1)
    print(">Proceeding with MPILEUP:")
    time.sleep(1)
    print("> ENTER TOOL shown here:mpileup ") ## need to see whether to add it or not cz unlike gatk toolsonly one is required
    tool=input()
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
    for k,v in pairs.items():
        print("> Writing Slurm for %s"%k)
        time.sleep(1)
        #start_time=time.time()
        slurm=open("%s/%s.slurm"%(slurm_scripts,k.split(sep)[0]),"w")
        if slurm_options=="y":
            options=["#!/bin/bash","\n","#SBATCH --mail-user=%s"%MAIL,"\n",
             "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%k.split(sep)[0],"\n","#SBATCH --time=80:00:00",
             "\n","#SBATCH --nodes=1","\n","#SBATCH --cpus-per-task=%s"%cpus,"\n","#SBATCH --output=%s.out"%k.split(sep)[0],
             "\n","#SBATCH --workdir=%s"%slurm_scripts,"\n","#SBATCH --mem=16gb","\n","#SBATCH -A %s"%ACCOUNT,"\n"]
        else:
            options=["#!/bin/bash","\n","#SBATCH --mail-user=%s"%MAIL,"\n",
             "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%k.split(sep)[0],"\n","#SBATCH --time=%s"%t,
             "\n","#SBATCH --nodes=%s"%nodes,"\n","#SBATCH --cpus-per-task=%s"%cpus,"\n","#SBATCH --output=%s.out"%k.split(sep)[0],
             "\n","#SBATCH --workdir=%s"%slurm_scripts,"\n","#SBATCH --mem=%sgb"%mem,"\n","#SBATCH -A %s"%ACCOUNT,"\n"]
        out_mpileup = k.split("_")[0]+".mpileup" ### need the extension option to be inputted.
        mpileup_file="/".join([OUT_DIR_MPILEUP,out_mpileup]) # out_dir for baserecalibrator and -BQSR for Print Reads
        for j in range(len(options)):
         #print(options[j])
            slurm.write(options[j])
        slurm.write("\n")
        slurm.write("%s mpileup -B -d 16000 -f %s %s %s > %s"%(samtools,tool,REF_GENOME,"/".join([INIT_DIR_MPILEUP,v]), "/".join([INIT_DIR_MPILEUP,k]),mpileup_file))
        slurm.write("\n")
        slurm.write("wait")
        time.sleep(1)

    slurm.close()
    print(">EXECUTION OF MPILEUP SCRIPTS")
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
    time.sleep(5)
    #Done, no_msg =fp.Success(slurm_scripts,slurm_out_list)
    #'''fp.forOutputLog(OUT_DIR_ALIGN, slurm_scripts,INIT_FASTQ_FILES, success_out, "align", ".sam")'''
    #''' THIS HAS TO BE MODIFIED '''
    #no_m=[k for k, v in no_msg.items()]
    #don=[k for k, v in Done.items()]
    print("length of completed samples %s"%len(os.listdir(OUT_DIR_MPILEUP)))
    #print("length of completed samples with no mesg %s"%len(no_m))
    return OUT_DIR_MPILEUP

def somatic(PROJECT,OUT_DIR_MPILEUP):
    PROC = "Somatic" #
    INIT_DIR_SOMATIC= INIT_FILES_PATH # MPILEUP files to be used 
    OUT_DIR_SOMATIC="/".join([PROJECT,PROC])# 
    if not os.path.exists(OUT_DIR_SOMATIC): 
        os.mkdir(OUT_DIR_SOMATIC)# 
    #files_in_directory= os.listdir(INIT_DIR_MUT)
    slurm_scripts = "/".join([OUT_DIR_SOMATIC, "somatic_SCRIPTS"])#
    if not os.path.exists(slurm_scripts):
        os.mkdir(slurm_scripts)#
    print("> Number of Mpileup files :")
    print([".mpileup files",len([f for f in os.listdir(INIT_DIR_SOMATIC) if(f.endswith(".mpileup"))])])
    time.sleep(1)
    print(">Proceeding with Varscan SOMATIC:")
    time.sleep(1)
    print(">ENTER full jar file location of Varscan:")
    Varscan = input()
    print("> ENTER TOOL shown here:somatic ") ## need to see whether to add it or not cz unlike gatk toolsonly one is required
    tool=input()
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
    files_in_directory = os.listdir(INIT_DIR_SOMATIC)
    for f in files_in_directory:
       if (f.endswith(".mpileup")):
        i=re.split("\.",f)[0]
        print("> Writing Slurm for %s"%i)
        time.sleep(1)
        slurm=open("%s/%s.slurm"%(slurm_scripts,i),"w")
        if slurm_options=="y":
            options=["#!/bin/bash","\n","#SBATCH --mail-user=%s"%MAIL,"\n",
             "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%i,"\n","#SBATCH --time=80:00:00",
             "\n","#SBATCH --nodes=1","\n","#SBATCH --cpus-per-task=%s"%cpus,"\n","#SBATCH --output=%s.out"%i,
             "\n","#SBATCH --workdir=%s"%slurm_scripts,"\n","#SBATCH --mem=16gb","\n","#SBATCH -A %s"%ACCOUNT,"\n"]
        else:
            options=["#!/bin/bash","\n","#SBATCH --mail-user=%s"%MAIL,"\n",
             "#SBATCH --mail-type=ALL","\n", "#SBATCH -J %s"%i,"\n","#SBATCH --time=%s"%t,
             "\n","#SBATCH --nodes=%s"%nodes,"\n","#SBATCH --cpus-per-task=%s"%cpus,"\n","#SBATCH --output=%s.out"%i,
             "\n","#SBATCH --workdir=%s"%slurm_scripts,"\n","#SBATCH --mem=%sgb"%mem,"\n","#SBATCH -A %s"%ACCOUNT,"\n"]
        out_somatic = i+"_output" ### need the extension option to be inputted.
        somatic_file="/".join([OUT_DIR_SOMATIC,out_somatic]) # out_dir for baserecalibrator and -BQSR for Print Reads
        for j in range(len(options)):
         #print(options[j])
            slurm.write(options[j])
        slurm.write("\n")
        slurm.write("%s -Xmx16g -jar %s somatic %s %s --mpileup 1"%(JAVA,Varscan,tool,"/".join([INIT_DIR_SOMATIC,f]),somatic_file))
        slurm.write("\n")
        slurm.write("wait")
        time.sleep(1)

    slurm.close()
    print(">EXECUTION OF Varscan SCRIPTS")
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
    '''time.sleep(10)
    Done, no_msg =fp.Success(slurm_scripts)
    no_m=[k for k, v in no_msg.items()]
    don=[k for k, v in Done.items()]
    print("length of completed samples %s"%len(don))
    print("length of completed samples with no mesg %s"%len(no_m))
    for i in no_m:
        if not i in don:
            don.append(i)
        else:
            continue'''
    don=[f for f in os.listdir(slurm_scripts) if f.endswith(".out")]
    out_log= open("/".join([slurm_scripts,"somatic.log"]), "w") ## making the combined file of last lines
    list_of_init_files = os.listdir(INIT_FILES_PATH)
    init_bam_files=[f for f in list_of_init_files if f.endswith(".mpileup")] ## initial bam files will be the same
    c=0
    #total_read_list = []
    #reads_left=[]
    final= []
    for k in sorted(don):
            print("writing log")
            for i in init_bam_files:
              if i.split(".")[0]==k.split(".")[0]:#for i in init_bam_files:
                  print([i.split(".")[0], k.split(".")[0]])
                  c+=1
                  init_path_1= "/".join([INIT_FILES_PATH,i])  # joining path to script CANNOT TAKE INTIAL PATH CZ THE THERE WILL BE TWO FILES
                  init_bam_size_1 =((os.path.getsize(init_path_1))/1000000000)
                  recent_output_snp= "/".join([OUT_DIR_SOMATIC,i.split(".")[0]+"_output.snp"])
                  out_snp_size=((os.path.getsize(recent_output_snp))/1000000)
                  recent_output_indel= "/".join([OUT_DIR_SOMATIC,i.split(".")[0]+"_output.indel"])
                  out_indel_size=((os.path.getsize(recent_output_indel))/1000000)
                 
                  out_log.write("%s. %s\n"%(c,i))
                 
                  out_log.write("INPUT_1: %s - %s gb\n"%(init_path_1,init_bam_size_1))
                  out_log.write("OUTPUT_snp: %s - %s mb\n"%(recent_output_snp,out_snp_size))
                  out_log.write("OUTPUT_indel: %s - %s mb\n"%(recent_output_indel,out_indel_size))
                  
                  slurm_out_file= "/".join([slurm_scripts,i.split(".")[0]+'.out'])
                  slurm_in_file = "/".join([slurm_scripts,i.split('.')[0]+'.slurm'])
                  lines_in = fp.GetLastNLines(3,slurm_in_file)
                  out_log.write("Program_Args: %s\n"%" ".join(lines_in[1]))
                
                  lines_out=fp.GetLastNLines(22,slurm_out_file)
                  out_log.write("%s\n"%(" ".join(lines_out[0])))
                  out_log.write("%s\n"%(" ".join(lines_out[3])))
                  out_log.write("%s\n"%(" ".join(lines_out[7])))
                  out_log.write("%s\n"%(" ".join(lines_out[9])))
                  
                  list_1=[k.split(".")[0]]
                  for i in lines_out[-5:]:
                      list_1.append(i[0])
                  final.append(list_1)
                  out_log.write("\n")
          
              else:
                  continue
               
    out_log.close()
    df = pd.DataFrame.from_records(final, columns =["sample","germline", "loh", "somatic", "unknown", "variant"])
    df.to_csv("/".join([slurm_scripts,"varscan_result.csv"], index= False))
