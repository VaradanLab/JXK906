# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
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
import FUNCTIONS_FOR_PIPELINE as fp
import NGSDATA_Processing_Pipeline as pl


print("RUN PROCESSES(a/s/n):")
ans = input()
if ans == "a":
    PROJECT, REF_GENOME, vcf, JAVA, b_mod, BWA, g_mod, GATK, PICARD, FASTQ, MAIL, ACCOUNT = pl.AskLocations()
    fastq = []
    bam = []
    sam = []
    for f in os.listdir(FASTQ):
        if (f.endswith(".fastq")):
            fastq.append(f)
            number_fastqfiles = len(fastq)
            actual_files_for_processing = int((number_fastqfiles)/2)
        elif (f.endswith(".sam")):
            sam.append(f)
        elif (f.endswith(".bam")):
            bam.append(f)
        else:
            continue
        
    list_proc = {}
    df_dict = {}
    ## 1 ALIGNMENT
    alignment=pl.alignment(PROJECT, FASTQ)
    time.sleep(5)
    list_proc[0]=alignment
    fp.Log_ASMRG(list_proc[0][1], list_proc[0][0],list_proc[0][2],"align", ".sam")
    data_frame_align = fp.samtools(list_proc[0][1],"align")
    df_dict[0]=data_frame_align
    
    ## 2 SORTING
    sortsam =pl.samsort(PROJECT, alignment[1])
    #if len(sortsam[1])==len(sam) :
    time.sleep(5)
    list_proc[1]=sortsam
    fp.Log_ASMRG(list_proc[1][1], list_proc[1][0],list_proc[1][2],"sortsam", ".bam")
    data_frame_samsort = fp.samtools(list_proc[1][1],"sortsam")
    df_dict[1]=data_frame_samsort
    
    ## 3 ADDREPLACEREADGROUP
    addrep = pl.addreplacegrp(PROJECT, sortsam[1])
    #if len(addrep[1])==len(bam):
    time.sleep(5)
    list_proc[2]=addrep
    fp.Log_ASMRG(list_proc[2][1], list_proc[2][0],list_proc[2][2],"add_rep", ".bam")
    data_frame_addrepgrp = fp.samtools(list_proc[2][1],"add_rep")
    df_dict[2]=data_frame_addrepgrp
    
    ## 4 MARKDUPLICATES
    markdup = pl.markduplicates(PROJECT, addrep[1])
    #if len(markdup[1])==len(bam):
    time.sleep(5)
    list_proc[3]=markdup
    fp.Log_ASMRG(list_proc[3][1], list_proc[3][0],list_proc[3][2],"markdup", ".bam")
    data_frame_markdup = fp.samtools(list_proc[3][1],"markdup")
    df_dict[3]=data_frame_markdup
    
    print("Remove SortSam file folder ? (y/n)")
    rem_sort = input()
    if rem_sort=="y":
        sub.call("rm -r %s"%sortsam[2],shell=True)
    
    ## 5 BASERECALIBRATION
    baserecalib = pl.baserecalibration(PROJECT, markdup[1])
    #if len(baserecalib[1])==len(bam):
    time.sleep(5)
    list_proc[4]=baserecalib
    fp.Log_BIRG(list_proc[4][1], list_proc[4][0],list_proc[4][2],"baserecalib", ".bam")
    data_frame_baserecalib = fp.samtools(list_proc[4][1],"baserecalib")
    df_dict[4]=data_frame_baserecalib
     
    print("Remove AddreplaceReadgroup file folder ? (y/n)")
    rem_add = input()
    if rem_add=="y":
        sub.call("rm -r %s"%addrep[1],shell=True)
    
    ## 6 INDEL-REALIGNMENT
    indelrealign = pl.indelrealignment(PROJECT, baserecalib[1])
    #if len(indelrealign[1])==len(bam):
    time.sleep(5)
    list_proc[5]=indelrealign
    fp.Log_BIRG(list_proc[5][1], list_proc[5][0],list_proc[5][2],"indelrealign", ".bam")
    data_frame_indelrealign = fp.samtools(list_proc[5][1],"indelrealign")
    df_dict[5]=data_frame_indelrealign
    
    print("Remove Markduplicates file folder ? (y/n)")
    rem_dup = input()
    if rem_dup=="y":
        sub.call("rm -r %s"%markdup[1],shell=True)
    
    print("Remove Baserecalibration file folder ? (y/n)")
    rem_baserecalib = input()
    if rem_baserecalib=="y":
        sub.call("rm -r %s"%baserecalib[1],shell=True)
    
    
    keys = sorted(list(df_dict.keys()))
    if len(keys)>2:
        df_total_merged = pd.merge(df_dict[keys[0]][0], df_dict[keys[1]][0], how="inner", on= "sample")
        df_mapped_merged = pd.merge(df_dict[keys[0]][1], df_dict[keys[1]][1], how="inner", on= "sample")
        for i in range(2,len(keys)):
                df_total_merged = pd.merge(df_total_merged, df_dict[keys[i]][0], how="inner", on= "sample")
                df_mapped_merged = pd.merge(df_mapped_merged, df_dict[keys[i]][1], how="inner", on= "sample")
    df_total_merged.to_csv("/".join([PROJECT, "TOTAL_READS_PER_SAMPLE_PER_PROCESS.csv"]), index = False)
    df_mapped_merged.to_csv("/".join([PROJECT, "MAPPED_READS_PER_SAMPLE_PER_PROCESS.csv"]), index = False)
  
elif ans == "s":
    PROJECT, REF_GENOME, vcf, JAVA, b_mod, BWA, g_mod, GATK, PICARD, FASTQ, MAIL, ACCOUNT = pl.AskLocations()
    fastq = []
    bam = []
    sam = []
    for f in os.listdir(FASTQ):
        if (f.endswith(".fastq")):
            fastq.append(f)
            number_fastqfiles = len(fastq)
            actual_files_for_processing = int((number_fastqfiles)/2)
        elif (f.endswith(".sam")):
            sam.append(f)
        elif (f.endswith(".bam")):
            bam.append(f)
        else:
            continue
    list_proc = {} 
    print("Proceed Alignent y/n")
    a = input()
    if a=="y":
        print("Require intial file path")
        init = input()
        alignment = pl.alignment(PROJECT, init)
        #if len(alignment[1])==len(actual_files_for_processing):
        list_proc[0]=alignment
    else:
        list_proc[0]="0"
    print("Proceed Sorting of sam files y/n")
    s = input()
    if s=="y":
        print("Require intial file path")
        init = input()
        sortsam =pl.samsort(PROJECT, init)
        #if len(sortsam[1])==len(sam) :
        list_proc[1]=sortsam
    else:
        list_proc[1]="0"
    print("Proceed Addreplacegrp in sorted bams y/n")
    ad = input()
    if ad=="y":
        print("Require intial file path")
        init = input()
        addrep = pl.addreplacegrp(PROJECT, init)
        #if len(addrep[1])==len(bam):
        list_proc[2]=addrep
    else:
        list_proc[2]="0"
    print("Proceed MarkDuplicates y/n")
    md = input()
    if md=="y":
        print("Require intial file path")
        init = input()
        markdup = pl.markduplicates(PROJECT, init)
        #if len(markdup[1])==len(bam):
        list_proc[3]=markdup
    else:
        list_proc[3]="0"
    print("Proceed base recalibration y/n")
    br = input()
    if br=="y":
        print("Require intial file path")
        init = input()
        baserecalib = pl.baserecalibration(PROJECT, init)
        #if len(baserecalib[1])==len(bam):
        list_proc[4]=baserecalib
    else:
        list_proc[4]="0"
    print("Proceed with indelrealignment y/n")
    ir = input()
    if ir=="y":
        print("Require intial file path")
        init = input()
        indelrealign = pl.indelrealignment(PROJECT, init)
        #if len(indelrealign[1])==len(bam):
        list_proc[5]=indelrealign
    else:
        list_proc[5]="0"
    
    dict_proc = {0:["align", ".sam"], 1:["sortsam", ".bam"], 2:["add_rep", ".bam"], 3:["markdup", ".bam"], 4:["baserecalib", ".bam"], 5:["indelrealign", ".bam"]}
    df_dict = {}
    for i in sorted(list_proc.keys()):
        if not list_proc[i]=="0":
            slurm_scripts = list_proc[i][0]
            outdir = list_proc[i][1]
            init_files = list_proc[i][2]
            if i==0 or i==1 or i==2 or i==3:
                #success_out=fp.Success(slurm_scripts)
                fp.Log_ASMRG(outdir, slurm_scripts,init_files,dict_proc[i][0], dict_proc[i][1])
                data_frame = fp.samtools(outdir,dict_proc[i][0])
                df_dict[i]=data_frame
            elif i==4 or i==5:
                #success_out=fp.Success(slurm_scripts)
                fp.Log_BIRG(outdir, slurm_scripts,init_files, dict_proc[i][0], dict_proc[i][1])
                data_frame = fp.samtools(outdir,dict_proc[i][0])
                df_dict[i]=data_frame
        else:
            continue
    keys = sorted(list(df_dict.keys()))
    #print(keys)
    #print(df_dict)
    if len(keys)>2:
        df_total_merged = pd.merge(df_dict[keys[0]][0], df_dict[keys[1]][0], how="inner", on= "sample")
        df_mapped_merged = pd.merge(df_dict[keys[0]][1], df_dict[keys[1]][1], how="inner", on= "sample")
        for i in range(2,len(keys)):
                df_total_merged = pd.merge(df_total_merged, df_dict[keys[i]][0], how="inner", on= "sample")
                df_mapped_merged = pd.merge(df_mapped_merged, df_dict[keys[i]][1], how="inner", on= "sample")
    elif len(keys)==2:
       df_total_merged = pd.merge(df_dict[keys[0]][0], df_dict[keys[1]][0], how="inner", on= "sample")
       df_mapped_merged = pd.merge(df_dict[keys[0]][1], df_dict[keys[1]][1], how="inner", on= "sample") 
    else:
        l = len(keys)
        c= keys[0:l]
        df_total_merged = df_dict[c[0]][0]
        df_mapped_merged = df_dict[c[0]][1]
    df_total_merged.to_csv("/".join([PROJECT, "TOTAL_READS_PER_SAMPLE_PER_PROCESS.csv"]), index = False)
    df_mapped_merged.to_csv("/".join([PROJECT, "MAPPED_READS_PER_SAMPLE_PER_PROCESS.csv"]), index = False)

else:
    print("Proceeding with HSMetrics")
    PROJECT, REF_GENOME, vcf, JAVA, b_mod, BWA, g_mod, GATK, PICARD, INIT = pl.AskLocations()
    hsmetrics = pl.hsmetrics(PROJECT,INIT)
    

    
    
