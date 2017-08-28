#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 14:44:46 2017

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

def knowJobStatus(uid):
    '''
    FUNCTION 03
    Name:           get the status of slurm jobs
    Description: return a list of squeue output and dictionary with job id as key and job name as value
    Keyword argument:
    uid -- user id 
    '''
    p=sub.Popen(['squeue','-u','%s'%uid], stdout=sub.PIPE)
    soutput,sinput=p.communicate()
    list_stat=[]
    dict_stat={}
    for line in soutput.decode("utf-8").split("\n"):
        list_stat.append(line.split())
    job_stat=list_stat[1:-1] # removing the heading and  last empty list as a result of split
    for i in range(len(job_stat)):
        dict_stat[job_stat[i][0]]=job_stat[i][2] # getting the keys = job id and values =job name so as to scancel these jobs in case of error
    return [job_stat, dict_stat]

def scancel(dict_stat):
        '''
        FUNCTION 04
        Description : Cancel's the jobs if user requests system exit on an error
        Keyword_Argument: 
        dict_stat -- dictionary with KEY as JOBID's and VALUES as JOB names (samples)
        '''
    
        list_1 = []
    #if len(err_samples)!=0:
        for k, v in dict_stat.items():
            list_1.append(k)
        mx = max(list_1)
        mn = min(list_1)
        sub.call("for i in {%s..%s}; do scancel $i;done"%(mn,mx), shell=True)
        
def GetErrors(script_dir,slurm_out_list):
    '''
    FUNCTION 05
    Name:           GetErrors
    Description:    Finds for error strings in all slurm.out files
    Output:         Until jobs are running it will check and give errors,
    As soon as it encounters an error in any output file, program is halted,
    If there is no error prog continues and returns a list of samples completed.
    Keyword argument:   File_path, slurm.out list & job_stat
    script_dir :         Is a path to slurm.out files
    slurm_out_list:     A list of slurm.out files
    '''
    job_stat=knowJobStatus("jxk906")[0] ## one time calling of job status to initiate loop
    dict_stat=knowJobStatus("jxk906")[1]
    while (len(job_stat)>0):
        print(job_stat)
        list_err=["slurmstepd", "error", "exception"]
        err_samples=[]
        ErrorsCalled = {}
        print("> LOOKING INTO ERRORS")
        for i in slurm_out_list:
          path="/".join([script_dir,i])
          Err_called=[]
          for j in list_err:
                 p=sub.Popen(['grep','-i',j,path], stdout=sub.PIPE)
                 soutput,sinput=p.communicate()
                 if soutput==b'':
                     continue
                 else:
                     Err_called.append(soutput)
                     err_samples.append(i)
                     
          ErrorsCalled[i]=Err_called
        time.sleep(10)
        if len(err_samples)!=0: 
            print("Error Found In Some Samples: Execution-HALTED")
            print(ErrorsCalled)
            time.sleep(3)
            print("Want to end jobs because of error? (y/n)")
            job_cancel = input()
            if job_cancel=="y":
                sys.exit() #just put up a user response to whether exit system or not  better control
                scancel(dict_stat)
            else:
                continue
        else:
            job_stat=knowJobStatus("jxk906")[0] ## calling job_status to continous monitor and keep loop updated 
            dict_stat=knowJobStatus("jxk906")[1]
            print("> NO ERROR OCUURED YET")
            #continue
    
    else:
          print("> Number of samples without error %s"%slurm_out_list)
          
    return dict_stat, err_samples


   
def Success(script_dir,slurm_out_list):
        '''
        FUNCTION 06
        Description : Search done and success mesgs in each samples output and provide
                      dictionary of samples
        Keyword Argument:
         script_dir -- location where scripts are stored
         slurm_out_list -- list of slurm output files
        '''
        list_success=["done", "main", "called"]
        Done={}
        no_msg={}
        print("> LOOKING INTO SAMPLES SUCCESSFULLY DONE")
        for k in slurm_out_list:
              #print(k)
              path="/".join([script_dir,k])
              Success_called=[]
              no_err=[]
              for m in list_success:
                   #print(k[0])
                   p1=sub.Popen(['grep','-i',m,path], stdout=sub.PIPE)
                   soutput,sinput=p1.communicate()
                   
                   if soutput==b'':
                       no_err.append("Neither Error nor Done message found")
                   else:
                       Success_called.append(soutput)

              print("> APPENDING SAMPLES THAT ARE DONE WITHOUT ANY ERROR")
              Done[k]=Success_called
              no_msg[k]=no_err
        return Done, no_msg

def GetLastNLines(n,script_dir):
    '''
    FUNCTION# 05
    Name:           Get LastNLines
    Description:        Gets last n lines using Unix tail
    Output:         returns last n lines of a file
    Keyword argument:
    n -- number of last lines to return
    filename -- Name of the out file you need to tail into. 
    '''
    p=sub.Popen(['tail','-n',str(n),script_dir], stdout=sub.PIPE)
    soutput,sinput=p.communicate()
    list_lines=[]
    for line in soutput.decode("utf-8").split("\n"):
        list_lines.append(line.split())
    return list_lines[:-1]

def forOutputLog(sample_dir, script_dir, initial_path, Sucess_out, out_file_name, ext_of_recent_output):
        '''
       FUNCTION# 06
       Name:   forOutputLog
       Description:  returns number of reads left and total reads reviewed after subtracting, Also builds log file from output scripts along with initial and final file size and location
       Output:  log file with few lines decribing number of reads successfully added to the bam or files, will have location of initial and current files,
       along with size(MB)
       Keyword Argument:
           initial_path: path o fprevios directory _files
           sample_dir:  Directory where new files are created
           script_dir:  path to .out files
           out_file_name: combined log file name
           ext_of_recent_output: file extension in initial_path or recent output files.
           success_out: is the output of success 
        '''
        out= Sucess_out
        don = out[0]## getting list of errors
        no_m=out[1]
        init_bam_files=[]
        print("Writing  Log with last few lines per sample")
         ## for combining samples already done.
        no_msg=[k for k, v in no_m.items()]
        done=[k for k, v in don.items()]
        for i in no_msg:
            if not i in done:
               done.append(i)
            else:
                continue
        out_log= open("/".join([script_dir,"%s.log"%out_file_name]), "w") ## making the combined file of last lines
        list_of_init_files = os.listdir(initial_path)
        init_bam_files=[f for f in list_of_init_files if f.endswith(".bam")] ## initial bam files will be the same
        c=0
        total_read_list = []
        reads_left=[]
        for k in done:
            for i in init_bam_files:
              if i.split(".")[0]==k.split(".")[0]:
                  c+=1
                  init_path= "/".join([initial_path,i])  # joining path to script
                  init_bam_size =((os.path.getsize(init_path))/1000000000)
                  recent_output = "/".join([sample_dir,re.split("\.", i)[0]+"%s"%ext_of_recent_output])
                  out_file_size=((os.path.getsize(recent_output))/1000000000)
                  out_log.write("%s %s\n"%(c,i))
                  #out_log.write("\n")
                  out_log.write("INPUT: %s - %s gb\n"%(init_path,init_bam_size))
                 # out_log.write("\n")
                  out_log.write("OUTPUT: %s - %s gb\n"%(recent_output,out_file_size))
                  #out_log.write("\n")
                  slurm_out_file= "/".join([script_dir,k])
                  p3= sub.Popen(['grep','-i','Program Args',slurm_out_file],stdout=sub.PIPE)
                  soutput,sinput=p3.communicate()
                  out_log.write("%s\n"%(" ".join(((soutput.decode("utf-8")).split())[4:])))
                  lines=GetLastNLines(20,slurm_out_file)
                  for i in lines[-8:-1]:
                      out_log.write("%s"%(" ".join([j for j in i])))
                      out_log.write("\n")
                  #out_log.write("%s"%lines)
                  #for j in lines[-6]:
                  print([c,i,lines[-6][-4:-1]])
                  reads_left.append([c,i,lines[-6][-4]-lines[-6][4]])
                  total_read_list.append([c,i,lines[-6][-4]])
                  #out_log.write("%s,%s,%s"%(c,i,lines[-6][-4:-1][0]))
                     #out_log.write("%s"%i[-4:-1])
                     #out_log.write("\n")
                  out_log.write("\n")
              else:
                  continue
        out_log.close()
        df_total_reads = pd.DataFrame.from_records(total_read_list, columns=["sample", "total_reads_%s"%out_file_name], index=False)
        df_reads_left = pd.DataFrame.from_records(reads_left, columns=["samples", "reads_left_%s"%out_file_name], index=False)
        return df_reads_left, df_total_reads

''' we cannot deploy this function for target creator and recal grp file in indel realigner and base reclibration steps '''
''' ALSO not for mutect and varscan'''

''' need to provide two variables for the assignment'''