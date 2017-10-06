#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 12:14:38 2017

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
import CALLER_PIPELINE as cp


global PROJECT, REF_GENOME, vcf, JAVA, g_mod, GATK, sam_mod, samtools, target_intervals,INIT, MAIL, ACCOUNT 

PROJECT = "/mnt/pan/Data16/jxk906/PROJECTS/TNBC"
REF_GENOME="/mnt/pan/Data16/jxk906/REFERENCE/hg19.fa"
vcf="/mnt/pan/Data16/jxk906/split_vcf/Merged_dbsnp.vcf"
JAVA="/mnt/pan/Data16/Tools/jdk1.8.0_91/bin/java"
GATK = "/mnt/pan/Data16/Tools/GATK.jar"
samtools = "samtools location"
target_intervals = "/mnt/pan/Data16/jxk906/target_file"
INIT = "/mnt/pan/Data16/jxk906/PROJECTS/TNBC/INDELREALIGNER"
MAIL = "jxk906@case.edu"
ACCOUNT = "lxh259"

pairs = cp.get_pairs_from_file()
cp.mutect2(PROJECT,INIT,pairs)