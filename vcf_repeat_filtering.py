#!/usr/bin/python

import subprocess
import os,sys
import re
import glob

suffix = ".PASS.GQ20.DP100.VAF90.deepvariant.vcf.gz" 
sys.path.insert(0, '/netscratch/dep_psl/grp_rgo/taklee/Henriette_Laessle/deepvariant')
from further_variant_filtering import statsgen,get_rank,merge_chkrep_vcfs,getflank,filter_Ancestral
#Qth = 30
#DPth = 10
workingpath = "."
hardfiltpath = "hardfiltered"
repeatfiltpath = "repeatfiltered"
firstinputfiles = glob.glob(f"{hardfiltpath}/*{suffix}")
firstinputs = ""
if len(firstinputfiles) > 0:
    firstinputs = f"{hardfiltpath}/*{suffix}"
    unzipping = "parallel bgzip -k -d {} ::: "+firstinputs
    print(1)
    subprocess.run(unzipping,shell=True)
else:
    suffix2 = suffix.replace(".vcf.gz",".vcf")
    firstinputs = f"{hardfiltpath}/*{suffix2}"
    zipping = "parallel bgzip -k {} ::: "+firstinputs
    print(2)
    subprocess.run(zipping,shell=True)

zipping = "parallel bgzip -k {} ::: "+hardfiltpath+"/*.vcf"
subprocess.run(zipping,shell=True)
indexing = "parallel tabix -p vcf {} ::: "+hardfiltpath+"/*.vcf.gz"
subprocess.run(indexing,shell=True)
merge_chkrep_vcfs(hardfiltpath,repeatfiltpath)
