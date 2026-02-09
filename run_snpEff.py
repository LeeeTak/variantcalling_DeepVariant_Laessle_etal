#!/usr/bin/python
import subprocess
import os,sys
import glob
#first build database: 
#newest version: java -Xmx1G -jar /home/tlee/anaconda3/envs/variantcall/share/snpeff-5.2-2/snpEff.jar build -gff3 araport11
#subprocess.run("mkdir -p "+tool,shell=True)
indir = "/netscratch/dep_psl/grp_rgo/taklee/Henriette_Laessle/deepvariant/repeatfiltered"
infiles=glob.glob(f"{indir}/*_repeatfilt.vcf")
print(infiles)
#sampnames = {}
#namings = {"Ath":"Ath","microSly":"mSly","Hvu":"Hvu","Toothpick":"Tth","system":"Sys","Ancestral":"Ancestral"}
#for l in open("/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/sample_map_table.txt","r"):
#    la= l.strip().split("\t")
#    if la[1] in namings:
#        if la[1] == "Ancestral":
#            sampnames[la[0]] = "Ancestral"
#        else:
#            sampnames[la[0]] = namings[la[1]]+"_G"+la[2]+"_T"+la[3]

#print(sampnames)
for f in infiles:
    f=f.strip()
    fa=f.split("/")
    outn = fa[-1].split(".")
    name = f"{outn[0]}"
    outfile=f"{name}_snpEff.vcf"
    cmd = "snpEff -c snpEff.config araport11 "+f+" > "+outfile
    subprocess.run(cmd,shell=True)
    #print(cmd)
