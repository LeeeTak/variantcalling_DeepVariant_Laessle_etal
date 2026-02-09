#!/usr/bin/python
import os,sys
import glob
import numpy as np
import subprocess
reffile = "Col0.deepvariant.vcf"

files = ["d1.deepvariant.vcf"]
outdir = "hardfiltered"
for f in files:
    f = f.strip()
    samp = f.split(".")[0]
    ref = reffile.split(".")[0]
    subfile = f"{samp}subt{ref}.deepvariant.vcf"
    subcmd = f"bedtools subtract -A -a {f} -b {reffile} > {subfile}"
    print(subcmd)
    subprocess.run(subcmd,shell=True)
    outfn = outdir+"/"+subfile.replace(".deepvariant.",".PASS.GQ20.DP100.VAF90.deepvariant.")
    with open(subfile,'r') as vcf, open(outfn,'w+') as outf, open("vcfheader","r") as header:
        for h in header:
            outf.write(h)
        for l in vcf:
            l=l.strip()
            la=l.split("\t")
            print(l)
            if la[6] == "PASS":
                fields = la[8].split(":")
                vals = la[9].split(":")
                fvals = dict(zip(fields,vals))
                gq = float(fvals["GQ"])
                dp = float(fvals["DP"])
                vaf = np.mean([float(x) for x in fvals["VAF"].split(",")])
                if (gq > 20) and (dp > 90) and (vaf > 0.9):
                    outf.write(l+"\n")
