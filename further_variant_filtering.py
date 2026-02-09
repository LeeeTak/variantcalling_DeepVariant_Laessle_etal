#!/usr/bin/python

import subprocess
import os,sys
import re
import glob

def get_rank(sval,th):
    total=len(sval)
    i=0
    for x in sval:
        if x >= th:
            per=i/total*100
            return([i,total,per])
        i=i+1

def statsgen(path,Qth,DPth,name):
    inputvcfs = glob.glob(path+"/*.vcf")
    if path+"/all.merged.vcf" in inputvcfs:
        print(inputvcfs)
        inputvcfs.remove(path+"/all.merged.vcf")
    statsfn = path+"/"+name+"_QUAL_DP_stats.txt"
    reportfn = path+"/"+name+"_QUAL_DP_cutoffrate.txt"
    statsfile = open(statsfn,"w+")
    statsfile.write("id\tvalue\n")
    reportfile = open(reportfn,"w+")
    qual = []
    depth = []
    for f in inputvcfs:
        for l in open(f,"r"):
            l=l.strip()
            if not re.match("#",l):
                la=l.strip().split("\t")
                q = la[5]
                statsfile.write("QUAL\t"+q+"\n")
                qual.append(float(q))
                info = la[7].split(";")
                for x in info:
                        dp = x.split("=")
                        if dp[0] == "DP":
                            statsfile.write("DP\t"+dp[1]+"\n")
                            depth.append(float(dp[1]))
    sorted_q = sorted(qual)
    sorted_dp = sorted(depth)
    Qfiltp = get_rank(sorted_q,Qth)
    DPfiltp = get_rank(sorted_dp,DPth)
    reportfile.write("QUAL:"+"\t".join(str(x) for x in Qfiltp)+"\n"+"DP:"+"\t".join(str(x) for x in DPfiltp))
    return(statsfn)

#statsgen(originalvcfs)
flnklen=100
###merging vcf files
def getflank(cm,p):
    #flnklen = 100
    p=int(p)
    ps=1
    pe=100000000000000000000
    if p < flnklen:
        ps=1
        pe=p+flnklen
    else:
        ps=p-flnklen
        pe=p+flnklen

    cmd = "samtools faidx /netscratch/dep_psl/grp_rgo/taklee/genomes/arabidopsis/Araport11/Athaliana_447_TAIR10.fa "+cm+":"+str(ps)+"-"+str(pe)
    result = subprocess.check_output(cmd, shell=True)
    flseq = result.decode('utf-8').split("\n")
    seqs = ''
    for x in flseq:
        if ">" not in x:
            seqs = seqs+x
    return(list(seqs.upper()))

def merge_chkrep_vcfs(inpath,outpath):
    #flnklen
    merged=outpath+"/all.merged.vcf"
    mergeinputs=inpath+"/*.vcf.gz"
    inputs=inpath+"/*.vcf"
    inputfile=glob.glob(inputs)
    cmd = "bcftools merge --force-samples "+mergeinputs+" -o "+merged
    print(cmd)
    subprocess.run(cmd,shell=True)
    outn = outpath+"/allvcf_cnts_flankingseq_norepeats.txt"
    rptn = outpath+"/allvcf_cnts_flankingseq_repeats.txt"

    if os.path.exists(merged):
        cnts_norpt = {}
        cnts_rpt = {}
        outfile = open(outn,"w+")
        rptout = open(rptn,"w+")
        for l in open(merged,"r"):
            l=l.strip()
            if not re.match("#",l):
                la = l.strip().split("\t")
                cnt = 0
                for x in la[9:]:
                    if ".:." not in x:
                        cnt+=1
                flseq = getflank(la[0],la[1])
                ref = la[3].split(",")[0]
                var = la[4].split(",")[0]
                reflen = len(ref)
                flseq[flnklen] = "\""+la[3]+"/"+la[4]+"\""
                flsall = "".join(k for k in flseq[:flnklen+1])+"".join(j for j in flseq[flnklen+reflen:])
                chnge = la[0]+"\t"+la[1]+"\t"+flsall
                rptchk = ""
                if len(var) > reflen: 
                    rptchk = var.replace(ref,"",1)
                elif reflen > len(var):
                    rptchk = ref.replace(var,"",1)
                elif reflen == len(var):
                    rptchk = "dontchk"
                if rptchk == "dontchk":
                    cnts_norpt[chnge] = cnt
                else:
                    followseq = "".join(j for j in flseq[flnklen+reflen:])
                    repN = len(re.findall(rptchk,"".join(r for r in list(followseq)[:len(rptchk)*5])))
                    if len(rptchk) > 1:
                        if repN > 1:
                            cnts_rpt[chnge] = cnt
                        else:
                            cnts_norpt[chnge] = cnt
                    else:
                        if repN > 4:
                            cnts_rpt[chnge] = cnt
                        else:
                            cnts_norpt[chnge] = cnt
        sorted_cnts_rpt = sorted(cnts_rpt.items(),key=lambda x:x[1],reverse=True)
        sorted_cnts_norpt = sorted(cnts_norpt.items(),key=lambda x:x[1],reverse=True)
        for y in sorted_cnts_norpt:
            outfile.write(y[0]+"\t"+str(y[1])+"\n")
        for yy in sorted_cnts_rpt:
            rptout.write(yy[0]+"\t"+str(yy[1])+"\n")
        rptout.close()
        outfile.close()

#removing repetitive sequences
        repeats = {}
        for l in open(rptn,"r"):
            la=l.strip().split("\t")
            allele=la[0]+"\t"+la[1]
            repeats[allele] = "x"
        for f in inputfile:
            f=f.strip()
            norepname = outpath+"/"+f.split("/")[-1].replace(".vcf","_repeatfilt.vcf")
            norepfile = open(norepname,"w+")
            for l in open(f,"r"):
                l=l.strip()
                if re.match("#",l):
                    #norepfile.write(l.strip().replace("_pilon","")+"\n")
                    norepfile.write(l.strip()+"\n")
                else:
                    la = l.split("\t")
                    al = la[0]+"\t"+la[1]
                    if al not in repeats:
                        #ll = l.replace("_pilon","")
                        norepfile.write(l.strip()+"\n")

#remove ancestral variants
def filter_Ancestral(inputfiles,outpath,ancestral):
    for f in inputfiles:
        idxcmd = "tabix -p vcf "+f
        subprocess.run(idxcmd,shell=True)
        output=outpath+"/"+f.split("/")[-1].replace(".vcf.gz","_Ancestralfilt.vcf")
        renamesamp = output.replace(".vcf","_samprename.vcf")
        cmd = "vcf-isec -f -c "+f+" "+ancestral+" > "+output
        subprocess.run(cmd,shell=True)
        sna = f.split("/")[-1].split("_")
        samplename = "Ancestral"
        if sna[0] != "Ancestral":
            samplename = sna[0]+"_"+sna[1]+"_"+sna[2]
        with open(output,"r") as afilt, open(renamesamp,"w+") as renamed:
            for l in afilt:
                la=l.strip().split("\t")
                if la[0] == "#CHROM":
                    la[-1] = samplename
                    renamed.write("\t".join(la)+"\n")
                else:
                    renamed.write(l.strip()+"\n")
        cleanup1 = f"rm {output}"
        cleanup2 = f"mv {renamesamp} {output}"
        subprocess.run(cleanup1,shell=True)
        subprocess.run(cleanup2,shell=True)


