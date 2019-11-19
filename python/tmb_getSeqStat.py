import sys
import pandas as pd
import numpy as np

def get_triminfo (f1,res):
    f = open(f1,"r")
    for line in f:
        wl = line.split(' ')
        if wl[0]=="Input":
            res.append(int(wl[3])*2)
            res.append(str(int((((int(wl[6])*2) / (int(wl[3])*2)) * 100)))+"%")


def get_aligninfo (f1,res):
    df = pd.read_csv(f1,skiprows=6,sep='\t')
    res.append(df[df["CATEGORY"]=="PAIR"]["TOTAL_READS"].values[0])

def get_dupinfo (f1,res):
    df = pd.read_csv(f1,skiprows=6,sep='\t',nrows=1)
    res.append(str(int(df["PERCENT_DUPLICATION"].values[0]*100))+"%")

def get_seqinfo (f1,res):
    df = pd.read_csv(f1,skiprows=6,sep='\t',nrows=1)
    res.append(str(int(df["TOTAL_READS"].values[0])))
    res.append(str(int(df["MEAN_TARGET_COVERAGE"].values[0]))+"x")
    res.append(str(int(df["PCT_SELECTED_BASES"].values[0] * 100))+"%")
    res.append(str(int(df["PCT_TARGET_BASES_10X"].values[0]* 100))+"%")
    res.append(str(int(df["PCT_TARGET_BASES_20X"].values[0]* 100))+"%")
    res.append(str(int(df["PCT_TARGET_BASES_50X"].values[0]* 100))+"%")
    res.append(str(int(df["PCT_TARGET_BASES_100X"].values[0]* 100))+"%")

def getstat(DIR,SAMPLE):
    result=[]
    file_trim=DIR+SAMPLE+".trimmomatic.summary.txt"
    get_triminfo(file_trim,result)

    file_align=DIR+SAMPLE+".sort.bam.alignmentMetrics.txt"
    get_aligninfo(file_align,result)

    file_dups=DIR+SAMPLE+".sort.rmdups.bam.metrics.txt"
    get_dupinfo(file_dups,result)

    file_seq=DIR+SAMPLE+".output_hs_metrics.txt"
    get_seqinfo(file_seq,result)

    return result


DIR=sys.argv[1]
Tumor=sys.argv[2]
DIR2=sys.argv[3]
Normal=sys.argv[4]
OUT_DIR=sys.argv[5]

df_stats = pd.DataFrame({'#Metrics':['Total-Reads', 'Q20', 'Total-Reads-AQC', 'Duplicate', 'Total-Reads-ADup','Coverage', 'Target-Coverage', 'Coverage-10X','Coverage-20X','Coverage-50X','Coverage-100X']})
df_stats['Tumor']=getstat(DIR,Tumor)
df_stats['Normal']=getstat(DIR2,Normal)

df_stats.to_csv(OUT_DIR+Tumor+"_"+Normal+".seq_stats",index=False)
