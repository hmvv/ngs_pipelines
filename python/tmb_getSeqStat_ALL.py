import sys
import pandas as pd
import numpy as np


DIR=sys.argv[1] 
#example--"/home/environments/ngs_test/exomeAnalysis/Single/190318_NS500761_0325_AH52LGBGXB/HD799_S1/Alignment/"
SAMPLE=sys.argv[2]
#example--"HD799_S1"

def get_triminfo (f1):
    f = open(f1,"r")
    for line in f:
        wl = line.split(' ')
        if wl[0]=="Input":
            print("Total")
            print (int(wl[3])*2)
            print(" ")
            print("Q20-ATW(10)M(30)")
            print (int(wl[6])*2)
            print("\n")


def get_aligninfo (f1):
    df = pd.read_csv(f1,skiprows=6,sep='\t')
    print("Alignment-mapQ20")
    print(df[df["CATEGORY"]=="PAIR"]["TOTAL_READS"].values[0])
    print("\n")



def get_dupinfo (f1):
    df = pd.read_csv(f1,skiprows=6,sep='\t',nrows=1)
    print("Duplicates")
    print(df["PERCENT_DUPLICATION"].values[0])
    print("\n")

def get_seqinfo (f1):
    df = pd.read_csv(f1,skiprows=6,sep='\t',nrows=1)
    print(df["TOTAL_READS"].values[0])
    print("\n")
    print("Coverage")
    print(str(int(df["MEAN_TARGET_COVERAGE"].values[0]))+"X")
    print("\n")
    print(str(int(df["PCT_SELECTED_BASES"].values[0] * 100))+"%")
    print(str(int(df["FOLD_ENRICHMENT"].values[0]))+"X")
    print(str(df["ZERO_CVG_TARGETS_PCT"].values[0]* 100)+"%")
    print(df["FOLD_80_BASE_PENALTY"].values[0])

    print(str(df["PCT_TARGET_BASES_2X"].values[0]* 100)+"%")
    print(str(df["PCT_TARGET_BASES_10X"].values[0]* 100)+"%")
    print(str(df["PCT_TARGET_BASES_20X"].values[0]* 100)+"%")
    print(str(df["PCT_TARGET_BASES_30X"].values[0]* 100)+"%")
    print(str(df["PCT_TARGET_BASES_40X"].values[0]* 100)+"%")
    print(str(df["PCT_TARGET_BASES_50X"].values[0]* 100)+"%")
    print(str(df["PCT_TARGET_BASES_100X"].values[0]* 100)+"%")


print(SAMPLE+"\n")
file_trim=DIR+SAMPLE+".trimmomatic.summary.txt"
get_triminfo(file_trim)

file_align=DIR+SAMPLE+".sorted.bam.alignmentMetrics.txt"
get_aligninfo(file_align)

file_dups=DIR+SAMPLE+".sorted.rmdups.bam.metrics.txt"
get_dupinfo(file_dups)

file_seq=DIR+SAMPLE+".output_hs_metrics.txt"
get_seqinfo(file_seq)
