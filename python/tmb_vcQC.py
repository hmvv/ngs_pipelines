import sys
import pandas as pd

def calcTiTvRatio(df):

    AtoG = df[ (df.REF=='A') & (df.ALT=='G') ].shape[0]
    GtoA = df[ (df.REF=='G') & (df.ALT=='A') ].shape[0]
    CtoT = df[ (df.REF=='C') & (df.ALT=='T') ].shape[0] ## very common
    TtoC = df[ (df.REF=='T') & (df.ALT=='C') ].shape[0]

    transition = AtoG + GtoA + CtoT + TtoC

    AtoC = df[ (df.REF=='A') & (df.ALT=='C') ].shape[0]
    CtoA = df[ (df.REF=='C') & (df.ALT=='A') ].shape[0]

    AtoT = df[ (df.REF=='A') & (df.ALT=='T') ].shape[0]
    TtoA = df[ (df.REF=='T') & (df.ALT=='A') ].shape[0]

    GtoC = df[ (df.REF=='G') & (df.ALT=='C') ].shape[0]
    CtoG = df[ (df.REF=='C') & (df.ALT=='G') ].shape[0]

    GtoT = df[ (df.REF=='G') & (df.ALT=='T') ].shape[0]
    TtoG = df[ (df.REF=='T') & (df.ALT=='G') ].shape[0]

    transversions = AtoC + CtoA + AtoT + TtoA + GtoC + CtoG + GtoT + TtoG

    return float(transition/transversions)


def calcIndelRatio(df):
    insertion = df[(df.REF.str.len()==1)& (df.ALT.str.len()>1)].shape[0]
    deletion = df[(df.REF.str.len()>1)& (df.ALT.str.len()==1)].shape[0]

    return float(insertion/deletion)


df1=pd.read_csv(sys.argv[1],sep='\t')
with open(sys.argv[2], 'w') as f: f.write("titv ratio: "+str(calcTiTvRatio(df1)));f.close()
