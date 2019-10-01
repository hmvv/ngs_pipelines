import sys
import matplotlib
import pandas as pd
import numpy as np

def addVCInfo(row,df):
    if(df[ ( \
            (df['CHROM']==row['CHROM']) & \
            (df['POS']==row['POS']) & \
            (df['REF']==row['REF']) & \
            (df['ALT']==row['ALT']) )].shape[0] == 1):
        return 1
    else:
        return 0

##################################
### input variables
##################################
SAMPLE=sys.argv[1]
OUT_DIR=sys.argv[2]
ENV=sys.argv[3]
type=sys.argv[4]
############################################################


file1 = OUT_DIR+"/varscan/"+SAMPLE+".varscan."+type
df_f1 = pd.read_csv(file1,sep='\t')
file2 = OUT_DIR+"/mutect/"+SAMPLE+".mutect."+type
df_f2 = pd.read_csv(file2,sep='\t')
file3 = OUT_DIR+"/strelka/"+SAMPLE+".strelka."+type
df_f3 = pd.read_csv(file3,sep='\t')



df_all = pd.concat([df_f1,df_f2,df_f3])
df_all.drop_duplicates(subset=['CHROM','POS','ID','REF','ALT'],keep='first',inplace=True)
df_all.reset_index(inplace=True)
print("Unique variants-",df_all.shape)

for indx,row in df_all.iterrows():
    temp=[]
    temp.append(addVCInfo(row,df_f1))
    temp.append(addVCInfo(row,df_f2))
    temp.append(addVCInfo(row,df_f3))
    if sum(temp) >=2:
        df_all.at[indx,'INFO'] += ';'+str(temp[0]) + ';'+str(temp[1]) + ';'+str(temp[2])
    else:
        df_all.drop(indx,inplace=True)

df_all['POS'] = df_all['POS'].astype(int)
df_all=df_all[['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']]
df_all.to_csv(OUT_DIR+"/"+SAMPLE+".variantcallers.combine."+type,index=False,header=True,sep='\t')
