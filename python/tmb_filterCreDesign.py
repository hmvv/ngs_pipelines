import sys
import pandas as pd


##################################
### input variables
##################################

OUT_DIR=sys.argv[1]
SAMPLE=sys.argv[2]
ENV=sys.argv[3]

refstd=OUT_DIR+SAMPLE
df_ref = pd.read_csv(refstd,sep='\t')

DESIGN_FILE= "/home/environments/ngs_"+ENV+"/assayCommonFiles/tmbAssay/cre_design.bed"
df_design = pd.read_csv(DESIGN_FILE,sep='\t')
df_design.columns=['CHROM','START','END']


filter=[]
for indx,row in df_ref.iterrows():
    if ( df_design[ (df_design['CHROM'] == row['CHROM']) & (df_design['START'] <= row['POS']) & (df_design['END'] >= row['POS']) ].shape[0]>=1):
        filter.append([row['CHROM'],row['POS'],row['ID'],row['REF'],row['ALT'],row['QUAL'],row['FILTER'],row['INFO']])

df_filter=pd.DataFrame(filter)
df_filter.columns=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']
df_filter.to_csv(OUT_DIR+SAMPLE+".crefilter",sep='\t',index=False)
