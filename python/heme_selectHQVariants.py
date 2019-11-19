import sys
import pandas as pd
##################################
### input variables
##################################
vc1_file=sys.argv[1]
vc2_file=sys.argv[2]
vc3_file=sys.argv[3]
out_file=sys.argv[4]
##################################S
### read files
##################################

df_vc1 = pd.read_csv(vc1_file,sep='\t',skiprows=1)
main_columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTE', 'INFO', 'FORMAT','sample']
df_vc1.columns=main_columns
df_vc2 = pd.read_csv(vc2_file,sep='\t')
df_vc3 = pd.read_csv(vc3_file,sep='\t')


##################################
### vc1 + vc2(mutect)
##################################
df_12 = pd.merge(df_vc1,df_vc2,on=['CHROM','POS','REF','ALT'],how='left',indicator=True)
df_12 = df_12[df_12['_merge']=="both"]
df_12 = df_12.iloc[:,0:10]
df_12.columns=main_columns


##################################
### vc1 + vc3 (freebayes)
##################################
df_13 = pd.merge(df_vc1,df_vc3,on=['CHROM','POS'],how='left',indicator=True)
df_13 = df_13[df_13['_merge']=="both"]
df_13 = df_13.iloc[:,0:10]
df_13.columns=main_columns


##################################
### vc2(mutect) + vc3 (freebayes)
##################################
df_23 = pd.merge(df_vc2,df_vc3,on=['CHROM','POS'],how='left',indicator=True)
df_23 = df_23[df_23['_merge']=="both"]
df_23 = df_23.iloc[:,0:10]

for index, row in df_23.iterrows():
    GT=row[7]
    REF_DP=int(row[8].split(',')[0])
    ALT_DP=int(row[8].split(',')[1])
    DP= REF_DP + ALT_DP
    VAF = round(row[9],4)

    #FIX format same as varscan output
    df_23.ix[index,7] = "DP=%s" % (DP)

    df_23.ix[index,8] = "GT:VF:DP:AD"

    df_23.ix[index,9] = "%s:%s:%s:%s,%s" % (GT,VAF,DP,REF_DP,ALT_DP)

df_23.columns=main_columns
df_23['FILTE']="PASS"

##################################
### merge all
##################################

df_all = pd.concat([df_12,df_13,df_23])
df_all.drop_duplicates(subset=['CHROM','POS','REF','ALT'],keep='first',inplace=True)
# df_all.reset_index(inplace=True)
df_all.to_csv(out_file,sep='\t',header=False,index=False)
