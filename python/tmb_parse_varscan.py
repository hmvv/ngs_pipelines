import sys
import pandas as pd
##################################
### input variables
##################################
SAMPLE=sys.argv[1]
OUT_DIR=sys.argv[2]
ENV=sys.argv[3]
DEPTH=sys.argv[4]
NALF=sys.argv[5]
TALF=sys.argv[6]

##################################
### filter parameters
##################################
VARSCAN_PVAL=0.01
NORMAL_DEPTH=int(DEPTH)
TUMOR_DEPTH=int(DEPTH)
NORMAL_ALT_ALLELE_FREQ=float(NALF)
TUMOR_ALT_ALLELE_FREQ=float(TALF)
NORMAL_STRAND_BIAS_PLUS=0
NORMAL_STRAND_BIAS_MINUS=0
TUMOR_STRAND_BIAS_1_MINUS=0
TUMOR_STRAND_BIAS_1_PLUS=0
TUMOR_STRAND_BIAS_2_MINUS=0
TUMOR_STRAND_BIAS_2_PLUS=0
FDR=0.01

##################################S
### read files
##################################
file_snp=OUT_DIR+SAMPLE+".varscan.output.snp"
file_indel=OUT_DIR+SAMPLE+".varscan.output.indel"
df_snp = pd.read_csv(file_snp,sep='\t')
df_snp['chrom'] =df_snp['chrom'].apply(str)
df_snp['group']='snp'

df_indel = pd.read_csv(file_indel,sep='\t')
df_indel['chrom'] =df_indel['chrom'].apply(str)
df_indel['group']='indel'
for indx, row in df_indel.iterrows():
    if row['var'][0] == '-':
        df_indel.set_value(indx, "ref", row['ref']+row['var'][1:])
        df_indel.set_value(indx, "var", row['ref'])
    elif row['var'][0] == '+':
        df_indel.set_value(indx, "var", row['ref'] + row['var'][1:] )

##################################
### combine snp + indel
##################################
df_both = pd.concat([df_snp,df_indel])
df_both.sort_values(by=['somatic_status'],inplace=True)

##################################
### remove unwanted chr
##################################
filter_indx = df_both[df_both['chrom'].str.contains("GL|gl")].index
df_both.drop(filter_indx, inplace=True)
filter_indx = df_both[df_both['chrom'].str.contains("chrM")].index
df_both.drop(filter_indx, inplace=True)

##################################
### filter somatic
##################################
df_both_pass_s0 = df_both[(df_both['somatic_status'].isin(['Somatic']) )]


##################################
### update additional columns
##################################
df_both_pass_s0.loc[:,'normal_var_freq_percent'] = df_both_pass_s0['normal_var_freq'].str.replace('%','' ).astype(float)
df_both_pass_s0.loc[:,'tumor_var_freq_percent'] = df_both_pass_s0['tumor_var_freq'].str.replace('%','' ).astype(float)
df_both_pass_s0.loc[:,'normal_depth'] = df_both_pass_s0['normal_reads1'] + df_both_pass_s0['normal_reads2']
df_both_pass_s0.loc[:,'tumor_depth'] = df_both_pass_s0['tumor_reads1'] + df_both_pass_s0['tumor_reads2']


##################################
### raw before filter
##################################
# df_both_pass_s0.to_csv(OUT_DIR+"/"+SAMPLE+"_raw_before_filter."+DEPTH+"_"+NALF+"_"+TALF+".txt",sep='\t',index=False)


##################################
### pval + FDR
##################################
def pvalAdj(p_val_rank,total,fdr):
    return (p_val_rank*fdr)/float(total)

## line 93 throws  SettingWithCopyWarning # WARNING: -- this is fine
df_both_pass_s0.sort_values(by=['somatic_p_value'],inplace=True)
df_both_pass_s0.reset_index(inplace=True)

# df_both_pass_s0 = df_both_pass_s0[df_both_pass_s0['somatic_p_value']<=VARSCAN_PVAL]

df_both_pass_s0.loc[:,'somatic_p_value_rank']= [ x for x in range(1,df_both_pass_s0.shape[0]+1)]
df_both_pass_s0.loc[:,'somatic_p_value_adj']= [ pvalAdj(x,df_both_pass_s0.shape[0],FDR) for x in df_both_pass_s0['somatic_p_value_rank']]
df_both_pass_s0.loc[:,'somatic_p_value_adj_filter']=[1 if x[0]<x[1] else 0 for x in zip(df_both_pass_s0['somatic_p_value'], df_both_pass_s0['somatic_p_value_adj']) ]
df_both_pass_s0 = df_both_pass_s0[df_both_pass_s0['somatic_p_value_adj_filter']==1]


#################################
## alt frequency
#################################
df_both_pass_s0 = df_both_pass_s0[df_both_pass_s0['normal_var_freq_percent']<NORMAL_ALT_ALLELE_FREQ]
df_both_pass_s0 = df_both_pass_s0[df_both_pass_s0['tumor_var_freq_percent']>=TUMOR_ALT_ALLELE_FREQ]

#################################
## depth
#################################
df_both_pass_s0 = df_both_pass_s0[ (df_both_pass_s0['normal_depth']>=NORMAL_DEPTH) & (df_both_pass_s0['tumor_depth']>=TUMOR_DEPTH) ]

##################################
### strand bias
##################################
df_both_pass_s0 = df_both_pass_s0[ ( (df_both_pass_s0['normal_reads1_plus']>=NORMAL_STRAND_BIAS_PLUS) &
                                     (df_both_pass_s0['normal_reads1_minus']>=NORMAL_STRAND_BIAS_MINUS) &
                                     (df_both_pass_s0['tumor_reads1_plus']>=TUMOR_STRAND_BIAS_1_PLUS) &
                                     (df_both_pass_s0['tumor_reads1_minus']>=TUMOR_STRAND_BIAS_1_MINUS)  &
                                      (df_both_pass_s0['tumor_reads2_plus']>=TUMOR_STRAND_BIAS_2_PLUS) &
                                      (df_both_pass_s0['tumor_reads2_minus']>=TUMOR_STRAND_BIAS_2_MINUS)  ) ]

# print('After strand bias:'+ str(df_both_pass_s0.shape[0]))
# df_both_pass_s0.to_csv(OUT_DIR+"/"+SAMPLE+"_filter_final."+DEPTH+"_"+NALF+"_"+TALF+".txt",sep='\t',index=False)
# print('Filter:'+ str(df_both_pass_s0.shape[0]))
# print('Filter:snp:'+ str(df_both_pass_s0[df_both_pass_s0['group']=='snp'].shape[0]))
# print('Filter:indel:'+ str(df_both_pass_s0[df_both_pass_s0['group']=='indel'].shape[0]))

##################################
### format for VEP
##################################
vcf=[]
for indx,row in df_both_pass_s0.iterrows():
    info=str(row['normal_depth'])+';'+str(row['tumor_depth'])
    vcf.append([row['chrom'],row['position'],'.',row['ref'],row['var'],'.','.',info])

df_vcf=pd.DataFrame(vcf)
df_vcf.columns=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']

##################################
### filter CRE Design region
##################################

if "COLO" in SAMPLE:
    DESIGN_FILE= "/home/environments/ngs_"+ENV+"/assayCommonFiles/tmbAssay/cre_design.bed"
    df_design = pd.read_csv(DESIGN_FILE,sep='\t')
    df_design.columns=['CHROM','START','END']


    filter=[]
    for indx,row in df_vcf.iterrows():
        if ( df_design[ (df_design['CHROM'] == row['CHROM']) & (df_design['START'] <= row['POS']) & (df_design['END'] >= row['POS']) ].shape[0]>=1):
            filter.append([row['CHROM'],row['POS'],row['ID'],row['REF'],row['ALT'],row['QUAL'],row['FILTER'],row['INFO']])

    df_filter=pd.DataFrame(filter)
    df_filter.columns=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']
    df_filter.to_csv(OUT_DIR+SAMPLE+".varscan."+DEPTH+"_"+NALF+"_"+TALF,sep='\t',index=False)

else:
    df_vcf.to_csv(OUT_DIR+SAMPLE+".varscan."+DEPTH+"_"+NALF+"_"+TALF,sep='\t',index=False)
