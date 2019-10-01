import sys
import pandas as pd

SAMPLEID=sys.argv[1]
SAMPLE=sys.argv[2]
OUT_DIR=sys.argv[3]
ENV=sys.argv[4]
type=sys.argv[5]

print(SAMPLEID)
print(SAMPLE)
print(OUT_DIR)
print(ENV)
print(type)

VEP_FILE=OUT_DIR+"/"+SAMPLE+".variantcallers.combinev2."+type+".vep"

header=[]
with open(VEP_FILE) as myfile:
    header = [next(myfile) for x in range(3)]
columns_line=str(header[2].strip().split(':')[1]).split('|')

df_vep = pd.read_csv(VEP_FILE,sep='\t',skiprows=4,header=None)
df_vep.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
df_vep[columns_line] = df_vep['INFO'].str.split('|',expand=True)
df_vep.drop(['INFO'],axis=1,inplace=True)
df_vep= df_vep[df_vep['GIVEN_REF']==df_vep['USED_REF']]
costum_line=['ND','TD','VARSCAN','MUTECT','STRELKA','CSQ']
df_vep[costum_line] = df_vep[' Allele'].str.split(';',expand=True)
df_vep.drop([' Allele'],axis=1,inplace=True)


### filter VOI i.e. variants of interest
df_vep_2 = df_vep[df_vep['Consequence'].str.contains("missense_variant|inframe_insertion|inframe_deletion|stop_gained|frameshift_variant|coding_sequence_variant")]
df_vep_2.to_csv(OUT_DIR+"/"+SAMPLE+".variantcallers.combinev2."+type+".vep.parse.txt",sep='\t',index=False)


### TMB Calculation
tmb_score=round(float(df_vep_2.shape[0])/35.50,3)

tmb_group=""
if tmb_score<2.0:
    tmb_group="LOW"
elif tmb_score>2.0 and tmb_score<3.0:
    tmb_group="INTERMEDIATE"
else:
    tmb_group="HIGH"

print(SAMPLEID+","+SAMPLE+","+str(df_vep_2.shape[0])+","+str(tmb_score)+","+tmb_group)
with open(OUT_DIR+"/"+SAMPLE+".tmb_result", 'w') as f:
    f.write(SAMPLEID+","+SAMPLE+","+str(df_vep_2.shape[0])+","+str(tmb_score)+","+tmb_group);f.close()
