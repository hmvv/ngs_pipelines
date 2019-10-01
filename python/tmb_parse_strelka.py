import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import warnings
from scipy.stats import norm,expon
from sklearn.mixture import GaussianMixture
from pomegranate import GeneralMixtureModel,ExponentialDistribution,NormalDistribution
warnings.filterwarnings("ignore")


def getGaussianMixCutOff(samples):
    mixture = GaussianMixture(n_components=2).fit(samples.reshape(-1, 1))
    if mixture.converged_:
        means_hat = mixture.means_.flatten()
        sds_hat = np.sqrt(mixture.covariances_).flatten()
        fp_index=means_hat.argmin()
        mu= means_hat[fp_index]
        sigma=sds_hat[fp_index]
        return norm.ppf(0.95,mu,sigma)
    else :
        return 0.0

def getMixtureModelCutOff(samples,alpha,mu,sigma):
    mixture_m = GeneralMixtureModel([ ExponentialDistribution(alpha), NormalDistribution(mu,sigma)] )
    model = mixture_m.fit(samples.reshape(-1,1))
    pred_alpha = model.distributions[0].parameters[0]
    return expon.ppf(0.95,0,1/pred_alpha)

def getSNPAlleleFreq(ref, a_count, t_count,g_count, c_count,tier):
    allele_count=0
    if ref=='A':
        allele_count=int(a_count.split(',')[tier])
    elif ref=='T':
        allele_count=int(t_count.split(',')[tier])
    elif ref=='G':
        allele_count=int(g_count.split(',')[tier])
    elif ref=='C':
        allele_count=int(c_count.split(',')[tier])

    total_allele= sum([int(x.split(',')[tier]) for x in [a_count, t_count,g_count, c_count]])

    if total_allele == 0.0:
        return 0.0
    else:
        return (allele_count/total_allele)

def getINDELAlleleFreq(type, normal_count, indel_count,tier):
    allele_count =0.0
    if type=="normal" :
        allele_count=int(normal_count.split(',')[tier])
    elif type=="tumor":
        allele_count=int(indel_count.split(',')[tier])

    total_allele= int(normal_count.split(',')[tier]) + int(indel_count.split(',')[tier])

    if total_allele == 0.0:
        return 0.0
    else:
        return (allele_count/total_allele)


##################################
### read files
##################################
SAMPLE=sys.argv[1]
OUT_DIR=sys.argv[2]
ENV=sys.argv[3]
DEPTH=sys.argv[4]
NALF=sys.argv[5]
TALF=sys.argv[6]
#########################

QSS_ALPHA = 1/50
QSS_MU = 600
QSS_SIGMA = 1000

NORMAL_DP_ALPHA = 1/70
NORMAL_DP_MU = 250
NORMAL_DP_SIGMA = 1000

TUMOR_DP_ALPHA = 1/50
TUMOR_DP_MU = 200
TUMOR_DP_SIGMA = 1000

NORMAL_REF_FREQ = 1.0-(float(NALF)/100.0)
NORMAL_ALT_FREQ = (float(NALF)/100.0)
TUMOR_ALT_FREQ =(float(TALF)/100.0)

MIN_DEPTH = int(DEPTH)


###### process snp
file_snp=OUT_DIR+"output.snvs.vcf.filter.txt"
df_snp = pd.read_csv(file_snp,sep='\t')
df_snp['GROUP']='snp'

###filtering:
filter_indx = df_snp[df_snp['CHROM'].str.contains("GL|gl")].index
df_snp.drop(filter_indx, inplace=True)
filter_indx = df_snp[df_snp['CHROM'].str.contains("chrM")].index
df_snp.drop(filter_indx, inplace=True)
df_snp = df_snp[ df_snp['FILTER']=="PASS"]


df_snp['NORMAL_REF_FREQ']= df_snp.apply(lambda row: getSNPAlleleFreq(row['REF'],row['NORMAL.AU'],row['NORMAL.TU'],row['NORMAL.GU'],row['NORMAL.CU'],0),axis=1)
df_snp['NORMAL_ALT_FREQ']= df_snp.apply(lambda row: getSNPAlleleFreq(row['ALT'],row['NORMAL.AU'],row['NORMAL.TU'],row['NORMAL.GU'],row['NORMAL.CU'],0),axis=1)
df_snp['TUMOR_ALT_FREQ']= df_snp.apply(lambda row: getSNPAlleleFreq(row['ALT'],row['TUMOR.AU'],row['TUMOR.TU'],row['TUMOR.GU'],row['TUMOR.CU'],0),axis=1)



df_snp= df_snp[['CHROM', 'POS', 'REF', 'ALT','SomaticEVS','NORMAL.DP','TUMOR.DP','NORMAL_REF_FREQ',\
                'NORMAL_ALT_FREQ','TUMOR_ALT_FREQ','QSS','QSS_NT','GROUP','SNVSB',\
                'NORMAL.FDP','NORMAL.SDP', 'NORMAL.SUBDP','TUMOR.FDP', 'TUMOR.SDP', 'TUMOR.SUBDP','MQ']]


### filtering parameters:
SomaticEVS_Cutoff_snp = getGaussianMixCutOff(df_snp['SomaticEVS'])
QSS_Cutoff_snp=getMixtureModelCutOff(df_snp['QSS'],QSS_ALPHA,QSS_MU,QSS_SIGMA)
Normal_DP_snp=getMixtureModelCutOff(df_snp['NORMAL.DP'],NORMAL_DP_ALPHA,NORMAL_DP_MU,NORMAL_DP_SIGMA)
Tumor_DP_snp=getMixtureModelCutOff(df_snp['TUMOR.DP'],TUMOR_DP_ALPHA,TUMOR_DP_MU,TUMOR_DP_SIGMA)



###filtering:
df_snp = df_snp[ df_snp['SomaticEVS']>=SomaticEVS_Cutoff_snp]
df_snp = df_snp[ df_snp['QSS']>=QSS_Cutoff_snp]
df_snp = df_snp[ df_snp['NORMAL.DP']>=Normal_DP_snp]
df_snp = df_snp[ df_snp['TUMOR.DP']>=Tumor_DP_snp]
df_snp = df_snp[ df_snp['NORMAL_REF_FREQ']>=NORMAL_REF_FREQ]
df_snp = df_snp[ df_snp['NORMAL_ALT_FREQ']<=NORMAL_ALT_FREQ]
df_snp = df_snp[ df_snp['TUMOR_ALT_FREQ']>=TUMOR_ALT_FREQ]



df_snp= df_snp[['CHROM', 'POS', 'REF', 'ALT','NORMAL.DP','TUMOR.DP','NORMAL_REF_FREQ',\
                'NORMAL_ALT_FREQ','TUMOR_ALT_FREQ','SomaticEVS','QSS','QSS_NT']]


############ processing indels
file_indel=OUT_DIR+"output.indels.vcf.filter.txt"
df_indel = pd.read_csv(file_indel,sep='\t')
df_indel['GROUP']='indel'


###filtering:
filter_indx = df_indel[df_indel['CHROM'].str.contains("GL|gl")].index
df_indel.drop(filter_indx, inplace=True)
filter_indx = df_indel[df_indel['CHROM'].str.contains("chrM")].index
df_indel.drop(filter_indx, inplace=True)
df_indel = df_indel[ df_indel['FILTER']=="PASS"]


df_indel['NORMAL_REF_FREQ']= df_indel.apply(lambda row : getINDELAlleleFreq("normal", row['NORMAL.TAR'], row['NORMAL.TIR'],0),axis=1)
df_indel['NORMAL_ALT_FREQ']= df_indel.apply(lambda row :  getINDELAlleleFreq("tumor", row['NORMAL.TAR'], row['NORMAL.TIR'],0),axis=1)
df_indel['TUMOR_ALT_FREQ']=  df_indel.apply(lambda row : getINDELAlleleFreq("tumor", row['TUMOR.TAR'], row['TUMOR.TIR'],0),axis=1)


df_indel = df_indel[['CHROM', 'POS', 'REF', 'ALT','SomaticEVS','NORMAL.DP','TUMOR.DP','NORMAL_REF_FREQ','NORMAL_ALT_FREQ','TUMOR_ALT_FREQ','QSI','QSI_NT','GROUP']]


### filtering parameters:

SomaticEVS_Cutoff_indel = getGaussianMixCutOff(df_indel['SomaticEVS'])

# ###filtering:
SomaticEVS_Cutoff_indel = getGaussianMixCutOff(df_indel['SomaticEVS'])
QSS_Cutoff_indel=getMixtureModelCutOff(df_indel['QSI'],QSS_ALPHA,QSS_MU,QSS_SIGMA)
Normal_DP_indel=getMixtureModelCutOff(df_indel['NORMAL.DP'],NORMAL_DP_ALPHA,NORMAL_DP_MU,NORMAL_DP_SIGMA)
Tumor_DP_indel=getMixtureModelCutOff(df_indel['TUMOR.DP'],TUMOR_DP_ALPHA,TUMOR_DP_MU,TUMOR_DP_SIGMA)


df_indel = df_indel[ df_indel['SomaticEVS']>=SomaticEVS_Cutoff_indel]
df_indel = df_indel[ df_indel['QSI']>=QSS_Cutoff_indel]
df_indel = df_indel[ df_indel['NORMAL.DP']>=Normal_DP_indel]
df_indel = df_indel[ df_indel['TUMOR.DP']>=Tumor_DP_indel]
df_indel = df_indel[ df_indel['NORMAL_REF_FREQ']>=NORMAL_REF_FREQ]
df_indel = df_indel[ df_indel['NORMAL_ALT_FREQ']<=NORMAL_ALT_FREQ]
df_indel = df_indel[ df_indel['TUMOR_ALT_FREQ']>=TUMOR_ALT_FREQ]


df_indel = df_indel[['CHROM', 'POS', 'REF', 'ALT','NORMAL.DP','TUMOR.DP','NORMAL_REF_FREQ','NORMAL_ALT_FREQ','TUMOR_ALT_FREQ','SomaticEVS','QSI','QSI_NT']]
df_indel.columns= ['CHROM', 'POS', 'REF', 'ALT','NORMAL.DP','TUMOR.DP','NORMAL_REF_FREQ','NORMAL_ALT_FREQ','TUMOR_ALT_FREQ','SomaticEVS','QSS','QSS_NT']


# ##################################
# ### combine snp + indel
# ##################################

df_both = pd.concat([df_snp,df_indel])
# df_both.to_csv(OUT_DIR+"strelka.filter.txt",sep='\t',index=False)

##################################
### format for VEP
##################################
vcf=[]
for indx,row in df_both.iterrows():
    if row['NORMAL.DP'] <= MIN_DEPTH or row['TUMOR.DP'] <= MIN_DEPTH:
        continue
    info=str(row['NORMAL.DP'])+';'+str(row['TUMOR.DP'])
    vcf.append([row['CHROM'],row['POS'],'.',row['REF'],row['ALT'],'.','.',info])

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
    df_filter.to_csv(OUT_DIR+SAMPLE+".strelka."+DEPTH+"_"+NALF+"_"+TALF,sep='\t',index=False)

else:
    df_vcf.to_csv(OUT_DIR+SAMPLE+".strelka."+DEPTH+"_"+NALF+"_"+TALF,sep='\t',index=False)
