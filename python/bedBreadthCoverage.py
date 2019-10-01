import sys
import pandas as pd


def updateSize(tab,chr,start,end,filter):
    if start >= tab[chr][1]:
        tab[chr][2] += (end - start)
        tab[chr][0] = start
        tab[chr][1] = end
        filter.append([chr,start,end,end-start])
    elif start <= tab[chr][1] and end > tab[chr][1]:
        tab[chr][0] = tab[chr][1]+1
        tab[chr][1] = end
        filter.append( [ chr,tab[chr][0],tab[chr][1],( tab[chr][1] - tab[chr][0])])
        tab[chr][2] += (tab[chr][1]-tab[chr][0])


def calSum(tab):
    sum=0
    for chr in tab:
        sum += tab[chr][2]
    return sum


file=sys.argv[1]
df = pd.read_csv(file,sep='\t',header=None)
df.columns = ['chrom','start','end']

tab={}
filter=[]
for indx,row in df.iterrows():
    if (row['chrom']) in tab:
        updateSize(tab,row['chrom'],row['start'],row['end'],filter)
    else:
        tab[row['chrom']]=[row['start'],row['end'],row['end']-row['start']]
        filter.append([row['chrom'],row['start'],row['end'],row['end']-row['start']])


print(calSum(tab))
pd.DataFrame(filter).to_csv(file+'_filter.csv',index=False,header=None)
