
import pandas as pd
import optparse

def filterVariants(inputfile, outputfile, assay):

    if assay == 'neuro':

        exclude = ["MUC17","MUC16"]

        df = pd.read_csv(inputfile,sep='\t',header=None)
        df = df[~df[0].isin(exclude)]
        df.to_csv(outputfile, sep='\t', index=False, header=False)

    elif assay == 'gene50':

        df = pd.read_csv(inputfile,sep='\t',header=None)
        df.to_csv(outputfile, sep='\t', index=False, header=False)

try:
    parser = optparse.OptionParser()
    parser.add_option('-i', '--infile', help = 'provide input file')
    parser.add_option('-o', '--outfile', help = 'provide output file')
    parser.add_option('-a', '--assay', help = 'provide assay')
    options,args = parser.parse_args()
    infile = options.infile
    outfile = options.outfile
    assay = options.assay
    filterVariants(infile, outfile, assay)

except TypeError:
	print ("python filterVEP.py -help for help")
