
import pandas as pd
import optparse

def loadIntoDB(inputfile, outputfile, sample, instrument):

    if instrument == 'proton':
        df_amplicon = pd.read_csv(inputfile, sep="\t", header=0)
        df_amplicon = df_amplicon[['contig_id', 'contig_srt', 'contig_end','region_id','attributes','total_reads']]
        df_amplicon['ampliconName'] = (df_amplicon['region_id'].astype(str) + '.' +
                        df_amplicon['contig_id'].astype(str) + '.' +
                        df_amplicon['contig_srt'].astype(str) + '.' +
                        df_amplicon['contig_end'].astype(str) + '.' +
                        df_amplicon['attributes'].astype(str) )
        df_amplicon = df_amplicon[['ampliconName','total_reads','region_id']]
        df_amplicon.columns = ['ampliconName','readDepth','gene']
        df_amplicon['sampleID'] = sample
        df_amplicon = df_amplicon[['sampleID','gene','ampliconName','readDepth']]
        df_amplicon.to_csv(outputfile, sep='\t', index=False, header=False)

    elif instrument == 'miseq':
        df_amplicon = pd.read_csv(inputfile, sep="\t", header=None, skiprows=1)
        df_amplicon.columns = ['ampliconName','readDepth']
        df_amplicon['gene'] = [x.split('.')[0] for x in df_amplicon['ampliconName'].values]
        df_amplicon['sampleID'] = sample
        df_amplicon = df_amplicon[['sampleID','gene','ampliconName','readDepth']]
        df_amplicon.to_csv(outputfile, sep='\t', index=False,header=False)

    elif instrument == 'nextseq':
        df_amplicon = pd.read_csv(inputfile, sep="\t", header=None)
        df_amplicon.columns = ['ampliconName','readDepth']
        df_amplicon['gene'] = [x.split('.')[0] for x in df_amplicon['ampliconName'].values]
        df_amplicon['sampleID'] = sample
        df_amplicon = df_amplicon[['sampleID','gene','ampliconName','readDepth']]
        df_amplicon.to_csv(outputfile, sep='\t', index=False,header=False)

try:
    parser = optparse.OptionParser()
    parser.add_option('-i', '--infile', help = 'provide input file')
    parser.add_option('-o', '--outfile', help = 'provide output file')
    parser.add_option('-s', '--sample', help = 'provide sampleID')
    parser.add_option('-n', '--instrument', help = 'provide instrumentID')
    options,args = parser.parse_args()
    infile = options.infile
    outfile = options.outfile
    sample = options.sample
    instrument = options.instrument
    loadIntoDB(infile, outfile, sample, instrument)

except TypeError:
	print ("python loadSampleAmplicons.py -help for help")
