import sys
import optparse
import re
import commonTool


def parseIonNewVarView(infile, outfile):
	infile = open(infile, 'r')
	outfile = open(outfile, 'w')
	for line in infile:
		if not line.startswith("#"):
			record = line.split("\t")
			if len(record[3]) == len(record[4]):
				type = "snv"
			else:
				type = "indel"
			vcf = commonTool.vcf(record)
			readDepth = vcf.parseInfo("FDP")
			altFreq = str(float(vcf.parseInfo("AF"))*100)
			altDepth = vcf.parseInfo("FAO")
			vep = vcf.parseInfo("CSQ").split('|')
			#parse dbsnp and cosmic from vep field
			IDs = vep[17].split("&")
			dbSNP = ""
			cosmic = ""
			for ID in IDs:
				if ID.startswith("rs"):
					if dbSNP == "":
						dbSNP = ID
					else:
						dbSNP += ",%s" %ID
				elif ID.startswith("COSM"):
					if cosmic == "":
						cosmic = ID
					else:
						cosmic += ",%s" %ID

			globalMinorAllele = ""
			globalMinorAleleFreq = ""
			alleleFreqAmr = ""
			alleleFreqAsn = ""
			alleleFreqAfr = ""
			alleleFreqEur = ""
			try:
				globalMinorAllele = vep[25].split(":")[0]
				globalMinorAleleFreq = str(float(vep[25].split(":")[1])*100)
				alleleFreqAmr = str(float(vep[27].split("&")[0].split(":")[1])*100)
				alleleFreqAsn = str(float(vep[28].split("&")[0].split(":")[1])*100)
				alleleFreqAfr = str(float(vep[26].split("&")[0].split(":")[1])*100)
				alleleFreqEur = str(float(vep[29].split("&")[0].split(":")[1])*100)
			except IndexError:
				pass
			#construct output
			output = [vep[3], vep[8], record[0], record[1], record[3], record[4], vep[2], type, record[5], altFreq,
			readDepth, altDepth, vep[1], vep[23], vep[24], vep[10], vep[11], dbSNP, vep[35]]

			outfile.write("\t".join(output))
			outfile.write("\n")
	infile.close()
	outfile.close()

def parseIlluminaNextseq(infile, outfile):
	infile = open(infile, 'r')
	outfile = open(outfile, 'w')
	for line in infile:
		if not line.startswith("#"):
			record = line.split("\t")
			if len(record[3]) == len(record[4]):
				type = "snv"
			else:
				type = "indel"
			vcf = commonTool.vcf(record)
			readDepth = vcf.parseInfo("DP")
			altFreq = str(float(vcf.parseDetail("VF", 10))*100)
			alleleDepth = vcf.parseDetail("AD", 10)
			altDepth = alleleDepth.split(",")[1].rstrip()#removes the newline character
			vep = vcf.parseInfo("CSQ").split('|')

			IDs = vep[17].split("&")

			dbSNP = ""
			cosmic = ""
			for ID in IDs:
				if ID.startswith("rs"):
					if dbSNP == "":
						dbSNP = ID
					else:
						dbSNP += ",%s" %ID
				elif ID.startswith("COSM"):
					if cosmic == "":
						cosmic = ID
					else:
						cosmic += ",%s" %ID
			#construct output
			output = [vep[3], vep[8], record[0], record[1], record[3], record[4], vep[2], type, record[5], altFreq,
			readDepth, altDepth, vep[1], vep[23], vep[24], vep[10], vep[11], dbSNP, vep[35]]

			outfile.write("\t".join(output))
			outfile.write("\n")
	infile.close()
	outfile.close()

try:

	if sys.argv[1] == "parseIonNewVarView":
		parser = optparse.OptionParser()
		parser.add_option('-I', '--vcfInfile',
		help = 'input vcf file')
		parser.add_option('-o', '--outfile',
		help = 'output file')
		options,args = parser.parse_args()
		infile = options.vcfInfile
		outfile = options.outfile
		parseIonNewVarView(infile, outfile)


	if sys.argv[1] == "parseIlluminaNextseq":
		parser = optparse.OptionParser()
		parser.add_option('-I', '--vcfInfile',
		help = 'input vcf file')
		parser.add_option('-o', '--outfile',
		help = 'output file')
		options,args = parser.parse_args()
		infile = options.vcfInfile
		outfile = options.outfile
		parseIlluminaNextseq(infile, outfile)
except IndexError:
	print "Usage:"
	print "python parseVEP.py parseIonNewVarView -h"
	print "python parseVEP.py parseIlluminaNextseq -h"
