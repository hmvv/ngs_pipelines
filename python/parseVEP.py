import sys
import optparse
import re
import commonTool


def parseHeme(infile, outfile):
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
			altDepth = vcf.parseDetail("AD", 10).split(",")[1]
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
						cosmic += ",%s" %cosmic
			
			globalMinorAllele = ""
			globalMinorAleleFreq = ""
			alleleFreqAmr = ""
			alleleFreqAsn = ""
			alleleFreqAfr = ""
			alleleFreqEur = ""
			try:
				globalMinorAllele = vep[27].split(":")[0]
				globalMinorAleleFreq = str(float(vep[27].split(":")[1])*100)
				alleleFreqAmr = str(float(vep[29].split("&")[0].split(":")[1])*100)
				alleleFreqAsn = str(float(vep[30].split("&")[0].split(":")[1])*100)
				alleleFreqAfr = str(float(vep[28].split("&")[0].split(":")[1])*100)
				alleleFreqEur = str(float(vep[31].split("&")[0].split(":")[1])*100)
			except IndexError:
				pass
			#construct output
			output = [vep[3], record[0], record[1], record[3], record[4], vep[2], type, vep[19], record[5], altFreq,
			readDepth, altDepth, vep[1], vep[24], vep[25], vep[10], vep[11], dbSNP, "", globalMinorAleleFreq, globalMinorAllele, 
			alleleFreqAmr, alleleFreqAsn, alleleFreqAfr, alleleFreqEur, cosmic, "", ""]
			
			outfile.write("\t".join(output))
			outfile.write("\n")
	infile.close()
	outfile.close()
			
			

def parseIon(infile, outfile):
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
						cosmic += ",%s" %cosmic
			
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
			readDepth, altDepth, vep[1], vep[22], vep[23], vep[10], vep[11], dbSNP, globalMinorAleleFreq, globalMinorAllele, 
			alleleFreqAmr, alleleFreqAsn, alleleFreqAfr, alleleFreqEur]
			
			outfile.write("\t".join(output))
			outfile.write("\n")
	infile.close()
	outfile.close()
	
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
						cosmic += ",%s" %cosmic
			
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
			readDepth, altDepth, vep[1], vep[22], vep[23], vep[10], vep[11], dbSNP, vep[34]]
			
			outfile.write("\t".join(output))
			outfile.write("\n")
	infile.close()
	outfile.close()

def parseIllumina(infile, outfile):
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
			altDepth = alleleDepth.split(",")[1]
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
						cosmic += ",%s" %cosmic
			#construct output
			output = [vep[3], vep[8], record[0], record[1], record[3], record[4], vep[2], type, record[5], altFreq,
			readDepth, altDepth, vep[1], vep[22], vep[23], vep[10], vep[11], dbSNP, vep[34]]
			
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
			
			#print vep
			
			#parse dbsnp and cosmic from vep field
			#CSQ=G|intron_variant|MODIFIER|CSF3R|ENSG00000119535|Transcript|ENST00000373103|protein_coding||15/16|ENST00000373103.1:c.1959-107T>C|||||||||-1|HGNC|2439|||||||||||||	GT:GQ:AD:VF:NL:SB:GQX	0/1:59:827,26:0.0304:20:-30.0099:59
			#CSQ=-|intron_variant|MODIFIER|UBE2N|ENSG00000177889|Transcript|ENST00000318066|protein_coding||1/3|ENST00000318066.2:c.30+3207delc|||||||||-1|HGNC|12492|||||||||||||	GT:VF:DP:AD	0/1:0.047619:42:40,2
			#CSQ=T|intergenic_variant|MODIFIER|||||||||||||||rs749868180||||||||T:0.0130|T:0.0047|T:0.0265|T:0.0000|T:0.0083|T:0.0000||||	GT:VF:DP:AD	0/1:1.000000:37:0,37
			#CSQ=C|intron_variant|MODIFIER|PTEN|ENSG00000171862|Transcript|ENST00000371953|protein_coding||5/8|ENST00000371953.3:c.492+14T>C|||||||COSM14247||1|HGNC|9588|||||||||||1|1|	GT:VF:DP:AD	0/1:0.010707:2802:2772,30
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
						cosmic += ",%s" %cosmic
			#construct output
			output = [vep[3], vep[8], record[0], record[1], record[3], record[4], vep[2], type, record[5], altFreq,
			readDepth, altDepth, vep[1], vep[22], vep[23], vep[10], vep[11], dbSNP, vep[34]]
			
			outfile.write("\t".join(output))
			outfile.write("\n")
	infile.close()
	outfile.close()
	
try:
	if sys.argv[1] == "parseHeme":
		parser = optparse.OptionParser()
		parser.add_option('-I', '--vcfInfile',
		help = 'input vcf file')
		parser.add_option('-o', '--outfile',
		help = 'output file')
		options,args = parser.parse_args()
		infile = options.vcfInfile
		outfile = options.outfile
		parseHeme(infile, outfile)
	
	if sys.argv[1] == "parseIon":
		parser = optparse.OptionParser()
		parser.add_option('-I', '--vcfInfile',
		help = 'input vcf file')
		parser.add_option('-o', '--outfile',
		help = 'output file')
		options,args = parser.parse_args()
		infile = options.vcfInfile
		outfile = options.outfile
		parseIon(infile, outfile)
	
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
		
	if sys.argv[1] == "parseIllumina":
		parser = optparse.OptionParser()
		parser.add_option('-I', '--vcfInfile',
		help = 'input vcf file')
		parser.add_option('-o', '--outfile',
		help = 'output file')
		options,args = parser.parse_args()
		infile = options.vcfInfile
		outfile = options.outfile
		parseIllumina(infile, outfile)
		
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
	print "python parseVEP.py parseHeme -h"
	print "python parseVEP.py parseIon -h"
	print "python parseVEP.py parseIonNewVarView -h"
	print "python parseVEP.py parseIllumina -h"
	print "python parseVEP.py parseIlluminaNextseq -h"
	
