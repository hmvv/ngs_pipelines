import sys
import optparse
import commonTool

def splitVcf(infile, outfile, fields):
	infile = open(infile, 'r')
	outfile = open(outfile, 'w')
	fields = fields.split(",")
	outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
	for line in infile:
		if not line.startswith("#"):
			record = line.split("\t")
			dic = fieldDic(record, fields)
			vcf = commonTool.vcf(record)
			if not "," in record[4]:
				output = record[0:7]
				info = ""
				for field in fields:
					info += "%s=%s;" %(field, dic[field])
				info = info.strip(";")
				output.append(info)
				outfile.write("\t".join(output))
				outfile.write("\n")
			else:
				altAlleles = record[4].split(",")
				n = 0
				while n < len(altAlleles):
					output = record[0:4] + [altAlleles[n]] + record[5:7]
					info = ""
					for field in fields:
						if "," in dic[field]:
							info += "%s=%s;" %(field, dic[field].split(",")[n])
						else:
							info += "%s=%s;" %(field, dic[field])
					info = info.strip(";")
					output.append(info)
					newVcf = commonTool.vcf(output)
					outputLeft = newVcf.leftAlign()
					outfile.write("\t".join(outputLeft))
					outfile.write("\n")
					n += 1
	infile.close()
	outfile.close()
					
def fieldDic(record, fields):
	dic = {}
	vcf = commonTool.vcf(record)
	for field in fields:
		value = vcf.parseInfo(field)
		dic[field] = value
	return dic


if len(sys.argv) > 1:
	parser = optparse.OptionParser()
	parser.add_option('-I', '--infile',
	help = 'input vcf file')
	parser.add_option('-o', '--outfile',
	help = 'output file')
	parser.add_option('-f', '--fields',
	help = 'fields to include in the output, seperated by ,')
	options,args = parser.parse_args()
	infile = options.infile
	outfile = options.outfile
	fields = options.fields
	splitVcf(infile, outfile, fields)
else:
	print "split multiple alt alleles into seperate lines.  Only the first 8 columns are retained in output file"
	print "For more Usage info: python splitVcf.py -h"