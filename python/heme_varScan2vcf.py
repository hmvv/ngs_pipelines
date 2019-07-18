import sys
import optparse
import re


def convert(input, output):
	infile = open(input, 'r')
	outfile = open(output, 'w')
	outfile.write('##fileformat=VCFv4.1\n')
	outfile.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTE\tINFO\tFORMAT\tsample\n')
	thisPosition = []
	preChr = ''
	prePos = '0'
	for line in infile:
		line = line.split()
		chr = line[0]
		pos = line[1]
		ref = line[2]
		alt = line[3].strip('-')
		refC = int(line[4])
		altC = int(line[5])
		refQ = int(line[6])
		altQ = int(line[7])
		depth = refC+altC
		altF = float(altC)/float(depth)
		AD = '%i,%i' %(refC, altC)
		quality = int(float(refC*refQ + altC*altQ)/float(depth))
		sampleInfo = '0/1:%f:%i:%s' %(altF, depth, AD)
		info = 'DP=%i' %depth
		out = '%s\t%s\t.\t%s\t%s\t%i\tPASS\t%s\tGT:VF:DP:AD\t%s' %(chr,pos,ref,alt,quality,info,sampleInfo)
		#print out
		if chr == preChr and pos == prePos:
			#print "same as previous chr: %s and previous position: %s" %(preChr, prePos)
			#print thisPosition
			if (chr,pos,ref,alt) in thisPosition:
				#print "same ref and alt"
				pass
			else:
				#print "different ref and alt"
				thisPosition.append((chr,pos,ref,alt))
				outfile.write(out)
				outfile.write('\n')
		else:
			#print "different from previous chr: %s and previous position: %s" %(preChr, prePos)
			thisPosition = [(chr, pos, ref, alt)]
			preChr = chr
			prePos = pos
			outfile.write(out)
			outfile.write('\n')
	infile.close()
	outfile.close()


try:
	
	parser = optparse.OptionParser()
	parser.add_option('-I', '--input',
	help = 'input varscan output file')
	parser.add_option('-o', '--output',
	help = 'output vcf file')
	options,args = parser.parse_args()
	input = options.input
	output = options.output
	convert(input, output)	
	
except TypeError:
	print 'python varScan2vcf.py -help for help'