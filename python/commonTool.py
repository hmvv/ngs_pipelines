import sys
import optparse
import re

aaDic = {"Ala":"Alanine", "Ile":"Isoleucine", "Leu":"Leucine", "Val":"Valine",
"Phe":"Phenylalanine", "Trp":"Tryptophan", "Tyr":"Tyrosine", "Asn":"Asparagine",
"Met":"Methionine", "Cys":"Cysteine", "Ser":"Serine", "Gln":"Glutamine", 
"Thr":"Threnine", "Asp":"Aspartic acid", "Glu":"Glutamic acid", "Arg":"Arginine",
"His":"Histidine", "Lys":"Lysine", "Gly":"Glycine", "Pro":"Proline"}

def hgvs_protein_parser(tag):
	#print tag
	if tag == "":
		output = "no protein coding change"
	else:
		original_aa = re.compile(r'p\.([A-Z][a-z][a-z])([0-9]+)')
		sub_aa = re.compile(r'p\.([A-Z][a-z][a-z])([0-9]+)([A-Z][a-z][a-z])$')
		truncation = re.compile(r'p\.([A-Z][a-z][a-z])([0-9]+)\*$')
		fs = re.compile(r'p\.([A-Z][a-z][a-z])([0-9]+)fs$')
		range_deletion = re.compile(r'p\.([A-Z][a-z][a-z])([0-9]+)_([A-Z][a-z][a-z])([0-9]+)del$')
		range_deletion_ins = re.compile(r'p\.([A-Z][a-z][a-z])([0-9]+)_([A-Z][a-z][a-z])([0-9]+)delins([A-Z][a-z][a-z])$')
		single_deletion = re.compile(r'p\.([A-Z][a-z][a-z])([0-9]+)del$')
		original_aa_name = original_aa.search(tag).group(1)
		original_aa_name = aaDic[original_aa_name]
		original_aa_number = original_aa.search(tag).group(2)
		if sub_aa.search(tag):
			#simple substitution
			sub_aa_name = sub_aa.search(tag).group(3)
			sub_aa_name = aaDic[sub_aa_name]
			output = "%s at amino acid %s changed to %s " % (original_aa_name, original_aa_number, sub_aa_name)
		elif fs.search(tag):
			#frame shift
			output = "frameshift starting from %s at amino acid %s" % (original_aa_name, original_aa_number)
		elif truncation.search(tag):
			#truncation
			output = "truncation at amino acid %s at %s" % (original_aa_name, original_aa_number)
		elif range_deletion.search(tag):
			#more than one aa deleted
			del_aa_name = range_deletion.search(tag).group(3)
			del_aa_name = aaDic[del_aa_name]
			del_aa_number = range_deletion.search(tag).group(4)
			output = "deletion from %s at %s to %s at %s" % (original_aa_name, original_aa_number, del_aa_name, del_aa_number)
		elif range_deletion_ins.search(tag):
			#more than one aa deleted and then aa inserted
			del_aa_name = range_deletion.search(tag).group(3)
			del_aa_name = aaDic[del_aa_name]
			del_aa_number = range_deletion.search(tag).group(4)
			ins_aa_name = aaDic[ins_aa_name]
			output = "deletion from %s at %s to %s at %s, insertion of %s in place" % (original_aa_name, original_aa_number, del_aa_name, del_aa_number, ins_aa_name)
		elif single_deletion.searth(tag):
			#single aa deletion
			output = "deletion of %s at %s" % (original_aa_name, original_aa_number)
		else:
			#not included in the above situations
			output = "protein change not recognized"
	return(output)
	
class fasta:
	def __init__(self, fastaFile):
		#take a fasta file as class parameter
		self.handle = open(fastaFile, 'r')
	def nextRecord(self):
		seq = ''
		header = None
		previous = 0
		line = self.handle.readline()
		while line:
			#print "current line: %s current pointer %i" %(line,previous)
			if line.startswith('>') and not header:
				#print "new header"
				#print line
				header = line.strip()
			elif not line.startswith('>') and header:
				#print "appending sequences"
				#print line
				seq += line.strip()
			elif line.startswith('>') and header:
				#print "next header"
				#print line
				self.handle.seek(previous)
				break
			previous = self.handle.tell()
			line = self.handle.readline()
		return(header, seq)
	def random(self, length, seq):
		import random
		totalLength = len(seq)
		start = random.randint(0, totalLength - length)
		return seq[start:start+length+1]
		
class dna:
	def __init__(self, sequence):
		self.sequence = sequence
	def reverseComplement(self):
		reverseSequence=self.sequence[::-1].upper()
		outSequence = ''
		for base in reverseSequence:
			if base == "A":
				outSequence += "T"
			elif base == "T":
				outSequence += "A"
			elif base == "C":
				outSequence += "G"
			elif base == "G":
				outSequence += "C"
			else:
				outSequence += "N"
		return outSequence
		
		
class vcf:
	def __init__(self, oneVcf):
		self.record = oneVcf
	def parseInfo(self, infoID):
		record = self.record
		info = record[7].split(';')
		infoDic={}
		for field in info:
			try:
				id = field.split('=')[0]
				value = field.split('=')[1].strip()
			except IndexError:
				continue
			if id not in infoDic:
				infoDic[id] = value
		try:
			return infoDic[infoID]
		except KeyError:
			#print "specified info field not in info column"
			return "null"
	def parseDetail(self, detailID, detailCol):
		record = self.record
		detailCol = int(detailCol) - 1
		detail = record[detailCol].split(':')
		format = record[8]
		formatList = format.split(':')
		if not detailID in formatList:
			#print "specified detail ID not exists"
			return "null"
		else:
			index = formatList.index(detailID)
			return detail[index]
	def leftAlign(self):
		record = self.record
		ref = record[3]
		alt = record[4]
		left = 0
		##truncate from the right end
		while len(ref) > 1 and len(alt) > 1:
			##truncate from the right end
			if ref[-1] == alt[-1]:
				ref = ref[0:-1]
				alt = alt[0:-1]
			else:
				break
		##truncate from the 
		while len(ref) > 1 and len(alt) > 1:
			##truncate from the right end
			if ref[0] == alt[0]:
				ref = ref[1:]
				alt = alt[1:]
				left += 1
			else:
				break
		record[1] = str(int(record[1]) + left)
		record[3] = ref
		record[4] = alt
		return record
			
			
			
		