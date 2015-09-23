##########################################################################################
## to use in binning pandakiller assembled samples into seqfilter parsed tag directories############################################################################################
import os
import sys
from Bio import SeqIO
import string

assemb = open(sys.argv[1], 'rU')
dict_assemb = {}
for records in SeqIO.parse(assemb, "fastq"):
	records.id = records.id.rsplit(':', 1)[0]
	records.description = records.description.rsplit(':', 1)[0]
	records.name = records.name.rsplit(':', 1)[0]
	#print records
	dict_assemb[records.id] = records

rootdir = "."
lfiles = []
for topdir, subdir, filelist in os.walk(rootdir):
	for files in filelist:
		if files.endswith("trimmed.fasta"):
			lfiles.append(files)

for each in lfiles:
	fout = '%s_assem.fastq' % each.strip()
        fp = open(fout, 'w')
	f = open(each, 'rU')
	l = []
	for records in SeqIO.parse(f, 'fasta'):
		name = records.id.split(' ')[0]
		if name in dict_assemb:
			#print each + ':' + name
			l.append(name)
#	print len(l)
	for item in l:
		SeqIO.write(dict_assemb[item], fp, 'fastq')
#		fp.write('>%s\n' % item) 
#		fp.write('%s\n' % dict_assemb[item])
#	for files in filelist:
#		if files.endswith("trimmed.fasta"):
#			lfiles.append(files)
#			for each in lfiles:
#				outfiles_keep = open(each + '_' + 'assem' + '.fasta', 'w')
				#outfiles_drop = open(each + '_' + 'dropped' + 'reads' + '.txt', 'w')
#				f = open(each, 'rU')
#				names = []
#				for records in SeqIO.parse(f, 'fasta'):
#					names.append(records.id.split(' ')[0])
#				for item in names:
#					if item in dict_assemb:	
#						outfiles_keep.write('>%s\n' % (item))
#						outfiles_keep.write('%s\n' % (dict_assemb[item]))
					#else:
					#	outfiles_drop. write('%s\n' % (item))
