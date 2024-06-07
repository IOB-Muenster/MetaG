#!/usr/bin/env python2


# SOURCE

# This script was originally obtained from
# https://github.com/adina/tutorial-ngs-2014/blob/master/ncbi/parse-genbank.py;
# commit: 1142c6e and was modified for this project. Consult the original source
# for copyright information.


import sys 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


n = 1

for record in list(SeqIO.parse(sys.argv[1], 'genbank')):
    id = record.id
    id = id.replace("_", "")
    for feat in record.features:
        if feat.type == "rRNA":
	    if 'product' in feat.qualifiers.keys():
            	if '16S ribosomal RNA' == feat.qualifiers['product'][0] or \
		'small subunit ribosomal RNA' in feat.qualifiers['product'][0] or \
		'Ribosomal RNA small subunit' == feat.qualifiers['product'][0] or \
		'16S ribosomal RNA (rr' in feat.qualifiers['product'][0] or \
		'16S RNA'  == feat.qualifiers['product'][0] or \
		'16S rRNA.' in feat.qualifiers['product'][0] or \
		'16SrRNA' in feat.qualifiers['product'][0] or \
		'Small Subunit Ribosomal RNA' in feat.qualifiers['product'][0] or \
		'ribosomal RNA 16S'  in feat.qualifiers['product'][0] or \
		'ribosomal RNA 16s'  in feat.qualifiers['product'][0] or \
		'Small subunit ribosomal RNA' in feat.qualifiers['product'][0]:
                	start = feat.location.start.position
                	end = feat.location.end.position
                	pos = [start, end]
                	print '>' + sys.argv[1].split('.')[0] + id + '_seq' + str(n)
                	print feat.extract(record.seq)
                	n += 1
