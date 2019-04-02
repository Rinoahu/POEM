#!usr/bin/env python
# this script is used to do cog annotation
import sys
import os
from Bio import SeqIO

##################
# the input files
##################
if len(sys.argv[1:]) < 2:
    print ('this script is used to convert gene identifier of prokka to genemark.\nusage:')
    print('python this.py foo.genome.fsa foo.faa foo.gff')
    raise SystemExit()

genome, qry, gff = sys.argv[1:4]

scf_dict = {}
for i in SeqIO.parse(genome, 'fasta'):
    qid = i.id
    qid_c = qid.replace('|', '_')
    scf_dict[qid_c] = qid


id_dict = {}
# parse gff
f = open(gff, 'r')
for i in f:
    j = i[:-1].split('\t')
    if len(j) != 9:
        continue
    scf, met, start, end, strand, qid = j[0], j[1], j[3], j[4], j[6], j[8]
    try:
        scf = scf_dict[scf]
    except:
        continue
    qid = qid.split(';')[0][3:].strip()
    length = str(abs(int(end) - int(start) // 3)) + '_aa'
    desc = ['>', qid, '|', met, '|', length, '|', strand, '|', start, '|', end, '\t', '>', scf]
    id_dict[qid] = ''.join(desc)

f.close()

seqs = SeqIO.parse(qry, 'fasta')

for i in seqs:
    if i.id in id_dict:
        print id_dict[i.id]
        print str(i.seq)












