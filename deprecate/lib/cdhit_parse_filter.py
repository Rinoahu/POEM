#!usr/bin/env python
# this script is used to do cd-hit cluster and parse the results
import sys
from Bio import SeqIO

if len(sys.argv[1:]) < 1:
    print 'python this.py foo.pep.fsa.cdhit.clstr foo.pep.fsa'
    raise SystemExit()


qry = sys.argv[1]
fasta = sys.argv[2]
f = open(qry, 'r')

# parse the cdhit cluster file
# the output is list
# 1st elem of list is cluster number
# rest of list is sequence length, id, etc


def cdhit_parse(f):
    output = []
    for i in f:
        if i.startswith('>'):
            if output:
                yield output
            output = [i[: -1]]
        else:
            output.append(i[: -1])

    if output:
        yield output

# filter the ortholog with 70 <= id <= 98

def get_id(elem):
    return float(elem.split('... at')[1].split('%')[0])

# build the cdhit filter table
cdhit_dict = {}
for i in cdhit_parse(f):
    output = []
    for j in i:
        if '... at' not in j or 75 <= get_id(j) <= 98:
            output.append(j)
        else:
            continue
    # keep gene have 3 ortholog at least
    if len(output) >= 1:

        group = output[0].replace(' ', '_')
        for k in output[1:]:
            qid, feature = k.split('>')[1].split('... ')
            feature = feature.replace(' ', '_')
            cdhit_dict[qid] = [group, feature]


# print cdhit_dict
# keep the filtered aa seq
seqs = SeqIO.parse(fasta, 'fasta')
for i in seqs:
    if i.id in cdhit_dict:
        group, feature = cdhit_dict[i.id]
        print '>' + i.description + '$$' + group + '$$' + feature
        print str(i.seq)
    else:
        print '>' + i.description + '$$' + 'Orphan' + '$$' + 'Orphan'
        print str(i.seq)

