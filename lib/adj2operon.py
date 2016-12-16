#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
# CreateTime: 2016-12-14 11:35:45

# convert the operonic adjacency to full operon
import sys

if len(sys.argv[1: ]) < 2:
    print 'python this.py foo.adjacency foo.cog'
    raise SystemExit()


# extract full operon from operonic adjacency
def extract_operon(f):
    operon = []
    for i in f:
        j = i[:-1].split('\t')
        qid, strand, sid, label = j[0], j[2], j[5], j[-1]
        if label == 'False':
            continue

        if operon and operon[-1][1] == qid and operon[-1][2] == strand:
            operon.append([qid, sid, strand])
        else:
            if operon:
                yield [elem[0] for elem in operon] + [operon[-1][1]], operon[0][2]
            operon = [[qid, sid, strand]]

    if operon:
        yield [elem[0] for elem in operon] + [operon[-1][1]], operon[0][2]



cog = sys.argv[2]
f = open(cog, 'r')
cog_annot_dict = {}
for i in f:
    j = i[: -1].split('\t')
    if len(j) == 3:
        qid = j[0].split('$$')[0]
        cog_annot_dict[qid] = j[1: ]

f.close()


qry = sys.argv[1]

header = ['gene_id', 'COG_id', 'annotation']
print '\t'.join(header)
f = open(qry, 'r')
for i, j in extract_operon(f):
    #print i, j
    sep = j == '+' and '-->' or '<--'
    cog_annot = [cog_annot_dict.get(elem.split('$$')[0], ['unknown', 'unknown::unknown']) for elem in i]
    cog = [elem[0] for elem in cog_annot]
    annot = [elem[1] for elem in cog_annot]
    print '\t'.join(map(sep.join, [i, cog, annot]))
f.close()
