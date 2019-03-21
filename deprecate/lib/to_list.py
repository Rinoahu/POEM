#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
# CreateTime: 2016-12-05 13:20:06

import sys

if len(sys.argv[1:]) < 1:
    print 'python this.py genemark_reid_aa.fasta'
    raise SystemExit()


qry = sys.argv[1]
f = open(qry, 'r')
for i in f:
    if i.startswith('>'):
        j = i[1: -1]
        loc, tax = j.split('$$')
        k = loc.split('|')
        strand, start, end = k[-3: ]
        if '-flag=' in tax and '-multi=' in tax and '-len=' in tax:
            tax = tax.split('-flag=')[0]
        print '\t'.join([j, tax, strand, start, end])

f.close()
