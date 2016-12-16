#!usr/bin/env python
import sys

if len(sys.argv[1:]) < 1:
    print('python this.py metagenemark_aa.fsa')
    raise SystemExit()


qry = sys.argv[1]
f = open(qry, 'r')
for i in f:
    if i.startswith('>'):
        print i[: -1].replace('\t>', '$$').replace(' ', '-')
    else:
        print i[: -1]
