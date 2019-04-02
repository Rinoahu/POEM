#!usr/bin/env python
# change the header of genome for prokka

import sys

qry = sys.argv[1]


flag = 0
f = open(qry, 'r')
for i in f:
    if i.startswith('>'):
        j = i[1:-1]
        print '>AAA_%d'%flag + ' xxxxxx' + j
    else:
        print i[:-1]

    flag += 1
f.close()

