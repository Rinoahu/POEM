#!usr/bin/env python
# this script is used to do cog annotation
import sys
import os
from Bio import SeqIO

##################
# the input files
##################
if len(sys.argv[1:]) < 2:
    print(
        'python this.py foo.fsa dir[contain prot2003-2014.fa, cog2003-2014.csv, cognames2003-2014.tab and fun2003-2014.tab] [blast]')
    raise SystemExit()


qry, pwd = sys.argv[1: 3]
try:
    blast_flag = sys.argv[3]
except:
    blast_flag = None
#ref, cids, names, funs = sys.argv[1:]
ref = pwd + '/prot2003-2014.fa'
cids = pwd + '/cog2003-2014.csv'
names = pwd + '/cognames2003-2014.tab'
funs = pwd + '/fun2003-2014.tab'

##############################
# build the ncbi to cog table
##############################
ncbi2cog = {}
if sys.version_info[0] > 2:
    f = open(cids, 'r', encoding='windows-1252')
else:
    f = open(cids, 'r')
for i in f:
    j = i[: -1].split(',')
    nbid = j[0]
    cgid = j[-3]
    ncbi2cog[nbid] = cgid

f.close()

# the cog to function table
#f = open(names, 'r')
if sys.version_info[0] > 2:
    f = open(names, 'r', encoding='windows-1252')
else:
    f = open(names, 'r')

cog2fun = {}
for i in f:
    j = i[: -1].split('\t')
    cog2fun[j[0]] = j[1:]

f.close()

###############################
# build the function to family
###############################
fun2fam = {}
f = open(funs, 'r')
for i in f:
    j = i[: -1].split('\t')
    fun2fam[j[0]] = j[1]

f.close()


#################################
# do blast and keep the best hit
#################################
blast_format = 'formatdb -i %s -p T'
blast_query = 'blastall -p blastp -i %s -d %s -o %s.blast8 -m 8 -a 16 -F F -e 1e-3 -b 1 -v 1'

diamond_format = 'diamond makedb --quiet --in %s -d %s'
diamond_query = 'diamond blastp --quiet -k 1 -q %s -d %s -e 1e-3 -p 16 -o %s.blast8'
if blast_flag == 'blast':
    # check the database for blastp
    if len(set([elem.split('.')[-1] for elem in os.listdir(pwd) if elem.endswith('.phr') or elem.endswith('.pin') or elem.endswith('.psq')])) != 3:
        os.system(blast_format%ref)
    os.system(blast_query%(qry, ref, qry))

elif blast_flag == 'diamond':
    # chekc the database for diamond
    if not os.path.isfile(ref + '.dmnd'):
        os.system(diamond_format%(ref, ref))
    os.system(diamond_query%(qry, ref, qry))
else:
    pass

# filter the blast results
os.system("awk '$3>=60&&$4>=60' %s.blast8 > %s.blast8.flt"%(qry, qry))

print('\t'.join(['gene', 'cog', 'annot']))

#f = open(qry + '.blast8', 'r')
f = open(qry + '.blast8.flt', 'r')
visit = set()
visit_qid = set()
outputs = set()
for i in f:
    m8 = i[: -1].split('\t')
    key = tuple(m8[: 2])
    j = m8[1].split('|')[1]
    if True:
        cgid = ncbi2cog[j]
        fun, fun_name = cog2fun.get(cgid, ['unknown', 'unknown'])

        if fun == fun_name == 'unknown':
            annot = ['unknown::unknown']
        else:
            annot = []
            for k in sorted(fun):
                fam = fun2fam[k]
                annot.append(k + '::' + fun_name + '::' + fam)

        ##############################################################
        # annot is the annotation for eache gene in following format:
        # K1::func1::fam1;;K2::func2:fam2
        ##############################################################
        annot = '$$'.join(annot)

        key = (m8[0], cgid)
        if key not in visit:
            print('\t'.join([m8[0], cgid, annot]))
            visit.add(key)
            visit_qid.add(m8[0])

f.close()

for i in SeqIO.parse(qry, 'fasta'):
    #if i.id not in visit:
    if i.id not in visit_qid:
        print '\t'.join([i.id, 'unknown', 'unknown::unknown' ])

print 'uniq gene', len(visit)
