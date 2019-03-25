#!usr/bin/env python
# this script is used to find operon by graph theory
import sys
import os
import networkx as nx
from Bio import SeqIO
from collections import Counter
import cPickle as pickle
#from sklearn import cluster
import numpy as np


if len(sys.argv[1:]) < 1:
    print('python this.py foo.cog foo.list.adjancy [threshold]')
    raise SystemExit()

# get the script path
scriptpath='/'.join(os.path.realpath(__file__).split('/')[: -1]) + '/'


################## the graph for cog pair #########################################################
# build cog pair and annotation graph
fuc_dict = {}
f = open(scriptpath + '../database/cog/cog2014/fun2003-2014.tab', 'r')
for i in f:
    if i.startswith('#'):
        continue

    j = i[: -1].split('\t')[: 2]
    qid, fuc = j
    fuc_dict[qid] = fuc

f.close()

#f = open(scriptpath + '../database/cog/cog2014/cognames2003-2014.tab', 'r')
if sys.version_info[0] == 3:
    f = open(scriptpath + '../database/cog/cog2014/cognames2003-2014.tab', 'r', encoding='windows-1252')
else:
    f = open(scriptpath + '../database/cog/cog2014/cognames2003-2014.tab', 'r')


cog_dict = {}
cog_name_dict = {}
for i in f:
    j = i[:-1].split('\t')
    cog, fuc, name = j
    for k in fuc:
        try:
            cog_dict[k].add(cog)
        except:
            cog_dict[k] = set([cog])

    cog_name_dict[cog] = name

f.close()

flag = 0
# if cog id A and B shared at least 1 cog catergory, then store A and B in G_cog
G_cog = nx.Graph()
for i, j in cog_dict.items():
    for n1 in j:
        for n2 in j:
            flag += 1
            if n1 != n2:
                G_cog.add_edge(n1, n2, fuc = i, name = fuc_dict[i])

###################################################################################################

cog = sys.argv[1]
adj = sys.argv[2]
try:
    threshold = int(eval(sys.argv[3]))
except:
    threshold = 12


# the adjacent gene graph
G_adj = nx.Graph()
f = open(adj, 'r')
for i in f:
    j = i[: -1].split('\t')
    qid, sid, group = j[0], j[5], j[-1]
    if group == 'True':
        G_adj.add_edge(qid, sid)

f.close()



#  the graph
G = nx.Graph()
f = open(cog, 'r')
# elimate the header
f.next()
gene_cog_dict = {}
for i in f:
    j = i[: -1].split('\t')
    if len(j) != 3 or j[0] == 'gene':
        continue
    qid, cogid, annot = j
    # get the gene id, scaffold inf, cluster number, identity to the seed
    # during clustering.
    # get the family name: ABCDEF...
    # add the qid to the graph
    G.add_node(qid)
    try:
        gene_cog_dict[qid].add(cogid)
    except:
        gene_cog_dict[qid] = set([cogid])



# build the graph
#locus = G.nodes()
locus = list(G.nodes())

def sort_fuc(qid):
    #print 'debug', qid
    geneid, scaf, clstr, idy = qid.split('$$')
    strand, st, ed = geneid.split('|')[-3:]
    return scaf, strand, int(st), int(ed)

# construct cog graph
locus.sort(key=sort_fuc)
N = len(locus)
# the cog-eron like graph
cog_graph = nx.Graph()
# the gene graph
'''
gene_graph = nx.Graph()
'''

cog_graph_node_freq = Counter()
for i in gene_cog_dict.values():
    for j in i:
        if j != 'unknown':
            cog_graph_node_freq[j] += 1

cog_graph_edge_freq = Counter()


# determine cog id A and B share the same cog catergory
def find_share(acogs, bcogs, G_cog):
    for acog in acogs:
        for bcog in bcogs:
            if G_cog.has_edge(acog, bcog):
                #fuc = G_cog.edge[acog][bcog].get('fuc', 'unknown')
                fuc = G_cog[acog][bcog].get('fuc', 'unknown')
                #name = G_cog.edge[acog][bcog].get('name', 'unknown')
                name = G_cog[acog][bcog].get('name', 'unknown')
                share = fuc + '::' + name
                return acog, bcog, name, share
    return None, None, 'unknown', 'unknown::unknown'


'''
print '# the predicted operon gene adjacency'
print '\t'.join(['geneA', 'geneB', 'annotAB'])
'''

visit_AB = set()
for i in xrange(N - 1):
    a, b = locus[i: i + 2]
    a_scaf, a_strand, a_st, a_ed = sort_fuc(a)
    b_scaf, b_strand, b_st, b_ed = sort_fuc(b)

    A = '$$'.join(a.split('$$')[: 2])
    B = '$$'.join(b.split('$$')[: 2])


    # add edge between a and b if:
    # 1. the dist of a-b is < 450
    # 2. a, b have same function id of cog
    acogs = gene_cog_dict.get(a, set())
    bcogs = gene_cog_dict.get(b, set())

    '''
    share = None

    #print 'a, b', a, b, G.node[a], G.node[b], 'a, b cog', acog, bcog
    for acog in acogs:
        for bcog in bcogs:
            # the a, b should have same cog catergory
            #if G_cog.has_edge(acog, bcog):
            if 1:
                try:
                    fuc = G_cog.edge[acog][bcog].get('fuc', 'unknown')
                except:
                    fuc = 'unknown'
                try:
                    name = G_cog.edge[acog][bcog].get('name', 'unknown')
                except:
                    name = 'unknown'
                share = fuc + '::' + name
                if G_adj.has_edge(A, B):
                    ab_funid = share
                    cog_graph.add_edge(acog, bcog, fuc = fuc, name = name, annot = share)
                    key = tuple(sorted([acog, bcog]))
                    cog_graph_edge_freq[key] += 1
                    gene_graph.add_edge(a, b)
    '''

    if G_adj.has_edge(A, B):
        acog, bcog, name, share = find_share(acogs, bcogs, G_cog)

        '''
        if (a, b) not in visit_AB:
            print '\t'.join([a, b, name])
            visit_AB.add((a, b))
        '''

        if acog and bcog:
            cog_graph.add_edge(acog, bcog, fuc = fuc, name = name, annot = share)
            key = tuple(sorted([acog, bcog]))
            cog_graph_edge_freq[key] += 1
        '''
        gene_graph.add_edge(a, b)
        '''

    else:
        continue

# stat of cog freq
edge_freq = Counter()
cog_freq = Counter()
gene_visit = set()
for edge in G.edges():
    a, b = edge
    a_cog, b_cog = G.node[a]['cog'], G.node[b]['cog']
    if a not in gene_visit:
        cog_freq[a_cog] += 1
        gene_visit.add(a)
    if b not in gene_visit:
        cog_freq[b_cog] += 1
        gene_visit.add(b)
    edge_cog = tuple(sorted([a_cog, b_cog]))
    edge_freq[edge_cog] += 1

total_edge = sum(cog_graph_edge_freq.values())

'''
print '#' * 30
print 'Graph Information'
print '#' * 30
print 'total_edge', total_edge, sum(cog_graph_edge_freq.values())
total_cog = sum(cog_graph_node_freq.values())
print 'total_cog', total_cog
'''

G_AB = nx.Graph()
print '#' * 30
print '# COG adjacency'
print '#' * 30

print 'geneA', 'geneB', 'cA', 'cB', 'cAB', 'cA+cB-cAB', 'cAB/(cA+cB-cAB)'
for i, j in cog_graph_edge_freq.items():
    a, b = i
    a_count = cog_graph_node_freq[a]
    b_count = cog_graph_node_freq[b]
    ab_count = j
    Pa_and_b = ab_count * 1. / total_edge
    Pa_or_b = (a_count + b_count) * 1. / total_edge - Pa_and_b
    if j >= threshold and a != 'unknown' and b != 'unknown':
        try:
            res = ab_count * 1.0 / (a_count + b_count - ab_count)
        except:
            res = 1.
        print a, b, a_count, b_count, ab_count, a_count + b_count - ab_count, res

        G_AB.add_edge(a, b, weight = ab_count)


if sys.version_info[0] == 3:
    f = open(scriptpath + '../database/operondb/operon_graph.txt', 'rb')
else:
    f = open(scriptpath + '../database/operondb/operon_graph.txt', 'r')


#paths, gene_cog_dict, cog_annot_dict = pickle.load(f)
if sys.version_info[0] > 2:
    paths, gene_cog_dict, cog_annot_dict = pickle.load(f, encoding='latin1')
else:
    paths, gene_cog_dict, cog_annot_dict = pickle.load(f)

f.close();

operons = [[gene_cog_dict.get(b, set()) for b in a ] for a in paths]


# find the best targe operon that share the most gene with predicted operon
def findbest(Cogs, operons, paths):
    best_c = 0
    best_o = 0
    cogs = list(Cogs)
    M = len(cogs)
    for op, path in zip(operons, paths):
        N = len(op)
        Ns = [0] * N
        Ms = [0] * M
        for a in xrange(M):
            for b in xrange(N):
                if cogs[a] in op[b]:
                    Ns[b] = 1
                    Ms[a] = 1
        score_c  = sum(Ms)
        score_o = sum(Ns)
        if score_o >= best_o:
            best_c = score_c
            best_o = score_o
            best_op = op
            best_path = '$$'.join(path)

    return best_op, best_c * 1. / len(cogs), best_o * 1. / len(best_op), best_path


# extract the component from graph
def get_comp(G_AB, init = 3):
    for i in nx.connected_components(G_AB):
        if len(i) > 30:
            print 'iteration init', init
            G_AB_sub = G_AB.subgraph(i)
            big = set(sum([[a, b] for a, b in G_AB_sub.edges() if G_AB_sub[a][b]['weight'] > threshold], []))
            if len(big) == 0:
                yield i
            elif 0 < len(big) < len(G_AB_sub.nodes()):
                G_AB_small = G_AB_sub.subgraph([elem for elem in G_AB_sub.nodes() if elem not in big])
                G_AB_big = G_AB_sub.subgraph(big)
                get_comp(G_AB_small, threshold)
                get_comp(G_AB_big, threshold + 1)
            else:
                get_comp(G_AB_sub, threshold + 1)
        else:
            yield i


# filter node with too many neigbors
#print 'G_AB is', G_AB, G_AB.degree()
#outdeg = G_AB.degree()
outdeg = dict(G_AB.degree())
x = outdeg.values()
max_dg = np.mean(x) + np.std(x) * 2
G_AB_sub = G_AB.subgraph([elem for elem in outdeg if outdeg[elem] <= max_dg])


#print '#' * 30
#print 'Core COG adjacency'
#print '#' * 30
header = [ 'geneA', 'geneB', 'cA', 'cB', 'cAB', 'cA+cB-cAB', 'cAB/(cA+cB-cAB)']
# the adjacency of core operon
_o = open(cog[: -len('.cog')] + '.core_cog_adjacency', 'w')
_o.write('\t'.join(header) + '\n')
for i in get_comp(G_AB_sub, threshold):
    g_sub = G_AB.subgraph(i)
    for a, b in g_sub.edges():
        a_count = cog_graph_node_freq[a]
        b_count = cog_graph_node_freq[b]
        key = a < b and (a, b) or (b, a)
        ab_count = cog_graph_edge_freq[key]
        #Pa_and_b = ab_count * 1. / total_edge
        #Pa_or_b = (a_count + b_count) * 1. / total_edge - Pa_and_b
        try:
            res = ab_count * 1.0 / (a_count + b_count - ab_count)
        except:
            res = 1.

        #print a, b, a_count, b_count, ab_count, a_count + b_count - ab_count, res, 'True'
        output = [a, b, a_count, b_count, ab_count, a_count + b_count - ab_count, res]
        output = map(str, output)
        _o.write('\t'.join(output) + '\n')

_o.close()



print '\n'
print '#' * 30
print '# evalute the core operon'
print '#' * 30
#print 'predict_operon\treal_operon\tSP\tSN\tF1'
print 'predict_operon\treal_operon_name\treal_operon_cog\tPrecise\tRecall\tF1\tpredict_annot\treal_annot'

for i in get_comp(G_AB_sub, threshold):
    if len(i) < 2:
        continue
    best_op, SP, SN, best_path = findbest(i, operons, paths)

    #print 'predict', i, 'real_operon', best_op, 'SP', SP, 'SN', SN
    #predict_annot = '$$'.join(sorted(set([cog_name_dict.get(elem, 'unknown') for elem in i])))
    predict_annot = '$$'.join([cog_name_dict.get(elem, 'unknown') for elem in i])

    #real_annot = '$$'.join(sorted(set(sum([[cog_name_dict.get(elem2, 'unknown') for elem2 in elem] for elem in best_op], []))))
    real_annot = '$$'.join(sum([[cog_name_dict.get(elem2, 'unknown') for elem2 in elem] for elem in best_op], []))


    #best_op = '$$'.join([';;'.join(elem) for elem in best_op])
    best_ops = [';;'.join(elem) for elem in best_op]
    best_ops = [elem and elem or '*' for elem in best_ops]
    best_op = '$$'.join(best_ops)

    if SN < .2 or SP < .2:
        best_path = best_op = 'not found'

    F1 =  SP * SN > 0 and 2. * SP * SN / (SP + SN) or 0
    output = '%s\t%s\t%s\t%f\t%f\t%f'%('$$'.join(i), best_path, best_op, SP, SN, F1)
    print output + '\t' + predict_annot + '\t' + real_annot

