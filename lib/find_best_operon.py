#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#  Copyright © XYM
# CreateTime: 2017-02-24 16:08:49

import pickle

f = open('./operon_graph.txt', 'r')
paths, gene_cog_dict, cog_annot_dict = pickle.load(f)

operons = [[gene_cog_dict.get(b, set()) for b in a ] for a in paths]


# find the best targe operon that share the most gene with predicted operon
def findbest1(Cogs, operons, paths):
    best_c = 0
    best_o = 0
    cogs = list(Cogs)
    M = len(cogs)
    output = []
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
        #if score_o >= best_o:
        #if score_o > best_o:
        if 1:
            best_c = score_c
            best_o = score_o
            best_op = op
            best_path = '$$'.join(path)
            out = [best_op, best_c * 1. / len(cogs), best_o * 1. / len(best_op), best_path, score_o]
            output.append(out)

    output.sort(key = lambda x: x[-1], reverse = True)

    return output
