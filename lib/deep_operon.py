#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
# CreateTime: 2016-09-21 16:51:48

import numpy as np
from Bio import SeqIO, Seq, SeqUtils
#from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex
from Bio.SeqUtils import GC
from Bio.SeqUtils.CodonUsage import SynonymousCodons
import math
from math import log, sqrt
from collections import Counter
import pickle

# from sklearn import cross_validation, metrics  # Additional scklearn
# functions
#from sklearn import metrics  # Additional scklearn functions
# from sklearn.grid_search import GridSearchCV  # Perforing grid search
#from sklearn.svm import LinearSVC as SVC

import os

# set backend
os.environ['KERAS_BACKEND'] = 'tensorflow'
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"   # see issue #152
os.environ["CUDA_VISIBLE_DEVICES"] = ""



from keras import backend as K
# K.set_image_data_format('channels_first')

try:
    import theano.ifelse
    from theano.ifelse import IfElse, ifelse
except:
    pass

try:
    from tensorflow import set_random_seed
except:
    set_random_seed = np.random.seed


try:
    import lightgbm as lgb
except:
    pass

from cPickle import load, dump

# global parameters:
global_kmr = 3
global_len = 128
global_split_rate = 1./3
global_epoch = 32

#global_kmr = 5
#global_len = 256

###############################################################################
# multiple-gpu support
###############################################################################
'''
from keras.layers import merge
from keras.layers.core import Lambda
from keras.models import Model
import tensorflow as tf

def get_slice(data, idx, parts):
    shape = tf.shape(data)
    size = tf.concat([shape[:1] // parts, shape[1:]], axis=0)
    stride = tf.concat([shape[:1] // parts, shape[1:] * 0], axis=0)
    start = stride * idx
    return tf.slice(data, start, size)


def make_parallel(model, gpu_count):

    outputs_all = []
    for i in range(len(model.outputs)):
        outputs_all.append([])

    # Place a copy of the model on each GPU, each getting a slice of the batch
    for i in range(gpu_count):
        with tf.device('/gpu:%d' % i):
            with tf.name_scope('tower_%d' % i) as scope:

                inputs = []
                # Slice each input into a piece for processing on this GPU
                for x in model.inputs:
                    input_shape = tuple(x.get_shape().as_list())[1:]
                    slice_n = Lambda(get_slice, output_shape=input_shape, arguments={
                                     'idx': i, 'parts': gpu_count})(x)
                    inputs.append(slice_n)

                outputs = model(inputs)

                if not isinstance(outputs, list):
                    outputs = [outputs]

                # Save all the outputs for merging back together later
                for l in range(len(outputs)):
                    outputs_all[l].append(outputs[l])

    # merge outputs on CPU
    with tf.device('/cpu:0'):
        merged = []
        for outputs in outputs_all:
            merged.append(merge(outputs, mode='concat', concat_axis=0))

        return Model(input=model.inputs, output=merged)
'''


#os.environ['CUDA_VISIBLE_DEVICES'] = '-1'


def f1_score(y_true, y_pred):
    def recall(y_true, y_pred):
        """Recall metric.

        Only computes a batch-wise average of recall.

        Computes the recall, a metric for multi-label classification of
        how many relevant items are selected.
        """
        true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
        possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
        recall = true_positives / (possible_positives + K.epsilon())
        return recall

    def precision(y_true, y_pred):
        """Precision metric.

        Only computes a batch-wise average of precision.

        Computes the precision, a metric for multi-label classification of
        how many selected items are relevant.
        """
        true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
        predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
        precision = true_positives / (predicted_positives + K.epsilon())
        return precision
    precision = precision(y_true, y_pred)
    recall = recall(y_true, y_pred)
    return 2 * ((precision * recall) / (precision + recall))

fbeta_score = f1_score


##########################################################################
# overlap of two gene
##########################################################################
overlap = lambda s0, e0, s1, e1: min(e0, e1) - max(s0, s1) + 1


##########################################################################
# share kmer between 2 sequence
##########################################################################
def pearson(x, y):
    N, M = len(x), len(y)
    assert N == M
    x_m, y_m = sum(x) * 1. / N, sum(y) * 1. / M
    a, b, c = 0., 0., 0.
    for i in xrange(N):
        xi, yi = x[i] - x_m, y[i] - y_m
        a += xi * yi
        b += xi ** 2
        c += yi ** 2
    try:
        return a / sqrt(b * c)
    except:
        return 0


def sharekmer(s1, s2):
    # SynonymousCodons
    n1, n2 = map(len, [s1, s2])
    k1 = [s1[elem: elem + 3] for elem in xrange(0, n1, 3)]
    k2 = [s1[elem: elem + 3] for elem in xrange(0, n2, 3)]
    fq1 = Counter(k1)
    fq2 = Counter(k2)
    flag = 0

    kmers = []
    for i in SynonymousCodons:
        j = SynonymousCodons[i]
        if len(j) < 2:
            continue
        c1 = [[fq1[elem], elem] for elem in j]
        best1 = max(c1, key=lambda x: x[0])
        c2 = [[fq2[elem], elem] for elem in j]
        best2 = max(c2, key=lambda x: x[0])

        #c1.sort(key = lambda x: x[0], reverse = True)
        #c2.sort(key = lambda x: x[0], reverse = True)

        if best1[1] == best2[1]:
            kmers.append(best1[1])

    # for val in SynonymousCodons.values():
    #	if len(val) > 5:
    #		kmers.extend(val)

    # print 'the kmer', kmers, len(kmers)
    vec1 = [fq1[elem] for elem in kmers]
    vec2 = [fq2[elem] for elem in kmers]

    return pearson(vec1, vec2)


##########################################################################
# the motif found
##########################################################################
box_up10 = ['TATAAT', [77, 76, 60, 61, 56, 82]]
box_up35 = ['TTGACA', [69, 79, 61, 56, 54, 54]]
# find the best region that may be a candidate of a motif


def find_motif(seq, motif, bg=None):
    if bg is None:
        bg = {}
    l = len(motif[0])
    #best = float('-inf')
    best = -100
    idx = -1
    for i in xrange(0, len(seq) - l + 1):
        lmer = seq[i: i + l]
        score = 0
        for a, b, c in zip(lmer, motif[0], motif[1]):
            if a == b:
                score += log(float(c) / bg.get(a, 1.))
            else:
                score += log((100. - c) / bg.get(a, 1.))
                # try:
                #       score += log((100. - c) / bg.get(a, 1.))
                # except:
                #       print c, bg.get(a, 1.)

        if score >= best:
            idx = i
            best = score

    return [seq[idx: idx + l], len(seq) - idx, best]


##########################################################################
# cai, from biopython
##########################################################################
index = Counter({'GCT': 1, 'CGT': 1, 'AAC': 1, 'GAC': 1, 'TGC': 1, 'CAG': 1, 'GAA': 1, 'GGT': 1, 'CAC': 1, 'ATC': 1, 'CTG': 1, 'AAA': 1, 'ATG': 1, 'TTC': 1, 'CCG': 1, 'TCT': 1, 'ACC': 1, 'TGG': 1, 'TAC': 1, 'GTT': 1, 'ACT': 0.965, 'TCC': 0.744, 'GGC': 0.724, 'GCA': 0.586, 'TGT': 0.5, 'GTA': 0.495, 'GAT': 0.434, 'GCG': 0.424, 'AGC': 0.41, 'CGC': 0.356, 'TTT': 0.296, 'CAT': 0.291, 'GAG': 0.259,
                 'AAG': 0.253, 'TAT': 0.239, 'GTG': 0.221, 'ATT': 0.185, 'CCA': 0.135, 'CAA': 0.124, 'GCC': 0.122, 'ACG': 0.099, 'AGT': 0.085, 'TCA': 0.077, 'ACA': 0.076, 'CCT': 0.07, 'GTC': 0.066, 'AAT': 0.051, 'CTT': 0.042, 'CTC': 0.037, 'TTA': 0.02, 'TTG': 0.02, 'GGG': 0.019, 'TCG': 0.017, 'CCC': 0.012, 'GGA': 0.01, 'CTA': 0.007, 'AGA': 0.004, 'CGA': 0.004, 'CGG': 0.004, 'ATA': 0.003, 'AGG': 0.002})


def cai(seq):
    if seq.islower():
        seq = seq.upper()

    N = len(seq)
    cai_value, cai_length = 0, 0
    for i in xrange(0, N, 3):
        codon = seq[i: i + 3]
        if codon in index:
            if codon not in ['ATG', 'TGG']:
                cai_value += math.log(index[codon])
                cai_length += 1
        elif codon not in ['TGA', 'TAA', 'TAG']:
            continue
        else:
            continue

    if cai_length > 0:
        return math.exp(cai_value / cai_length)
    else:
        return 0


##########################################################################
# get the features
##########################################################################
# convert ATCG based kmer number
#code = {'A': 1, 'a': 1, 'T': 2, 't': 2, 'G': 3, 'g': 3, 'C': 4, 'c': 4}
code = [0] * 256
code5 = [0] * 256
flag = 0
for i in 'ATGC':
    code[ord(i.lower())] = code[ord(i)] = flag
    code5[ord(i.lower())] = code5[ord(i)] = flag + 1
    flag += 1

# convert string to number


def s2n(s, code=code, scale=None):
    if scale == None:
        scale = max(code) + 1
    N = 0
    output = 0
    for i in s[::-1]:
        #output += code.get(i, 0) * scale ** N
        output += code[ord(i)] * scale ** N
        N += 1

    return output

# reverse of s2n


def n2s(n, length, alpha='ATGC', scale=None):
    if scale == None:
        scale = max(code) + 1
    N = n
    s = []
    for i in xrange(length):
        s.append(alpha[N % scale])
        N /= scale

    return ''.join(s[::-1])


# convert the dna sequence to kmer-position matrix.
# if length of dna < given, then add NNN in the center of the sequence.
# else if length of dna > given, then trim the center of the sequence.

# the new kpm, reshape
def kpm(S, d=64, k=3, code=code, scale=None):
    if scale == None:
        scale = max(code) + 1

    N = scale ** k
    assert isinstance(d, int)
    L = len(S)
    if d < L:
        F = d // 2
        R = d - F
        seq = ''.join([S[: F], S[-R:]])
    elif d > L:
        F = L // 2
        R = L - F
        seq = ''.join([S[: F], 'N' * (d - L), S[-R:]])
    else:
        seq = S

    mat = [[0] * (d // k) for elem in xrange(N * k)]
    for i in xrange(0, d - k + 1):
        kmer = seq[i: i + k]
        if 'N' in kmer or 'n' in kmer:
            continue
        R = s2n(kmer, code=code, scale=scale)
        mat[R + i % k * N][i // k] = 1

    mat = np.asarray(mat, 'int8')
    return mat


def kpm0(S, d=64, k=3, code=code, scale=None):
    if scale == None:
        scale = max(code) + 1

    N = scale ** k
    assert isinstance(d, int)
    L = len(S)
    if d < L:
        F = d // 2
        R = d - F
        seq = ''.join([S[: F], S[-R:]])
    elif d > L:
        F = L // 2
        R = L - F
        seq = ''.join([S[: F], 'N' * (d - L), S[-R:]])
    else:
        seq = S

    mat = [[0] * (d // 3) for elem in xrange(N * 3)]
    for i in xrange(0, d - k + 1):
        kmer = seq[i: i + k]
        if 'N' in kmer or 'n' in kmer:
            continue
        R = s2n(kmer, code=code, scale=scale)
        mat[R + i % 3 * N][i // 3] = 1

    mat = np.asarray(mat, 'int8')
    return mat


# get features by give loc1, start and end:
# get xx
def get_xx(j, seq_dict, kmer=2, dim=128, mode='train', context=False):
    loc1, scf1, std1, st1, ed1, loc2, scf2, std2, st2, ed2 = j[: 10]
    if scf1 != scf2 or std1 != std2:
        if context:
            X0 = np.ones((4 ** kmer * kmer, dim // kmer * kmer))
        else:
            X0 = np.ones((4 ** kmer * kmer, dim // kmer))
        X1 = [10**4] * 11
        X2 = [127] * dim
        return [X0], X1, X2

    # get the sequence
    st1, ed1, st2, ed2 = map(int, [st1, ed1, st2, ed2])
    st1 -= 1
    st2 -= 1

    if st1 > st2:
        loc1, scf1, std1, st1, ed1, loc2, scf2, std2, st2, ed2 = loc2, scf2, std2, st2, ed2, loc1, scf1, std1, st1, ed1

    seq1 = seq_dict[scf1][st1: ed1]
    seq1 = std1 == '+' and seq1 or seq1.reverse_complement()
    seq2 = seq_dict[scf2][st2: ed2]
    seq2 = std1 == '+' and seq2 or seq2.reverse_complement()

    start, end = ed1, st2
    seq12 = seq_dict[scf1][start: end]

    seq12 = std1 == '+' and seq12 or seq12.reverse_complement()
    seq1, seq2, seq12 = map(str, [seq1.seq, seq2.seq, seq12.seq])
    seq1, seq2, seq12 = seq1.upper(), seq2.upper(), seq12.upper()

    # 1D features such as gc, dist
    cai1, cai2, cai12 = map(cai, [seq1, seq2, seq12])
    dist = st2 - ed1
    distn = (st2 - ed1) * 1. / (ed2 - st1)
    ratio = math.log((ed1 - st1) * 1. / (ed2 - st2))
    ratio = std1 == '+' and ratio or -ratio
    idx = -100
    bgs = Counter(seq12[idx:])
    up10, up35 = find_motif(seq12[idx:], box_up10, bgs), find_motif(
        seq12[idx:], box_up35, bgs)
    if seq12[idx:]:
        gc = SeqUtils.GC(seq12[idx:])
        try:
            skew = SeqUtils.GC_skew(seq12[idx:])[0]
        except:
            skew = 0.
    else:
        gc = skew = 0.

    bias = sharekmer(seq1, seq2)
    if st1 == st2 == '+':
        X1 = [cai1, cai2, bias, distn, ratio, gc, skew] + up10[1:] + up35[1:]
    else:
        X1 = [cai2, cai1, bias, distn, ratio, gc, skew] + up10[1:] + up35[1:]

    # 2D features of kmer matrix
    if context:
        seqmat12 = kpm(seq12, d=dim, k=kmer, scale=4)
        seqmat1 = kpm(seq1, d=dim, k=kmer, scale=4)
        seqmat2 = kpm(seq2, d=dim, k=kmer, scale=4)
        seqmat = np.concatenate((seqmat1, seqmat12, seqmat2), 1)
    else:
        seqmat = kpm(seq12, d=dim, k=kmer, scale=4)

    if ed1 > st2:
        seqmat[:] = 0
    X0 = [seqmat]
    n12 = len(seq12)
    X2 = [s2n(seq12[elem: elem + kmer], code5)
          for elem in xrange(n12 - kmer + 1)]

    return X0, X1, X2


def get_xx0(j, seq_dict, kmer=2, dim=128, mode='train', context=False):
    loc1, scf1, std1, st1, ed1, loc2, scf2, std2, st2, ed2 = j[: 10]
    if scf1 != scf2 or std1 != std2:
        if context:
            X0 = np.ones((4 ** kmer * 3, dim // 3 * 3))
        else:
            X0 = np.ones((4 ** kmer * 3, dim // 3))
        X1 = [10**4] * 11
        X2 = [127] * dim
        return [X0], X1, X2

    # get the sequence
    st1, ed1, st2, ed2 = map(int, [st1, ed1, st2, ed2])
    st1 -= 1
    st2 -= 1

    if st1 > st2:
        loc1, scf1, std1, st1, ed1, loc2, scf2, std2, st2, ed2 = loc2, scf2, std2, st2, ed2, loc1, scf1, std1, st1, ed1

    seq1 = seq_dict[scf1][st1: ed1]
    seq1 = std1 == '+' and seq1 or seq1.reverse_complement()
    seq2 = seq_dict[scf2][st2: ed2]
    seq2 = std1 == '+' and seq2 or seq2.reverse_complement()

    start, end = ed1, st2
    seq12 = seq_dict[scf1][start: end]

    seq12 = std1 == '+' and seq12 or seq12.reverse_complement()
    seq1, seq2, seq12 = map(str, [seq1.seq, seq2.seq, seq12.seq])
    seq1, seq2, seq12 = seq1.upper(), seq2.upper(), seq12.upper()

    # 1D features such as gc, dist
    cai1, cai2, cai12 = map(cai, [seq1, seq2, seq12])
    dist = st2 - ed1
    distn = (st2 - ed1) * 1. / (ed2 - st1)
    ratio = math.log((ed1 - st1) * 1. / (ed2 - st2))
    ratio = std1 == '+' and ratio or -ratio
    idx = -100
    bgs = Counter(seq12[idx:])
    up10, up35 = find_motif(seq12[idx:], box_up10, bgs), find_motif(
        seq12[idx:], box_up35, bgs)
    if seq12[idx:]:
        gc = SeqUtils.GC(seq12[idx:])
        try:
            skew = SeqUtils.GC_skew(seq12[idx:])[0]
        except:
            skew = 0.
    else:
        gc = skew = 0.

    bias = sharekmer(seq1, seq2)
    if st1 == st2 == '+':
        X1 = [cai1, cai2, bias, distn, ratio, gc, skew] + up10[1:] + up35[1:]
    else:
        X1 = [cai2, cai1, bias, distn, ratio, gc, skew] + up10[1:] + up35[1:]

    # 2D features of kmer matrix
    if context:
        seqmat12 = kpm(seq12, d=dim, k=kmer, scale=4)
        seqmat1 = kpm(seq1, d=dim, k=kmer, scale=4)
        seqmat2 = kpm(seq2, d=dim, k=kmer, scale=4)
        seqmat = np.concatenate((seqmat1, seqmat12, seqmat2), 1)
    else:
        seqmat = kpm(seq12, d=dim, k=kmer, scale=4)

    if ed1 > st2:
        seqmat[:] = 0
    X0 = [seqmat]
    n12 = len(seq12)
    X2 = [s2n(seq12[elem: elem + kmer], code5)
          for elem in xrange(n12 - kmer + 1)]

    return X0, X1, X2


# get single line of features
def get_xx_one(j, seq_dict, kmer=2, dim=128, mode='train'):
    X0, X1, X2 = get_xx(j, seq_dict, kmer, dim, mode)
    x0, x1, x2 = map(np.asarray, [[X0], [X1], [X2]])
    return x0, x1, X2

# generate training and testing data


def get_xxy(f, seq_dict, kmer=2, dim=128):
    # get the training data
    X0, X1, X2, y = [], [], [], []

    for i in f:
        j = i[:-1].split('\t')
        x0, x1, x2 = get_xx(j, seq_dict, kmer, dim)
        X0.append(x0)
        X1.append(x1)
        X2.append(x2)
        y.append(j[-1] == 'True' and 1 or 0)

    X0 = np.asarray(X0, 'int8')
    X1 = np.asarray(X1, 'float32')
    X2 = np.asarray(X2)
    y = np.asarray(y, 'int8')
    return X0, X1, X2, y

# split the X0, X1, y data to training and testing


def split_xxy(X0, X1, X2, y, train_size=1. / 3, seed=42):
    N = X0.shape[0]
    idx = np.arange(N)
    np.random.seed(seed)
    np.random.shuffle(idx)
    start = int(train_size * N)
    idx_train, idx_test = idx[: start], idx[start:]
    X0_train, X1_train, X2_train, y_train = X0[
        idx_train], X1[idx_train], X2[idx_train], y[idx_train]
    X0_test, X1_test, X2_test, y_test = X0[idx_test], X1[
        idx_test], X2[idx_test], y[idx_test]

    return X0_train, X1_train, X2_train, y_train, X0_test, X1_test, X2_test, y_test




##########################################################################
# the lstm based
##########################################################################

# get lstm features by give loc1, start and end:
def get_lstm_xx(j, seq_dict, kmer=2, dim=128, mode='train'):
    loc1, scf1, std1, st1, ed1, loc2, scf2, std2, st2, ed2 = j[: 10]
    if scf1 != scf2 or std1 != std2:
        #X0 = np.ones((4 ** kmer, dim))
        X0 = [127] * dim
        #X0 = None
        X1 = [10**4] * 11
        return X0, X1

    # get the sequence
    st1, ed1, st2, ed2 = map(int, [st1, ed1, st2, ed2])
    st1 -= 1
    st2 -= 1

    if st1 > st2:
        loc1, scf1, std1, st1, ed1, loc2, scf2, std2, st2, ed2 = loc2, scf2, std2, st2, ed2, loc1, scf1, std1, st1, ed1
    seq1 = seq_dict[scf1][st1: ed1]
    seq1 = std1 == '+' and seq1 or seq1.reverse_complement()
    seq2 = seq_dict[scf2][st2: ed2]
    seq2 = std1 == '+' and seq2 or seq2.reverse_complement()

    start, end = ed1, st2
    seq12 = seq_dict[scf1][start: end]

    # if len(seq12) > dim:
    #    seq12 = seq12[: dim // 2] + seq12[-dim // 2: ]

    seq12 = std1 == '+' and seq12 or seq12.reverse_complement()
    seq1, seq2, seq12 = map(str, [seq1.seq, seq2.seq, seq12.seq])
    seq1, seq2, seq12 = seq1.upper(), seq2.upper(), seq12.upper()

    # 1D features such as gc, dist
    cai1, cai2, cai12 = map(cai, [seq1, seq2, seq12])
    dist = st2 - ed1
    distn = (st2 - ed1) * 1. / (ed2 - st1)
    ratio = math.log((ed1 - st1) * 1. / (ed2 - st2))
    ratio = std1 == '+' and ratio or -ratio
    idx = -100
    bgs = Counter(seq12[idx:])
    up10, up35 = find_motif(seq12[idx:], box_up10, bgs), find_motif(
        seq12[idx:], box_up35, bgs)
    if seq12[idx:]:
        gc = SeqUtils.GC(seq12[idx:])
        try:
            skew = SeqUtils.GC_skew(seq12[idx:])[0]
        except:
            skew = 0.
    else:
        gc = skew = 0.

    bias = sharekmer(seq1, seq2)
    if st1 == st2 == '+':
        X1 = [cai1, cai2, bias, distn, ratio, gc, skew] + up10[1:] + up35[1:]
    else:
        X1 = [cai2, cai1, bias, distn, ratio, gc, skew] + up10[1:] + up35[1:]
        #X1 = [cai1, cai2, bias, distn, ratio, gc, skew] + up10[1: ] + up35[1: ]

    # 1D features of lstm
    n12 = len(seq12)

    '''
    L = dim // 2
    R = dim - L
    if n12 > dim:
        seq12 = seq12[: L] + seq12[-R: ]
    else:
        seq12 = seq12[: L] +'N' * (dim - n12)  + seq12[-R: ]
    '''
    lstm_seq = [s2n(seq12[elem: elem + kmer], code5)
                for elem in xrange(n12 - kmer + 1)]
                
    #X0 = lstm_seq[::kmer]
    #for i in xrange(kmer):
    #    X0.extend(lstm_seq[i::kmer])
    #X0 = lstm_seq
    X0 = [-1] * dim
    ndim = len(lstm_seq)
    if ndim == dim:
        X0 = lstm_seq
        #print 'X0 0', len(X0)

    elif 2 <= ndim < dim:
        ndim = ndim // 2
        X0[:ndim] = lstm_seq[:ndim]
        X0[-ndim:] = lstm_seq[-ndim:]
        #print 'X0 1', len(X0)
    elif ndim > dim:
        ndim = dim // 2
        X0[:ndim] = lstm_seq[:ndim]
        X0[-ndim:] = lstm_seq[-ndim:]
        #print 'X0 2', len(X0)
    else:
        pass

    return X0, X1

# get single line of lstm features


def get_lstm_xx_one(j, seq_dict, kmer=2, dim=128, mode='train'):
    X0, X1 = get_lstm_xx(j, seq_dict, kmer, dim, mode)
    x0, x1 = map(np.asarray, [[X0], [X1]])
    return x0, x1

# generate training and testing data


def get_lstm_xxy(f, seq_dict, kmer=2, dim=128):
    # get the training data
    X0, X1, y = [], [], []
    for i in f:
        j = i[:-1].split('\t')
        #print 'line', j
        x0, x1 = get_lstm_xx(j, seq_dict, kmer, dim)
        #print 'X0 is', len(x0)
        # if len(x0) < 1 or x0[0] == -1:
        #    continue
        cat = j[-1]
        X0.append(x0)
        X1.append(x1)
        y.append(cat == 'True' and 1 or 0)

    #print max(map(len, X0)), min(map(len, X0)), X0[: 2]

    X0 = np.asarray(X0)
    print 'X0 shape', X0.shape
    #print 'cat1, cat2, cai12, gc, skew, dist, distn, ratio, up10, up35'
    #print 'X1', X1[0]
    X1 = np.asarray(X1, 'float32')
    y = np.asarray(y, 'int8')

    #print 'X0', X0, X1, y
    return X0, X1, y


# split lstm xxy
def split_lstm_xxy(X0, X1, y, train_size=1. / 3, seed=42):
    #X = np.asarray(X)
    #y = np.asarray(y)
    N = X0.shape[0]
    idx = np.arange(N)
    np.random.seed(seed)
    np.random.shuffle(idx)
    start = int(train_size * N)
    idx_train, idx_test = idx[: start], idx[start:]
    X0_train, X1_train, y_train = X0[idx_train], X1[idx_train], y[idx_train]
    X0_test, X1_test, y_test = X0[idx_test], X1[idx_test], y[idx_test]

    return X0_train, X1_train, y_train, X0_test, X1_test, y_test



##########################################################################
# the CNN class
##########################################################################
class CNN:

    def __init__(self, nb_filter=32, nb_pool=3, nb_conv=2, nb_epoch=10, batch_size=64, maxlen=128, save_path='./weights.hdf5'):
        # def __init__(self, nb_filter=64, nb_pool=3, nb_conv=2, nb_epoch=10,
        # batch_size=32, maxlen=128, save_path='./weights.hdf5'):

        self.nb_filter = nb_filter
        self.nb_pool = nb_pool
        self.nb_conv = nb_conv
        self.nb_epoch = nb_epoch
        self.batch_size = batch_size
        self.maxlen = maxlen
        self.opt = Adam(lr=5e-4, beta_1=0.995, beta_2=0.999, epsilon=1e-09)
        #self.checkpointer = [ModelCheckpoint(filepath=save_path, verbose=1, save_best_only=True, mode='max', monitor='val_fbeta_score')]
        self.checkpointer = [ModelCheckpoint(
            filepath=save_path, verbose=1, save_best_only=True, mode='max', monitor='val_f1_score')]
        #self.metric = keras.metrics.fbeta_score
        self.metric = f1_score
        self.cross_val = 1. / 3

    def fit_2d(self, X_train, y_train, X_test=None, y_test=None):
        np.random.seed(42)
        set_random_seed(42)
        Y_train = np_utils.to_categorical(y_train).astype('float32')
        if type(y_test) == type(None):
            Y_test = None
        else:
            Y_test = np_utils.to_categorical(y_test).astype('float32')

        nb_classes = Y_train.shape[1]

        # set parameter for cnn
        loss = nb_classes > 2 and 'categorical_crossentropy' or 'binary_crossentropy'
        print 'loss function is', loss
        # number of convolutional filters to use
        nb_filters = self.nb_filter
        # size of pooling area for max pooling
        nb_pool = self.nb_pool
        # convolution kernel size
        nb_conv = self.nb_conv
        # traning iteration
        nb_epoch = self.nb_epoch
        batch_size = self.batch_size
        a, b, img_rows, img_cols = X_train.shape

        # set the conv model
        model = Sequential()
        #model.add(Convolution2D(nb_filters, (nb_conv, nb_conv), border_mode='same', input_shape=(b, img_rows, img_cols), activation='relu', name='conv1_1', data_format='channels_first'))
        model.add(Conv2D(nb_filters, (nb_conv, nb_conv), border_mode='same', input_shape=(b, img_rows, img_cols), activation='relu', name='conv1_1', data_format='channels_first'))
        #model.add(Conv2D(nb_filters, 1, border_mode='same', input_shape=(b, img_rows, img_cols), name='conv1_1', data_format='channels_first'))

        #model.add(Conv2D(64, (3, 3), activation='relu', name='conv1_2', data_format='channels_first'))
        #model.add(Conv2D(64, (3, 3), activation='relu', name='conv1_2', data_format='channels_first'))
        model.add(MaxPooling2D((3, 3), strides=(2, 2)))

        #model.add(Conv2D(128, (3, 3), activation='relu', name='conv2_1', data_format='channels_first'))
        #model.add(MaxPooling2D((2, 2), strides=(2, 2)))

        #model.add(Conv2D(256, (3, 3), activation='relu', name='conv3_1', data_format='channels_first'))
        #model.add(MaxPooling2D((2, 2), strides=(2, 2)))

        #model.add(Conv2D(512, (3, 3), activation='relu', name='conv4_1', data_format='channels_first'))

        model.add(Flatten())

        model.add(Dense(128, activation='relu'))
        #model.add(Dense(512, activation='relu'))
        model.add(Dropout(0.85))
        model.add(Dense(nb_classes, activation='sigmoid'))

        # try:
        #    model = make_parallel(model, 2)
        # except:
        #    pass

        opt = self.opt
        model.compile(loss=loss, optimizer='adam', metrics=[self.metric])

        model.summary()

        # set the check pointer to save the best model
        if type(X_test) != type(None) and type(Y_test) != type(None):
        #if 0:
            model.fit(X_train, Y_train, batch_size=batch_size, epochs=nb_epoch, verbose=1,
                      validation_data=(X_test, Y_test), shuffle=True, validation_split=1e-4, callbacks=self.checkpointer)
        else:
            model.fit(X_train, Y_train, batch_size=batch_size,
                      epochs=nb_epoch, verbose=1, shuffle=True, validation_split=self.cross_val, callbacks=self.checkpointer)

        self.model_2d = model

    def retrain_2d(self, X_train, y_train, X_test=None, y_test=None):

        Y_train = np_utils.to_categorical(y_train).astype('float32')
        if type(y_test) == type(None):
            Y_test = None
        else:
            Y_test = np_utils.to_categorical(y_test).astype('float32')

        if type(X_test) != type(None) and type(Y_test) != type(None):
        #if 0:
            self.model_2d.fit(X_train, Y_train, batch_size=self.batch_size, epochs=self.nb_epoch, verbose=1,
                      validation_data=(X_test, Y_test), shuffle=True, validation_split=1e-4, callbacks=self.checkpointer)
        else:
            self.model_2d.fit(X_train, Y_train, batch_size=self.batch_size,
                      epochs=self.nb_epoch, verbose=1, shuffle=True, validation_split=self.cross_val, callbacks=self.checkpointer)

    def predict_2d(self, X):
        return self.model_2d.predict(X).argmax(1)

    def fit_lstm(self, X_train, y_train, X_test=None, y_test=None):
        self.max_features = 2**12
        #print X_train.shape, y_train.shape
        N, D = X_train.shape
        model = Sequential()
        model.add(Embedding(self.max_features, D))
        #model.add(LSTM(D, dropout=0.2, recurrent_dropout=0.2))
        model.add(Bidirectional(CuDNNGRU(D, return_sequences=True)))
        #model.add(Bidirectional(GRU(D, return_sequences=True)))
        #model.add(BatchNormalization(axis=-1, momentum=0.99, epsilon=0.001, center=True, scale=True, beta_initializer='zeros', gamma_initializer='ones', moving_mean_initializer='zeros', moving_variance_initializer='ones', beta_regularizer=None, gamma_regularizer=None, beta_constraint=None, gamma_constraint=None))
        model.add(Dropout(0.2))

        model.add(Bidirectional(CuDNNGRU(D)))
        #model.add(Bidirectional(GRU(D)))

        #model.add(BatchNormalization(axis=-1, momentum=0.99, epsilon=0.001, center=True, scale=True, beta_initializer='zeros', gamma_initializer='ones', moving_mean_initializer='zeros', moving_variance_initializer='ones', beta_regularizer=None, gamma_regularizer=None, beta_constraint=None, gamma_constraint=None))
        model.add(Dropout(0.2))

        model.add(Dense(1, activation='sigmoid'))

        # try using different optimizers and different optimizer configs
        nb_classes = len(set(y_train))
        loss = nb_classes > 2 and 'categorical_crossentropy' or 'binary_crossentropy'
        #model.compile(loss='binary_crossentropy', optimizer='adam', metrics=[self.metric])
        model.compile(loss=loss, optimizer='adam', metrics=[self.metric])
        print('Train..., loss is %s %s'%(loss, D))
        if type(X_test) != type(None) and type(y_test) != type(None):
            model.fit(X_train, y_train, batch_size=self.batch_size, epochs=self.nb_epoch, validation_data=(X_test, y_test), shuffle=True, callbacks=self.checkpointer)
        else:
            model.fit(X_train, y_train, batch_size=self.batch_size, epochs=self.nb_epoch, validation_data=(X_test, y_test), verbose=1, shuffle=True, validation_split=self.cross_val, callbacks=self.checkpointer)
        score, acc = model.evaluate(X_test, y_test, batch_size=self.batch_size)
        print('Test score:', score)
        print('Test accuracy:', acc)
        self.model_2d = model


    def retrain_lstm(self, X_train, y_train, X_test=None, y_test=None):
        if type(X_test) != type(None) and type(y_test) != type(None):
            self.model_2d.fit(X_train, y_train, batch_size=self.batch_size, epochs=self.nb_epoch, validation_data=(X_test, y_test), shuffle=True, callbacks=self.checkpointer)
        else:
            self.model_2d.fit(X_train, y_train, batch_size=self.batch_size, epochs=self.nb_epoch, validation_data=(X_test, y_test), verbose=1, shuffle=True, validation_split=self.cross_val, callbacks=self.checkpointer)
        score, acc = self.model_2d.evaluate(X_test, y_test, batch_size=self.batch_size)




    def predict_lstm(self, X):
        return self.model_2d.predict_classes(X).flatten()

    # load an training model
    def load(self, name, mode='2d'):
        #model = keras.models.load_model(name, custom_objects={'f1_score': f1_score, 'tf': tf, 'fbeta_score': f1_score})
        model = keras.models.load_model(name, custom_objects={'f1_score': f1_score, 'fbeta_score': f1_score})
        if mode == '2d' or mode == 'lstm':
            self.model_2d = model
        else:
            pass

    # save the model
    def save(self, name, model='2d'):
        if model == '2d' or model == 'lstm':
            self.model_2d.save(name + '_' + model + '_%d_%d.hdf5' % (global_kmr, global_len))
        else:
            pass


# run training
def run_train(train, seq_dict, clf, mode='2d', rounds=0):
    # get the training data
    split_rate = global_split_rate

    if mode == '2d':
        print 'CNN training'
        f = open(train, 'r')
        #X, X1, X2, y = get_xxy(f, seq_dict, 3, 128)
        #X, X1, X2, y = get_xxy(f, seq_dict, 4, 64)
        X, X1, X2, y = get_xxy(f, seq_dict, global_kmr, global_len)
        X_train, X1_train, X2_train, y_train, X_test, X1_test, X2_test, y_test = split_xxy(
            X, X1, X2, y, split_rate)
        f.close()

        #X_train = np.asarray(X_train, 'float32')
        #X_test = np.asarray(X_test, 'float32')

        #clf.fit_2d(X_train, y_train, X_test, y_test)
        clf.fit_2d(X_train, y_train)

        for idx in xrange(rounds):
            print 'retraining cnn', idx
            #X_train, X1_train, X2_train, y_train, X_test, X1_test, X2_test, y_test = split_xxy(X, X1, X2, y, split_rate, seed=43+idx)
            #clf.retrain_2d(X_train, y_train, X_test, y_test)
            X_train0, X1_train0, X2_train0, y_train0, X_test0, X1_test0, X2_test0, y_test0 = split_xxy(X_train, X1_train, X2_train, y_train, 3./4, seed=43+idx)
            clf.retrain_2d(X_train0, y_train0, X_test0, y_test0)
            #clf.retrain_2d(X_train0, y_train0, X_test, y_test)
        # the test score
        Y_test = np_utils.to_categorical(y_test).astype('float32')
        score = clf.model_2d.evaluate(X_test, Y_test, verbose=0)

        print('   Test score:', score[0])
        print('Test accuracy:', score[1])

        # validate
        y_test_pred = clf.predict_2d(X_test)
        y_train_pred = clf.predict_2d(X_train)
    elif mode == 'lstm':

        f = open(train, 'r')
        #X, X1, X2, y = get_xxy(f, seq_dict, 3, 128)
        #X, X1, X2, y = get_xxy(f, seq_dict, 4, 64)
        X, X1, y = get_lstm_xxy(f, seq_dict, global_kmr, global_len)
        X_train, X1_train, y_train, X_test, X1_test, y_test = split_lstm_xxy(X, X1, y, split_rate)

        f.close()

        print 'LSTM training'
        print 'shape', X.shape, X_train.shape
        clf.fit_lstm(X_train, y_train, X_test, y_test)
        for idx in xrange(rounds):
            print 'retraining lstm', idx
            #X_train, X1_train, y_train, X_test, X1_test, y_test = split_lstm_xxy(X, X1, y, split_rate, seed=43+idx)
            #clf.retrain_lstm(X_train, y_train, X_test, y_test)

            X_train0, X1_train0, y_train0, X_test0, X1_test0, y_test0 = split_lstm_xxy(X_train, X1_train, y_train, 3./4, seed=43+idx)
            #clf.retrain_lstm(X_train0, y_train0, X_test, y_test)
            clf.retrain_lstm(X_train0, y_train0, X_test0, y_test0)

        # the test score
        #Y_test = np_utils.to_categorical(y_test).astype('float32')
        #score = clf.model_2d.evaluate(X_test, Y_test, verbose=0)
        #score = clf.model_2d.evaluate(X_test, y_test, verbose=0)
        y_test_pred = clf.predict_lstm(X_test)
        y_train_pred = clf.predict_lstm(X_train)
        print 'test', y_test_pred, y_test


    clf.save(train, mode)
    precise = metrics.precision_score(y_test, y_test_pred)
    recall = metrics.recall_score(y_test, y_test_pred)
    f1 = metrics.f1_score(y_test, y_test_pred)
    print 'test Precise:', precise
    print 'test  Recall:', recall
    print 'test      F1:', f1

    precise = metrics.precision_score(y_train, y_train_pred)
    recall = metrics.recall_score(y_train, y_train_pred)
    f1 = metrics.f1_score(y_train, y_train_pred)
    print 'train Precise:', precise
    print 'train Recall:', recall
    print 'train     F1:', f1



# train by lightgdm
def run_train_lgb(train, seq_dict, clf):
    # get the training data
    split_rate = 1. / 3

    f = open(train, 'r')
    #X, X1, X2, y = get_xxy(f, seq_dict, 3, 128)
    #X, X1, X2, y = get_xxy(f, seq_dict, 4, 64)
    X, X1, X2, y = get_xxy(f, seq_dict, global_kmr, global_len)
    X_train, X1_train, X2_train, y_train, X_test, X1_test, X2_test, y_test = split_xxy(
        X, X1, X2, y, split_rate)
    f.close()

    print 'data shape 1', X_train.shape, y_train.shape
    a0, a1, a2, a3 = X_train.shape
    X_train = X_train.reshape(a0, a1 * a2 * a3)
    print 'data shape 2', X_train.shape, y_train.shape

    clf.fit(X_train, y_train)
    # clf.save(train+'.lgb')

    _o = open(train + '.lgb', 'wb')
    dump(clf, _o)
    _o.close()

    a0, a1, a2, a3 = X_test.shape
    X_test = X_test.reshape(a0, a1 * a2 * a3)
    y_test_pred = clf.predict(X_test)

    precise = metrics.precision_score(y_test, y_test_pred)
    recall = metrics.recall_score(y_test, y_test_pred)
    f1 = metrics.f1_score(y_test, y_test_pred)
    print 'Precise:', precise
    print ' Recall:', recall
    print '     F1:', f1


# run the adjacent prediction
def run_adjacent_predict(adjacent, seq_dict, model, clf, mode='2d'):
    adjacent, model = sys.argv[3: 5]
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
    #print 'try loading lstm model', model
    clf.load(model, mode)

    # get the locus of genes
    f = open(adjacent, 'r')
    for i in f:
        j = i[:-1].split('\t')
        #x0, x1, x2 = get_xx_one(j, seq_dict, 3, 128, 'test')
        #x0, x1, x2 = get_xx_one(j, seq_dict, 4, 64, 'test')
        #x0, x1, x2 = get_xx_one(j, seq_dict, global_kmr, global_len, 'test')
        # print 'data shape', x0.shape, x1.shape
        if mode == '2d':
            x0, x1, x2 = get_xx_one(j, seq_dict, global_kmr, global_len, 'test')
            res = clf.predict_2d(x0)[0]
        else:
            x0, x1 = get_lstm_xx_one(j, seq_dict, global_kmr, global_len, 'test')
            res = clf.predict_lstm(x0)[0]

        res = res == 1 and 'True' or 'False'
        print i[: -1] + '\t' + str(res)

    f.close()


# run the adjacent prediction by lightgdm
def run_adjacent_predict_lgb(adjacent, seq_dict, model, clf, mode='2d'):
    adjacent, model = sys.argv[3: 5]
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
    #clf.load(model, mode)
    with open(model, 'rb') as md:
        clf = load(md)

    # get the locus of genes
    f = open(adjacent, 'r')
    for i in f:
        j = i[:-1].split('\t')
        #x0, x1, x2 = get_xx_one(j, seq_dict, 3, 128, 'test')
        #x0, x1, x2 = get_xx_one(j, seq_dict, 4, 64, 'test')
        x0, x1, x2 = get_xx_one(j, seq_dict, global_kmr, global_len, 'test')
        a0, a1, a2, a3 = x0.shape
        x0 = x0.reshape(a0, a1 * a2 * a3)
        # print 'data shape', x0.shape, x1.shape
        res = clf.predict(x0)
        res = res == 1 and 'True' or 'False'
        print i[: -1] + '\t' + str(res)

    f.close()


# run the whole genome prediction

# generate adjacent gene pairs from the gene list
def adjacent_genes(f):
    locus_list = []
    for i in f:
        j = i[: -1].split('\t')
        if len(j) < 7:
            j.extend([0] * 7)
        locus, scaf, strand, start, end = j[: 5]
        start, end = map(int, [start, end])
        locus_list.append([locus, scaf, strand, start, end])

    locus_list.sort(key=lambda x: x[1: 5])
    return locus_list


def run_genome_predict(genome, seq_dict, model, clf, mode='2d'):
    genome, model = sys.argv[3: 5]
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
    clf.load(model, mode)

    # get the locus of genes
    f = open(genome, 'r')
    locus_list = adjacent_genes(f)
    f.close()
    for a, b in zip(locus_list[: -1], locus_list[1:]):
        j = a + b
        #x0, x1, x2 = get_xx_one(j, seq_dict, 3, 128, 'test')
        x0, x1, x2 = get_xx_one(j, seq_dict, global_kmr, global_len, 'test')
        if mode == '2d':
            if a[1] == b[1] and a[2] == b[2]:
                res = clf.predict_2d(x0)[0]
            else:
                res = 0
        else:
            pass

        res = res == 1 and 'True' or 'False'
        i = '\t'.join(map(str, j))
        print i + '\t' + str(res)


if __name__ == '__main__':
    import sys
    if len(sys.argv[1:]) < 3:
        print '#' * 79
        print '# To train a model:'
        print '#' * 79
        print 'python this.py train foo.fasta foo.train.txt [mode]\n'
        print 'foo.train.txt is the gene location in the format:'
        print '       locus1\tscf1\tstrand1\tstart1\tend1\tlocus2\tscf2\tstrand2\tstart2\tend2\tcat\n'

        print '#' * 79
        print '# To make a adjacent genes prediction'
        print '#' * 79
        print 'python this.py adjacent foo.fasta foo.adjacent.txt foo.model [mode]\n'
        print 'foo.adjacent.txt is the gene location in the format:'
        print '       locus1\tscf1\tstrand1\tstart1\tend1\tlocus2\tscf2\tstrand2\tstart2\tend2\n'

        print '#' * 79
        print '# To make a whole genome prediction'
        print '#' * 79
        print 'python this.py genome foo.fasta foo.genome.txt foo.model [mode]'
        print 'foo.genome.txt is the gene location in the format:'
        print '       locus1\tscf1\tstrand1\tstart1\tend1'

        print ''
        print '#' * 79
        print 'start1/2: start of the gene in the genome, start > 0 need be adjust in the program'
        print '     cat: indicate whether  operon or not'
        print '    mode: 2d'
        raise SystemExit()

    import keras
    from keras.models import Sequential
    from keras.preprocessing import sequence
    from keras.layers import Dense, Dropout, Activation, Flatten, Embedding, BatchNormalization
    #from keras.layers import Input, Merge, LSTM, GRU, Bidirectional, UpSampling2D, InputLayer, CuDNNGRU
    from keras.layers import Input, LSTM, GRU, Bidirectional, UpSampling2D, InputLayer
    from keras.optimizers import SGD, Adam, RMSprop
    from keras.layers.convolutional import Convolution2D, MaxPooling2D, ZeroPadding2D, Conv2D
    from keras.utils import np_utils
    from keras.callbacks import ModelCheckpoint, TensorBoard
    from keras.models import Model
    from keras import backend as K
    from keras import objectives
    from keras.layers import Input, Dense, Lambda
    import numpy as np

    model, fasta = sys.argv[1: 3]

    # save the genome to an dict
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))

    if model.startswith('train'):
        train = sys.argv[3]
        try:
            mode = sys.argv[4]
        except:
            mode = '2d'

        clf = CNN(nb_epoch=global_epoch, maxlen=128, save_path=train + '_' +
                  mode + '_%d_%d' % (global_kmr, global_len) + '.hdf5')
        #clf = lgb.LGBMClassifier(n_estimators=50, max_bin=350, boosting_type='dart', learning_rate=.05)
        #clf = lgb.LGBMClassifier(n_estimators=150, learning_rate=.15)
        #clf = SVC(C=2)
        run_train(train, seq_dict, clf, mode)
        #run_train(train, seq_dict, clf, 'lstm')
        #run_train_lgb(train, seq_dict, clf)

    elif model.startswith('predict'):
        if len(sys.argv[1:]) < 4:
            print '#' * 79
            print '# To make a adjacent genes prediction'
            print '#' * 79
            print 'python this.py predict foo.fasta foo.adjacent.txt foo.model [mode]\n'
            print 'foo.adjacent.txt is the gene location in the format:'
            print '       locus1\tscf1\tstrand1\tstart1\tend1\tlocus2\tscf2\tstrand2\tstart2\tend2\n'
            raise SystemExit()

        test, model = sys.argv[3: 5]
        try:
            mode = sys.argv[5]
        except:
            mode = '2d'

        clf = CNN(nb_epoch=32, maxlen=128)

        # determine the number of col
        f = open(test, 'r')
        if sys.version_info[0] == 2:
            header = f.next().split('\t')
        else:
            header = next(f).split('\t')
        f.close()
        if header.count('+') + header.count('-') > 1:
            #print 'try loading lstm or 2d model'
            run_adjacent_predict(test, seq_dict, model, clf, mode)
            #run_adjacent_predict(test, seq_dict, model, clf, 'lstm')
            #run_adjacent_predict_lgb(test, seq_dict, model, clf, mode)

        else:
            run_genome_predict(test, seq_dict, model, clf, mode)
    else:
        pass
