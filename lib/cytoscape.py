#!usr/bin/env python
# the script to convert the cog PPI result to cytoscape format
# this script will generate 3 file:
# 1 is node node interaction
# 2 is node features
# 3 is edge features
import sys
import os
import random
from itertools import cycle

if len(sys.argv[1:]) < 1:
    #print('python this.py cog.txt cognames2003-2014.tab')
    print('python this.py cog.txt')
    raise SystemExit()

scriptpath='/'.join(os.path.realpath(__file__).split('/')[: -1]) + '/'


# colors
colors = ['AliceBlue', 'AntiqueWhite', 'Aqua', 'Aquamarine', 'Azure', 'Beige', 'Bisque', 'BlanchedAlmond', 'Blue', 'BlueViolet', 'Brown', 'BurlyWood', 'CadetBlue', 'Chartreuse', 'Chocolate', 'Coral', 'CornflowerBlue', 'Cornsilk', 'Crimson', 'Cyan', 'DarkBlue', 'DarkCyan', 'DarkGoldenRod', 'DarkGray', 'DarkGreen', 'DarkKhaki', 'DarkMagenta', 'DarkOliveGreen', 'DarkOrange', 'DarkOrchid', 'DarkRed', 'DarkSalmon', 'DarkSeaGreen', 'DarkSlateBlue', 'DarkSlateGray', 'DarkTurquoise', 'DarkViolet', 'DeepPink', 'DeepSkyBlue', 'DimGray', 'DodgerBlue', 'FireBrick', 'FloralWhite', 'ForestGreen', 'Fuchsia', 'Gainsboro', 'GhostWhite', 'Gold', 'GoldenRod', 'Gray', 'Green', 'GreenYellow', 'HoneyDew', 'HotPink', 'IndianRed', '', 'Indigo', '', 'Ivory', 'Khaki', 'Lavender', 'LavenderBlush', 'LawnGreen', 'LemonChiffon', 'LightBlue', 'LightCoral', 'LightCyan', 'LightGoldenRodYellow', 'LightGray', 'LightGreen', 'LightPink',
          'LightSalmon', 'LightSeaGreen', 'LightSkyBlue', 'LightSlateGray', 'LightSteelBlue', 'LightYellow', 'Lime', 'LimeGreen', 'Linen', 'Magenta', 'Maroon', 'MediumAquaMarine', 'MediumBlue', 'MediumOrchid', 'MediumPurple', 'MediumSeaGreen', 'MediumSlateBlue', 'MediumSpringGreen', 'MediumTurquoise', 'MediumVioletRed', 'MidnightBlue', 'MintCream', 'MistyRose', 'Moccasin', 'NavajoWhite', 'Navy', 'OldLace', 'Olive', 'OliveDrab', 'Orange', 'OrangeRed', 'Orchid', 'PaleGoldenRod', 'PaleGreen', 'PaleTurquoise', 'PaleVioletRed', 'PapayaWhip', 'PeachPuff', 'Peru', 'Pink', 'Plum', 'PowderBlue', 'Purple', 'RebeccaPurple', 'Red', 'RosyBrown', 'RoyalBlue', 'SaddleBrown', 'Salmon', 'SandyBrown', 'SeaGreen', 'SeaShell', 'Sienna', 'Silver', 'SkyBlue', 'SlateBlue', 'SlateGray', 'Snow', 'SpringGreen', 'SteelBlue', 'Tan', 'Teal', 'Thistle', 'Tomato', 'Turquoise', 'Violet', 'Wheat', 'White', 'WhiteSmoke', 'Yellow', 'YellowGreen']

colors = cycle(colors)

col_dict = {}
cog_col_dict = {}
#cogname = sys.argv[2]
cogname = scriptpath + '../database/cog/cog2014/cognames2003-2014.tab'

if sys.version_info[0] > 2:
    f = open(cogname, 'r', encoding='windows-1252')
else:
    f = open(cogname, 'r')
if sys.version_info[0] > 2:
    print next(f)[:-1] + '\t' + 'color'
else:
    print f.next()[:-1] + '\t' + 'color'
for i in f:
    qid, j, annot = i[:-1].split('\t')[: 3]

    # random select color for each func id of cog
    if j in col_dict:
        pass
    else:
        #col_dict[j] = random.sample(colors, 1)[0]
        col_dict[j] = colors.next()

    col = col_dict[j]
    cog_col_dict[qid] = [j, annot, col]


f.close()


qry = sys.argv[1]
prefix = qry[: -len('_cog_adjacency')]

f = open(qry, 'r')

_o1 = open('%s_network.sif' % prefix, 'w')


# deprecate edge features
# add the header
#_o2 = open('%s_edge.attrs' % prefix, 'w')
#_o2.write('Concurrency\n')

_o3 = open('%s_node.tab' % prefix, 'w')
_o3.write('\t'.join(['node', 'cog_category', 'function', 'color\n']))
visit = set()

# deprecate the header
f.next()
for i in f:
    j = i[: -1].split('\t')
    qid, sid, frqs, coo = j[0], j[1], j[4], j[-1]
    #_o1.write('%s nn %s\n' % (qid, sid))
    _o1.write('%s %s %s\n' % (qid, frqs, sid))
    #_o2.write('%s (nn) %s = %s\n' % (qid, sid, coo))
    if qid not in visit and qid in cog_col_dict:
        fuc, annot, col = cog_col_dict[qid]
        _o3.write('%s\t%s\t%s\t%s\n' % (qid, fuc, annot, col))
        visit.add(qid)
    if sid not in visit and sid in cog_col_dict:

        fuc, annot, col = cog_col_dict[sid]
        _o3.write('%s\t%s\t%s\t%s\n' % (sid, fuc, annot, col))
        visit.add(sid)
