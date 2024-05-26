import pandas as pd
from typing import Dict


def getsegmentGraph(base_ref, can_trans, non_trans, orf1ab_trans, ref):
    splice_edges = {}
    adjacency_edges = {}
    breakpoints = set()

    for key in can_trans.keys():
        breakpoints.add(key[0])
        breakpoints.add(key[1])
    for key in orf1ab_trans.keys():
         breakpoints.add(key[0])
         breakpoints.add(key[1])
    for key in non_trans.keys():
        breakpoints.add(key[0])
        breakpoints.add(key[1])





    for key, value in can_trans.items():
        if (key[0], key[1]) not in splice_edges.keys():
            splice_edges[(key[0],key[1])] = value['splice']
            splice_edges[(84,key[0])] = value['splice']
        else:
            splice_edges[(key[0], key[1])] += value['splice']
            splice_edges[(84, key[0])] += value['splice']

    for key, value in non_trans.items():
        if (key[0], key[1]) in splice_edges.keys():
            splice_edges[(key[0], key[1])] += value['non']
        else:
            splice_edges[(key[0], key[1])] = value['non']



    for key, value in orf1ab_trans.items():
        if (key[0], key[1]) not in splice_edges.keys():
            splice_edges[(key[0], key[1])] = value['orf1ab']
            splice_edges[(84, key[0])] = value['orf1ab']




    for pos in breakpoints:
        if pos != 49 and pos != 84:
            # adjacency_edges[(pos, pos)] = {ref.fetch(ref.references[0], pos-1, pos): 0}
            adjacency_edges[(pos, pos)] = 0

    # print(splice_edges)
    return breakpoints, splice_edges, adjacency_edges


