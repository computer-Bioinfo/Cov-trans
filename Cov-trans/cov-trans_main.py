# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
# import gurobipy
import pysam
import argparse
import sys
import sys
sys.setrecursionlimit(5000)
import re
import os
import time
import pandas as pd
from typing import Dict, Tuple, List
from bam_io import SpanningReads
from solvetrans import MILP_solve
from segment_graph import *
import numpy as np



class segmentGraph():

    def __init__(self, breakpoints=None, splice_edges=None, adjacency_edges=None, ref=None):
        self.breakpoints_custom = sorted(list(breakpoints))
        self.nodes = {}
        self.splice_edges = splice_edges if splice_edges else []
        self.adjacency_edges = adjacency_edges if adjacency_edges else []
        breakpoints.add(ref.lengths[-1]-1)
        breakpoints = sorted(list(breakpoints))
        self.num_node = 0
        if len(breakpoints) > 0:
            for i in range(len(breakpoints) - 1):
                key = (breakpoints[i], breakpoints[i + 1])
                self.nodes[key] = self.num_node
                self.num_node += 1
        self.ref = ref
        self.numEdges = len(self.splice_edges)
        self.numSpliceEdges = len(splice_edges)
        self.egdesnum=len(self.splice_edges) + len(self.adjacency_edges)
        self.edges={}
        # print("splice_edges:", self.numSpliceEdges

    def writeGraph(self, fname) -> None:
        with open(fname, 'w') as output:
            output.write(f">nodes\t{self.num_node}\n")
            for key, node_id in self.nodes.items():
                output.write(f"{node_id}\t{key[0]}\t{key[1]}\n")

            output.write(f">edges\t{self.egdesnum}\n")
            for edge in self.splice_edges:
                # 假设 edge 是一个元组 (source, target)
                output.write(f"{edge[0]}\t{edge[1]}\n")
            for edge in self.adjacency_edges:
                # 假设 edge 是一个元组 (source, target)
                output.write(f"{edge[0]}\t{edge[1]}\n")

    def dfs(node, path, adjacency_matrix, max_crossings, current_crossings):
        # Add current node to the path
        path.append(node)

        all_nonpaths=[]
        valid_path = 1 < current_crossings <= max_crossings


        for neighbor, relation in enumerate(adjacency_matrix[node]):
            if relation == 1:  # 1 represents adjacent edge
                if neighbor not in path:
                    paths=segmentGraph.dfs(neighbor, path.copy(), adjacency_matrix, max_crossings, current_crossings)
                    all_nonpaths.extend(paths)
            elif relation == 2:  # 2 represents crossing edge
                if current_crossings + 1 <= max_crossings and neighbor not in path:
                    paths=segmentGraph.dfs(neighbor, path.copy(), adjacency_matrix, max_crossings, current_crossings + 1)
                    all_nonpaths.extend(paths)
        if valid_path and not all_nonpaths:
            # This is a valid path that ends at the last node
            return [path]
        else:
            return all_nonpaths
        return all_nonpaths

    def dfs01(node, path, adjacency_matrix, max_crossings, current_crossings, first_move=False):
        # Add current node to the path
        path.append(node)

        all_nonpaths = []
        valid_path = 1 <= current_crossings <= max_crossings

        for neighbor, relation in enumerate(adjacency_matrix[node]):
            if relation == 1:  # 1 represents adjacent edge
                if neighbor not in path and not first_move:
                    paths = segmentGraph.dfs01(neighbor, path.copy(), adjacency_matrix, max_crossings, current_crossings)
                    all_nonpaths.extend(paths)
            elif relation == 2:  # 2 represents crossing edge
                if neighbor not in path:
                    if first_move:
                        # Allow only crossing edges for the first move
                        paths = segmentGraph.dfs01(neighbor, path.copy(), adjacency_matrix, max_crossings, current_crossings + 1,
                                    first_move=False)
                    else:
                        paths = segmentGraph.dfs01(neighbor, path.copy(), adjacency_matrix, max_crossings, current_crossings + 1)
                    all_nonpaths.extend(paths)

        if valid_path and not all_nonpaths:
            # This is a valid path that ends at the last node
            return [path]
        return all_nonpaths



    def calculate_point_averages(nodes, pos_coverage):
        point_averages = {}

        # 计算每个 (start, end) 范围的平均覆盖度
        for (start, end) in nodes.keys():
            # 计算每个位置的索引
            # 注意 pos_coverage 的第一个元素对应位置 49
            # 所以索引应为 position - 49
            if start==49:
                positions = range(start, end + 1)
                indexes = [pos - 49 for pos in positions]
                # print(positions,indexes)
            else:
                positions = range(start+1, end + 1)
                # print(positions)
                indexes = [pos - 49 for pos in positions]



            # 提取这些位置的覆盖度值
            thicknesses = [pos_coverage[i] if 0 <= i < len(pos_coverage) else 0 for i in indexes]

            # 计算平均覆盖度
            if thicknesses:
                average_thickness = sum(thicknesses) / len(thicknesses)
            else:
                average_thickness = 0

            point_averages[(start, end)] = average_thickness

        return point_averages





def main(args):
    ref = pysam.FastaFile(args.fasta)

    Spanningcounter = SpanningReads(args.fasta, contig=ref.references[0])

    pos_coverage = Spanningcounter.POScoverage(args.bam)

    df_sj, base_ref, can_trans, non_trans_splice, orf1ab_trans = Spanningcounter.get_spanning_reads(args.bam,
                                                                                                    pos_coverage)

    breakpoints, splice_edges, adjacency_edges, cal_continue_edges = getsegmentGraph(base_ref, can_trans,
                                                                                     non_trans_splice, orf1ab_trans,
                                                                                     ref)

    sampleGraph = segmentGraph(breakpoints=breakpoints, splice_edges=splice_edges, adjacency_edges=adjacency_edges,
                               ref=ref)
    if args.Graph:
        sampleGraph.writeGraph(args.Graph)

    nodes = {}
    # print(adjacency_edges)
    edges = dict(list(splice_edges.items()) + list(adjacency_edges.items()))
    print("edges_num", len(edges))

    breakpoints.add(49)
    breakpoints = sorted(list(breakpoints))

    num_node = 0
    if len(breakpoints) > 0:
        for i in range(len(breakpoints) - 1):
            key = (breakpoints[i], breakpoints[i + 1])
            nodes[key] = num_node
            num_node += 1
    print("nodes_num", num_node)
    # print('nodes:',nodes)

    # cal
    calpath_num = 0
    all_calpaths = []
    for (left, right) in can_trans.keys():
        for key, value in nodes.items():
            if key[0] == left:
                node_index1 = value

        for key, value in nodes.items():
            if key[1] == right:
                node_index2 = value
        if node_index2 - 1 != node_index1:
            path_between = list(range(node_index1, node_index2 + 1))
            all_calpaths.append([0] + path_between)
        else:
            all_calpaths.append([0, node_index1, node_index2])

    for (left, right) in orf1ab_trans.keys():
        for key, value in nodes.items():
            if key[0] == left:
                node_index1 = value

        for key, value in nodes.items():
            if key[1] == right:
                node_index2 = value
        if node_index2 - 1 == node_index1:
            all_calpaths.append([0, node_index1, node_index2])
        else:
            path_between = list(range(node_index1, node_index2 + 1))
            all_calpaths.append([0] + path_between)

    # 将路径中的点序号变成点区域
    reverse_nodes = {value: key for key, value in nodes.items()}
    all_calpaths_with_keys = []
    for path in all_calpaths:
        try:
            path_with_keys = [reverse_nodes[node] for node in path[0:]]
            all_calpaths_with_keys.append(path_with_keys)
        except KeyError as e:
            print(f"KeyError: {e}, Path: {path}")

    # Print the paths
    for path in all_calpaths_with_keys:
        calpath_num += 1
        # print("cal-Path:", path)
    print("Candidate_cal_path_num", calpath_num)


    cal_transcripts_all=[]
    for path in  all_calpaths_with_keys:
        contig = "NC_045512.2"  # 假设所有转录本都属于同一个contig
        transcript = {
            "contig": contig,
            "path": path,
            "type": "cal"
        }
        cal_transcripts_all.append(transcript)
    # for transcript in transcripts_all:
    #     print(transcript)





    Node_edge = [[0] * num_node for i in range(num_node)]

    for i in range(num_node):
        if i + 1 < num_node:
            Node_edge[i][i + 1] = 1

    # 处理 splice_edges 中的节点连接
    for (left, right) in non_trans_splice.keys():
        node_index1 = None
        node_index2 = None

        # 寻找与 left 匹配的节点索引
        for key, value in nodes.items():
            if key[1] == left:
                node_index1 = value
                break  # 找到后立即跳出循环

        # 寻找与 right 匹配的节点索引
        for key, value in nodes.items():
            if key[0] == right:
                node_index2 = value
                break  # 找到后立即跳出循环

        # 如果找到匹配的索引，更新 Node_edge
        if node_index1 is not None and node_index2 is not None:
            Node_edge[node_index1][node_index2] = 2

    start_node = 0
    max_crossings_allowed = 2
    initial_crossings = 0
    initial_path = []
    # all_nonpaths = segmentGraph.dfs(start_node, initial_path, Node_edge, max_crossings_allowed, initial_crossings)
    all_nonpaths = segmentGraph.dfs01(start_node, initial_path, Node_edge, max_crossings_allowed, initial_crossings, first_move=True)
    nonpath_num = 0

    all_nonpaths_with_keys = []
    for path in all_nonpaths:

        try:
            path_with_keys = [reverse_nodes[node] for node in path[0:]]
            all_nonpaths_with_keys.append(path_with_keys)
        except KeyError as e:
            print(f"KeyError: {e}, Path: {path}")

    for path in all_nonpaths_with_keys:
        nonpath_num += 1
        # print("non-Path:", path)
    print("Candidate_nonpath_num", nonpath_num)
    # for path in all_nonpaths:
    #     nonpath_num += 1
        # print("non-Path:", path)
    # print("noncal_path_num",nonpath_num)

    non_transcripts_all = []

    for path in all_nonpaths_with_keys:
        contig = "NC_045512.2"
        transcript = {
            "contig": contig,
            "path": path,
            "type": "non"
        }
        non_transcripts_all.append(transcript)
        # print(non_transcripts_all)

    # 创建 Path_node1 和 Path_node2
    Path_node1 = [[0] * len(nodes) for _ in range(len(all_nonpaths))]
    Path_node2 = [[0] * len(nodes) for _ in range(len(all_calpaths))]

    for row, path in enumerate(all_nonpaths):
        for node in path:
            Path_node1[row][node] = 1

    for row, path in enumerate(all_calpaths):
        for node in path:
            Path_node2[row][node] = 1

    # 构建边的索引字典以优化查找速度
    edge_index_map = {edge: idx for idx, edge in enumerate(edges.keys())}

    # 创建 Path_edge1 和 Path_edge2
    Path_edge1 = [[0] * len(edges) for _ in range(len(all_nonpaths_with_keys))]
    Path_edge2 = [[0] * len(edges) for _ in range(len(all_calpaths_with_keys))]

    for row, path in enumerate(all_nonpaths_with_keys):
        path_edges = [(node1[1], node2[0]) for node1, node2 in zip(path[:-1], path[1:])]
        for edge in path_edges:
            if edge in edge_index_map:
                Path_edge1[row][edge_index_map[edge]] = 1

    for row, path in enumerate(all_calpaths_with_keys):
        path_edges = [(node1[1], node2[0]) for node1, node2 in zip(path[:-1], path[1:])]
        for edge in path_edges:
            if edge in edge_index_map:
                Path_edge2[row][edge_index_map[edge]] = 1



    Path = [[0] for _ in range(len(all_nonpaths))]

    point_averages= segmentGraph.calculate_point_averages(nodes, pos_coverage)
    # print(point_averages)
    node1_flow = [[0]* len(nodes)for _ in range(1)]
    for index, (key, value) in enumerate(point_averages.items()):
        if index < len(nodes):
            node1_flow[0][index] = value
    # print(node1_flow)

    sorted_dict = dict(sorted(edges.items(), key=lambda item: item[0][1]))



    edge_flow = [[0] * len(sorted_dict) for _ in range(1)]  # Initialize a 2D list with zeros

    keys_list = list(sorted_dict.keys())
    # print(sorted_dict)
    # Iterate over sorted_dict items
    for edge, value1 in sorted_dict.items():
        edge_index1 = keys_list.index(edge)
        edge_flow[0][edge_index1] += value1
        for edge2, value in sorted_dict.items():
            if edge[1]==edge2[0] and edge2[0]==edge2[1]:
                edge_index2 = keys_list.index(edge2)
                edge_flow[0][edge_index2] -= value1
                break

        if edge[0]==edge[1]:
            edge_index3 = keys_list.index(edge)
            edge_flow[0][edge_index3] = pos_coverage[edge[0]-49]




    # print("edge_flow:",edge_flow)




    solve = MILP_solve(graph=sampleGraph, threads=args.threads, ref_length=ref.lengths[0])
    out_dir = ""
    outfile = out_dir + "milp_solver-trans.py"

    non_results=solve.solve_milp(edges, Path_edge2, Path_edge1, Path_node1, Path_node2, all_calpaths, all_nonpaths, Path,
                        Node_edge, node1_flow, edge_flow, outfile)

    if args.calTransGTF:
            solve.writecalgtf(args.calTransGTF, cal_transcripts_all, ref.references[0])
    if args.noncalTransGTF:
            solve.writenongtf(args.noncalTransGTF, non_transcripts_all, non_results, ref.references[0])



def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def get_command():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam", type=str)
    parser.add_argument("-f", "--fasta", type=str,help="fasta file", required=True)
    parser.add_argument("--Graph", type=str, help="output graph")
    parser.add_argument("--calTransGTF", type=str,help="output canonical transcripts file in GTF format")
    parser.add_argument("--noncalTransGTF", type=str, help="output nnon-canonical transcripts file in GTF format")
    parser.add_argument("--threads", type=int, default=1,help='number of threads allowed to be used [1]')
    parser.set_defaults(verbose=True)
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    args.vcf = None
    return args



if __name__ == '__main__':
    print("============================================================")
    print("=====================CoV-trans==============================")
    print("============================================================")

    arguments = get_command()
    main(arguments)
