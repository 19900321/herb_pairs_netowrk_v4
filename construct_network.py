#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 6.6.2019 14.53
# @Author  : YINYIN
# @Site    : 
# @File    : construct_network.py
# @Software: PyCharm

## read data to network
import networkx, random, copy
import os, pickle, numpy

try:
    from scipy.stats import rankdata
except:
    print("scipy is not installed, rank-based distance methods wont work")
# get the shortest path for those have been calculated.

from functools import wraps
def memo_2(func):
    distance = {}
    @wraps(func)
    def _inner_wraps_2(G, source_id, target_id):
        if (source_id, target_id) in distance:
            result = distance[(source_id, target_id)]
        elif (target_id,source_id) in distance:
            result = distance[(target_id,source_id)]
        else:
            result = func(G, source_id, target_id)
            distance[(source_id, target_id)] = result
        return result
    return _inner_wraps_2

@memo_2
def get_shortest_path_length_between(G, source_id, target_id):
    return networkx.shortest_path_length(G, source_id, target_id)


class Construct_Network:

    def __init__(self, file_name):
        self.filename = file_name
        self.setNode, self.setEdge, self.dictNode, self.dictEdge = self.get_nodes_and_edges_from_sif_file(store_edge_type=False, delim=None, data_to_float=True)
        self.g = self.create_network_from_sif_file( use_edge_data = False, delim = None, include_unconnected=True)
        self.G = self.get_network( only_lcc=True)
        #self.all_shoretst_pairs = dict(networkx.all_pairs_shortest_path_length(self.G))

    def get_nodes_and_edges_from_sif_file(self, store_edge_type = False, delim=None, data_to_float=True):
        """
        Parse sif file into node and edge sets and dictionaries
        returns setNode, setEdge, dictNode, dictEdge
        store_edge_type: if True, dictEdge[(u,v)] = edge_value
        delim: delimiter between elements in sif file, if None
         all whitespaces between letters are considered as delim
        """
        self.setNode = set()
        self.setEdge = set()
        self.dictNode = {}
        self.dictEdge = {}
        flag = False
        f=open(self.filename)
        for line in f:
            if delim is None:
                words = line.rstrip("\n").split()
            else:
                words = line.rstrip("\n").split(delim)
            id1 = 'T' + str(words[0])
            self.setNode.add(id1)
            if len(words) == 2:
                if data_to_float:
                    score = float(words[1])
                else:
                    score = words[1]
                self.dictNode[id1] = score
            elif len(words) >= 3:
                if len(words) > 3:
                    flag = True
                id2 = 'T' + str(words[2])
                self.setNode.add(id2)
                self.setEdge.add((id1, id2))
                if store_edge_type:
                    if data_to_float:
                        self.dictEdge[(id1, id2)] = float(words[1])
                    else:
                        self.dictEdge[(id1, id2)] = words[1]
        f.close()
        if len(self.setEdge) == 0: self.setEdge = None;
        if len(self.dictNode) == 0: self.dictNode = None;
        if self.dictEdge == 0:  self.dictEdge = None;
        if flag:
            print("Warning: Ignored extra columns in the file!")
        return self.setNode, self.setEdge, self.dictNode, self.dictEdge


    def create_network_from_sif_file(self, use_edge_data = False, delim = None, include_unconnected=True):
        self.g = networkx.Graph()
        if include_unconnected:
            self.g.add_nodes_from(self.setNode)
        if use_edge_data:
            for e,w in self.dictEdge.items():
                u,v = e
                self.g.add_edge(u,v,w=w) #,{'w':w})
        else:
            self.g.add_edges_from(self.setEdge)
        return self.g

    def get_network(self, only_lcc):
        if only_lcc and not self.filename.endswith(".lcc"):
            print("Shrinking network to its LCC", len(self.g.nodes()), len(self.g.edges()))

            """
            Finds (strongly in the case of directed network) connected components of graph
            returnAsGraphList: returns list of graph objects corresponding to connected components (from larger to smaller)
            otherwise returns list of node list corresponding nodes in connected components
            """
            self.result_list = [c for c in sorted(networkx.connected_components(self.g), key=len, reverse=True)]
            self.G = networkx.subgraph(self.g, self.result_list[0]) #  NetworkX subgraph method wrapper
            print("Final shape:", len(self.G.nodes()), len(self.G.edges()))
        return self.G

    def save_lcc(self):
        network_lcc_file = self.filename + ".lcc"
        if not os.path.exists(network_lcc_file):
            f = open(network_lcc_file, 'w')
            for u, v in self.G.edges():
                f.write("%s 1 %s\n" % (u, v))
            f.close()

    if __name__ == '__main__':
        pass

    def get_degree_binning(self, bin_size):
        degree_to_nodes = {}

        for node, degree in self.G.degree():  # .iteritems(): # iterator in networkx 2.0
            degree_to_nodes.setdefault(degree, []).append(node)

        values = degree_to_nodes.keys()
        values = sorted(values)
        bins = []
        i = 0
        while i < len(values):
            low = values[i]
            val = degree_to_nodes[values[i]]
            while len(val) < bin_size:
                i += 1
                if i == len(values):
                    break
                val.extend(degree_to_nodes[values[i]])
            if i == len(values):
                i -= 1
            high = values[i]
            i += 1
            # print(i, low, high, len(val) )
            if len(val) < bin_size:
                low_, high_, val_ = bins[-1]
                bins[-1] = (low_, high, val_ + val)
            else:
                bins.append((low, high, val))
        self.bins = bins

    def get_degree_equivalents(self, node):
        d = self.G.degree(node)
        for l, h, nodes in self.bins:
            if int(l) <= d and int(h) >= d:
                mod_nodes = list(nodes)
                mod_nodes.remove(node)
        return mod_nodes

    def get_random_equivalents_set(self, nodes_selected, random_times, seed):
        if seed is not None:
            random.seed(seed)
        random_set_list = []
        for node in nodes_selected:
            equivalent_nodes_all = self.get_degree_equivalents(node)
            nodes_random = random.sample(equivalent_nodes_all, random_times)
            random_set_list.append(nodes_random)
        random_set = list(zip(*random_set_list))
        return random_set

#file_name = "example/toy.sif"
#a = Construct_Network("example/toy.sif")