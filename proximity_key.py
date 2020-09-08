#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 6.6.2019 15.13
# @Author  : YINYIN
# @Site    : 
# @File    : proximity_key.py
# @Software: PyCharmlength_dict

import numpy as np
import pandas as pd
import numpy

import networkx
from functools import wraps



def memo_2(func):
    distance = {}

    @wraps(func)
    def _inner_wraps_2(G, source_id, target_id):

        if (source_id, target_id) in distance:
            result = distance[(source_id, target_id)]
        elif (target_id, source_id) in distance:
            result = distance[(target_id, source_id)]
        else:
            result = func(G, source_id, target_id)
            distance[(source_id, target_id)] = result
        return result

    return _inner_wraps_2


@memo_2
def drugs_shortest_path(G, source_id, target_id):
    return networkx.shortest_path_length(G, source_id, target_id)


def get_center_one(nodes_from, lengths_dict):
    if len(nodes_from) == 1:
        return nodes_from
    else:
        min_value = float("inf")
        center_source = []
        for node_from in nodes_from:
            inner_values = sum(list(lengths_dict[node_from].values()))
            if inner_values > min_value:
                continue
            elif inner_values == min_value:
                center_source.append(node_from)
            else:
                min_value = inner_values
                center_source = [node_from]
        return center_source


class Sets_Lengths:
    def __init__(self, nodes_from, nodes_to):
        self.nodes_from = nodes_from
        self.nodes_to = nodes_to

    def target_lengths(self, G):
        length_dict_AB = {node_from: {node_to: drugs_shortest_path(G, node_from, node_to)
                                   for node_to in self.nodes_to} for node_from in self.nodes_from}
        length_dict_BA = {node_to: {node_from: drugs_shortest_path(G, node_to, node_from)
                                        for node_from in self.nodes_from} for node_to in self.nodes_to}
        length_dict_AA = {node_from: {node_from2: drugs_shortest_path(G, node_from, node_from2)
                                   for node_from2 in self.nodes_from if node_from != node_from2}
                                        for node_from in self.nodes_from}
        length_dict_BB = {node_to: {node_to2: drugs_shortest_path(G, node_to, node_to2)
                                                    for node_to2 in self.nodes_to if node_to != node_to2}
                                        for node_to in self.nodes_to}
        return length_dict_AB, length_dict_BA, length_dict_AA, length_dict_BB

    def ingre_length(self, length_fuc, distance_method):
        length_dict_AB = {node_from: {node_to: length_fuc(node_from, node_to, distance_method)
                                           for node_to in self.nodes_to} for node_from in self.nodes_from}
        length_dict_BA = {node_to: {node_from: length_fuc(node_to, node_from, distance_method)
                                         for node_from in self.nodes_from} for node_to in self.nodes_to}
        length_dict_AA = {node_from: {node_from2: length_fuc(node_from, node_from2, distance_method)
                                           for node_from2 in self.nodes_from if node_from != node_from2}
                               for node_from in self.nodes_from}
        length_dict_BB = {node_to: {node_to2: length_fuc(node_to, node_to2, distance_method)
                                         for node_to2 in self.nodes_to if node_to != node_to2}
                               for node_to in self.nodes_to}
        return length_dict_AB, length_dict_BA, length_dict_AA, length_dict_BB



class Network_Distance:

    ''''if network is drug level, only give G is enough, if herb level, give score, length fuction,
     distance methods for ingredients '''

    def __init__(self, nodes_from, nodes_to, length_dict):
        self.nodes_from = nodes_from
        self.nodes_to = nodes_to
        self.length_dict = length_dict
        self.sets_lengths = Sets_Lengths(self.nodes_from, self.nodes_to)
        self.length_dict_AB, self.length_dict_BA, self.length_dict_AA, self.length_dict_BB = self.length_dict

    # closest is also the sAB in algorithm
    def cal_closest(self):
        values_outer = []
        for node_from in self.nodes_from:
            d = min(list(self.length_dict_AB[node_from].values()))
            values_outer.append(d)
        for node_to in self.nodes_to:
            d = min(list(self.length_dict_BA[node_to].values()))
            values_outer.append(d)
        closest_dis = numpy.mean(values_outer)
        return closest_dis

    def cal_separation_from(self):
        values_outer = []
        if len(self.nodes_from) == 1:
            return 0
        else:
            for node_from in self.nodes_from:
                d = min(list(self.length_dict_AA[node_from].values()))
                values_outer.append(d)
            separation_from_dis = numpy.mean(values_outer)
            return separation_from_dis

    def cal_separation_to(self):
        values_outer = []
        if len(self.nodes_to) == 1:
            return 0
        else:
            for node_to in self.nodes_to:
                d = min(list(self.length_dict_BB[node_to].values()))
                values_outer.append(d)
            separation_to_dis = np.mean(values_outer)
            return separation_to_dis

    def cal_separation_AB(self):
        separation_dis = self.cal_closest() - (self.cal_separation_to() +
                                                          self.cal_separation_from())/2
        return separation_dis

    #shortest
    def cal_shortest(self):
        """
        Calculate d shortest
        """
        values_outer = []
        for node_from in self.nodes_from:
            d = list(self.length_dict_AB[node_from].values())
            values_outer += d
        for node_to in self.nodes_to:
            d = list(self.length_dict_BA[node_to].values())
            values_outer += d
        cal_shortest = numpy.mean(values_outer)
        return cal_shortest

    ##kerel
    def kenel_process(self, value_list):
        inner_values = [np.exp(-i - 1) for i in value_list]
        value_e = np.log(np.mean(inner_values))
        return value_e

    def cal_kernel(self):
        """
        Calculate kernal
        """
        values_outer = []
        for node_from in self.nodes_from:
            value_e = self.kenel_process(list(self.length_dict_AB[node_from].values()))
            values_outer.append(value_e)
        for node_to in self.nodes_to:
            value_e = self.kenel_process(list(self.length_dict_BA[node_to].values()))
            values_outer.append(value_e)
        kernel_dis = -np.mean(values_outer)
        return kernel_dis

    ## center

    def get_center_one(self, nodes, lengths_dict):
        if len(nodes) == 1:
            return nodes
        else:
            min_value = float("inf")
            center_source = []
            for node_from in nodes:
                inner_values = sum(list(lengths_dict[node_from].values()))
                if inner_values > min_value:
                    continue
                elif inner_values == min_value:
                    center_source.append(node_from)
                else:
                    min_value = inner_values
                    center_source = [node_from]
            return center_source

    def cal_center(self):
        center_from = self.get_center_one(self.nodes_from, self.length_dict_AA)
        center_to = self.get_center_one(self.nodes_to, self.length_dict_BB)

        mean_from = np.mean([self.length_dict_AB[node_from][node_to]
                        for node_to in center_to for node_from in center_from])
        mean_to = np.mean([self.length_dict_BA[node_to][node_from]
                        for node_from in center_from for node_to in center_to])

        center_dis = np.mean([mean_from, mean_to])
        return center_dis

    # sepcific parameter about method
    def network_distance(self, distance_method):
        if distance_method == 'separation':
            distance = self.cal_separation_AB()
        elif distance_method == 'closest':
            distance = self.cal_closest()
        elif distance_method == 'shortest':
            distance = self.cal_shortest()
        elif distance_method == 'kernel':
            distance = self.cal_kernel()
        elif distance_method == 'center':
            distance = self.cal_center()
        return distance

    def cal_z_score(self, distance_method, random_time, network, seed):
        rand_from_list = network.get_random_equivalents_set(self.nodes_from, random_time, seed)
        rand_to_list = network.get_random_equivalents_set(self.nodes_to, random_time, seed)
        d = self.network_distance(distance_method)
        random_d_s = []
        for i in range(random_time):
            nodes_from_new = rand_from_list[i]
            nodes_to_new = rand_to_list[i]
            length_dict = Sets_Lengths(nodes_from_new, nodes_to_new).target_lengths(network.G)
            distance = Network_Distance(nodes_from_new, nodes_to_new, length_dict).network_distance(distance_method)
            random_d_s.append(distance)
        pval = float(sum(random_d_s <= d)) / random_time,
        m, s = np.mean(random_d_s), np.std(random_d_s)
        if s == 0:
            z = 0.0
        else:
            z = (d - m) / s

        return d, z, (m, s), pval

    if __name__ == '__main__':
        pass




