#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 30.7.2019 13.07
# @Author  : YINYIN
# @Site    : 
# @File    : analysis_fiter_0_7.py
# @Software: PyCharm


import network_used
import numpy
import proximity_key
import pandas as pd
import herb_herb_pairs
import construct_network
import herb_ingre_tar
import os
import importlib
import random

from itertools import permutations
# PPI network
file_name = "example/toy.sif"

net_final = construct_network.get_network(file_name, only_lcc = True)
target_in_network = list(net_final.nodes.keys())
#bins_target = network_used.get_degree_binning(net_final, 100, lengths=None)


#get her ingre dictionary and ingre target dictionary
def get_basic_element_data(target_in_network,score):
    ingre_targets_dict, ingre_targets_pd = herb_ingre_tar.ingre_target_dict_filter_score(target_in_network,
                                                                                         score)
    herb_ingre_dict = herb_ingre_tar.herb_ingre_dict_generation(ingre_targets_dict, ingre_targets_pd)
    herb_herb_pairs = herb_ingre_tar.herb_herb_pairs(herb_ingre_dict)
    herb_herb_pairs_ranked = herb_herb_pairs.sort_values(by=['frequency'], ascending=False)
    return ingre_targets_dict,herb_ingre_dict,herb_herb_pairs_ranked,score

def top_random(herb_herb_pairs_ranked, number_calculate, net_final, herb_ingre_dict,
                                                           ingre_targets_dict, method_number, lengths, filter_score):
    pd_d_values_top = herb_ingre_tar.get_d_values(herb_herb_pairs_ranked, number_calculate, net_final, herb_ingre_dict,
                                                  ingre_targets_dict, method_number, lengths, random_it=False)
    pd_d_values_top['class'] = 'top'
    # get normal distribution of 200 pairs
    pd_d_values_normal = herb_ingre_tar.get_d_values(herb_herb_pairs_ranked, number_calculate, net_final, herb_ingre_dict,
                                                     ingre_targets_dict, method_number, lengths,
                                                     random_it=True)
    pd_d_values_normal['class'] = 'random'
    pd_value_com = pd.concat([pd_d_values_top, pd_d_values_normal], sort=False)
    pd_value_com['filter_score'] =  filter_score
    return pd_value_com

def random_only(herb_herb_pairs_ranked, number_calculate, net_final, herb_ingre_dict,ingre_targets_dict, method_number, lengths, filter_score):
    # get normal distribution of 200 pairs
    pd_d_values_normal = herb_ingre_tar.get_d_values(herb_herb_pairs_ranked, number_calculate, net_final, herb_ingre_dict, ingre_targets_dict, method_number, lengths, random_it=True)
    pd_d_values_normal['class'] = 'random'
    pd_d_values_normal['filter_score'] =  filter_score
    return pd_d_values_normal

pd_all_result = pd.DataFrame()


for i in [0, 0.7, 0.9]:
    ingre_targets_dict, herb_ingre_dict, herb_herb_pairs_ranked, score = get_basic_element_data(target_in_network, i)
    for i in [herb_herb_pairs]:
        importlib.reload(i)
    for lengths in [None, 'herb']:
        if lengths == None:
            method_numbers = [0]
        else:
            method_numbers = [0, 4, 5, 6, 7]
        for method_number in method_numbers:
            for j in range(10):
                pd_value_com = random_only(herb_herb_pairs_ranked, 1000, net_final, herb_ingre_dict,
                                           ingre_targets_dict,
                                           method_number, lengths, score)
                pd_value_com['shuffle_num'] = j
                pd_all_result = pd_all_result.append(pd_value_com)
            print(pd_all_result.shape)

pd_all_result.to_csv('example/all_result_random_only_repeat_50.csv')

