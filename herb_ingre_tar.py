#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 6.6.2019 14.58
# @Author  : YINYIN
# @Site    : 
# @File    : herb_ingre_tar.py
# @Software: PyCharm

import pandas as pd
from proximity_key import *
import numpy as np
from collections import defaultdict
from functools import reduce
import MySQLdb
import MySQLdb.cursors


def qc_herb(data):
    herb_ingre_dict = data.groupby('herb_id')['ingredients_id'].apply(list).to_dict()
    herb_ingre_dict = {key: [m for m in value if isinstance(m, str)] for key, value in herb_ingre_dict.items() if
                       len(value) != 0}
    herb_ingre_dict_2 = {key: [m for m in value if m[1:].isdigit() == True] for
                         key, value in herb_ingre_dict.items()}
    herb_ingre_dict3 = {key: value for key, value in herb_ingre_dict_2.items() if len(value) != 0}
    return herb_ingre_dict3


class Ingredients:
    def __init__(self, filename, score):
        self.filename = filename
        self.data = self.read_file()
        self.score = score
        self.ingredients = self.data.index.tolist()
        self.ingretar_dict_all = self.ingredients_target_dict_all()
        self.entry_ensymbl_dict = self.target_enrty_ensymbl()
        self.ingre_id_name_dict_value = self.ingre_id_name_dict()


    def ingredients_info(self, ingredient):
        ingredient_dict = defaultdict()
        ingredient_dict['stitch_name'] = self.data.loc[ingredient, 'stitch_name']
        ingredient_dict['name'] = self.data.loc[ingredient, 'Ingredient Name']
        ingredient_dict['target_ensymble'] = self.data.loc[ingredient, 'ensymble'].split(',')
        ingredient_dict['target'] = self.data.loc[ingredient, 'targets'].split(',')
        ingredient_dict['score'] = self.data.loc[ingredient, 'score'].split(',')
        return ingredient_dict


    def ingre_id_name_dict(self):
        self.data['ingre_2'] = 'I' + self.data['Ingredient id'].astype(str)
        ingre_id_name_dict_value = dict(zip(self.data['ingre_2'], self.data['Ingredient Name']))
        return ingre_id_name_dict_value


    def ingredients_target_dict_all(self):
        ingretar_dict_all = defaultdict()
        for ingredient in self.ingredients:
            ingretar_dict_all[ingredient] = self.ingredients_info(ingredient)['target']
        return ingretar_dict_all


    def ingredients_target(self, ingredient, G_nodes):
        ingredient_dict = self.ingredients_info(ingredient)
        if len(ingredient_dict['target']) == len(ingredient_dict['score']):
            target_left = ['T' + str(ingredient_dict['target'][m]) for m,n in enumerate(ingredient_dict['score'])
                        if float(n) > self.score]
            target_left = [t for t in target_left if t in G_nodes]
            ingredient_dict['left_target'] = target_left
        else:
            ingredient_dict['left_target'] = ['T' + str(t) for t in ingredient_dict['target'] if 'T' + str(t) in G_nodes]
        return ingredient_dict


    def target_enrty_ensymbl(self):
        entry_ensymbl_dict = {}
        for ingredient in self.ingredients:
            ingredient_dict = self.ingredients_info(ingredient)
            ensymble = ingredient_dict['target_ensymble']
            targets_changed = ['T' + t for t in ingredient_dict['target']]
            if len(ensymble) == len(targets_changed):
                entry_ensymbl_dict_one = dict(zip(targets_changed, ensymble))
                entry_ensymbl_dict.update(entry_ensymbl_dict_one)
        return entry_ensymbl_dict

    def ingredients_target_dict(self, G_nodes):
        self.ingre_tar_dict = defaultdict()
        for ingredient in self.ingredients:
            targets = self.ingredients_target(ingredient, G_nodes)['left_target']
            if len(targets) != 0:
                self.ingre_tar_dict[ingredient] = targets
        return self.ingre_tar_dict


    def read_file(self):
        data = pd.read_csv(self.filename, sep = '\t')
        data = data.dropna(axis=0, how='any')
        data['ingredients_id'] = 'I' + data['ingredients_id'].astype(str)
        data = data.set_index('ingredients_id')
        return data


    def ingre_ingre_dis(self, ingre_from, ingre_to, network, distance_method):
        if any(ingre not in self.ingre_tar_dict.keys() for ingre in [ingre_from, ingre_to]):
            print('{} or {} not in dictionary'.format(ingre_from, ingre_to))
            return None
        else:
            nodes_from = self.ingre_tar_dict[ingre_from]
            nodes_to = self.ingre_tar_dict[ingre_to]
            length_dict = Sets_Lengths(nodes_from, nodes_to).target_lengths(network)
            dis_obj = Network_Distance(nodes_from, nodes_to, length_dict)
            distance = dis_obj.network_distance(distance_method)
            return distance

    def ingre_ingre_dis_all(self, ingre_from, ingre_to, network):

        distance_method_list = ['separation', 'closest', 'shortest', 'kernel', 'center']

        return {method: self.ingre_ingre_dis(ingre_from, ingre_to, network, method)
                for method in distance_method_list}


class Herb:
    def __init__(self, filename):
        self.filename = filename

        self.data_all = self.read_data2()
        self.data = self.precess_data()
        self.herb_ingre_dict_all = qc_herb(self.data)
        self.herbs = self.herb_ingre_dict_all.keys()
        self.ingredients = self.herb_ingre_dict_all.values()

    def precess_data(self):
        data_left = self.data_all.dropna(subset= ['ingredients_id'], axis=0)
        return data_left


    def herb_ingre_dict_all(self):
        return qc_herb(self.data)

    def herb_ingre_dict(self, ingredients_target_dict):
        self.herb_ingre_dict = {k: [v_2 for v_2 in v if v_2 in ingredients_target_dict.keys() ]
                                for k, v in self.herb_ingre_dict_all.items()
                                }

        self.herb_ingre_dict = {k: v
                                for k, v in self.herb_ingre_dict.items() if len(v) != 0
                                }
        return self.herb_ingre_dict

    def herb_ingretargets(self, herb, ingredients_target_dict):
         ingre_target_list = [ingredients_target_dict[ingre]
                              for ingre in self.herb_ingre_dict_all[herb] if ingre in ingredients_target_dict.keys()]
         if len(ingre_target_list) == 0:
             return None
         else:
            return list(set((reduce(lambda x, y: x+y, ingre_target_list))))

    def herb_ingretargets_dic(self, ingredients_target_dict):
        self.herb_ingretargets_dic = defaultdict()
        for herb in list(self.herbs):
            self.herb_ingretargets_dic[herb] = self.herb_ingretargets(herb, ingredients_target_dict)
        return self.herb_ingretargets_dic

    def read_data2(self):
        data = pd.read_csv(self.filename, sep=',')
        data['ingredients_id'] = 'I' + data['ingredients_id'].astype(str)
        data['herb_id'] = 'H' + data['herb_id'].astype(str)
        return data

