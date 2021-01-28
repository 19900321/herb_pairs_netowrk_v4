# case study1: huangqi + gancao
from collections import defaultdict
from itertools import combinations
import pandas as pd
from proximity_key import *
from generate_objects import *
from disease import *
from permuation_herbs import *
import matplotlib.pyplot as plt
import networkx as nx
from process.drug_annotation_pipline import get_symble_from_entries
g_obj.get_degree_binning(1001)
from herb_ingre_tar import annotaion_herbs_from_mqsql, annotaion_ingredients_from_mqsql


# use herb:center, ingre:closest as final distance, 6.25
class Herb_Pair_network:
    def __init__(self, herb_distance_obj, herb_from_name, herb_to_name, ingre_method, herb_method, herb_info):
        import pandas as pd
        self.herb_from_name = herb_from_name
        self.herb_to_name = herb_to_name
        self.herb_from = herb_info.herb_pinyin_dic[self.herb_from_name]
        self.herb_to = herb_info.herb_pinyin_dic[self.herb_to_name]
        self.distance_network = herb_distance_obj.herb_herb_dis_all(self.herb_from, self.herb_to)
        self.ingre_method = ingre_method
        self.herb_method = herb_method
        self.ingre_distance_dict_list = self.distance_network[self.ingre_method]['two_level']['length_dict']
        self.herb_level_distance = self.get_herb_level_distance()
        self.pd_ingre_pairs_dis = self.get_ingre_distance_pd()
        self.herb_ingre_id_name = self.get_herb_ingre_id_name_dict()
        self.pd_herb_ingre = self.get_herb_ingre_pd()
        self.center_ingredients = self.get_center_ingredients()
        self.ingre_z_dict = defaultdict()
        self.herb_z_dict = defaultdict()
        self.herb_ingre_z_dict = defaultdict()

    def get_herb_level_distance(self):
        return self.distance_network[self.ingre_method]['two_level']['distances'][self.herb_method]

    def name_find(self, ingre_id):
        return ingredients_obj.ingredients_info(ingre_id)['name']

    def name_trans_herb(self, herb_id):
        return herb_info.pinyin_herbid_dic[herb_id]

    def get_ingre_distance_pd(self):
        pd_ingredients = pd.concat([pd.DataFrame.from_dict(i, orient='index').stack().reset_index()
                                    for i in self.ingre_distance_dict_list])
        pd_ingredients.columns = ['node_from', 'node_to', 'distance']
        pd_ingredients['node_from_name'] = pd_ingredients['node_from'].apply(self.name_find)
        pd_ingredients['node_to_name'] = pd_ingredients['node_to'].apply(self.name_find)
        return pd_ingredients

    def get_herb_ingre_id_name_dict(self):
        return {self.herb_from: {ingre: self.name_find(ingre) for ingre in self.ingre_distance_dict_list[2].keys()},
                self.herb_to: {ingre: self.name_find(ingre) for ingre in self.ingre_distance_dict_list[3].keys()}}

    def get_herb_ingre_pd(self):
        pd_herb_ingre = pd.DataFrame.from_dict({k: v.keys() for k, v in self.herb_ingre_id_name.items()},
                                               orient='index').stack().reset_index()
        pd_herb_ingre.columns = ['node_from', 'index_add', 'node_to']
        pd_herb_ingre = pd_herb_ingre.drop(['index_add'], axis=1)
        pd_herb_ingre = pd_herb_ingre.rename(columns={'node_from':'herb_id', 'node_to':'ingredient_id'})
        pd_herb_ingre['herb_name'] = pd_herb_ingre['herb_id'].apply(self.name_trans_herb)
        pd_herb_ingre['ingredient_name'] = pd_herb_ingre['ingredient_id'].apply(self.name_find)
        pd_ingre_anno = annotaion_ingredients_from_mqsql()
        pd_ingre_anno['Ingredient id'] = pd_ingre_anno['Ingredient id'].apply(lambda x:'I'+str(x))
        pd_herb_ingre = pd.merge(pd_herb_ingre,
                                 pd_ingre_anno,
                                 left_on='ingredient_id',
                                 right_on='Ingredient id')
        return pd_herb_ingre


    def get_center_ingredients(self):
        ingredients_from = get_center_one(self.ingre_distance_dict_list[2].keys(), self.ingre_distance_dict_list[2])
        ingredients_to = get_center_one(self.ingre_distance_dict_list[3].keys(), self.ingre_distance_dict_list[3])
        return {self.herb_from_name: {ingre: self.name_find(ingre) for ingre in ingredients_from},
                self.herb_to_name: {ingre: self.name_find(ingre) for ingre in ingredients_to}}


    def get_disease_herb_ingre_z(self, disease_obj, disease, herb, distance_method, herb_ingre_dict, ingre_tar_dict,
                                 random_time, seed):
        ingre_disease_dict = disease_obj.cal_disease_herb_ingre_z_score(disease, herb, distance_method, herb_ingre_dict, ingre_tar_dict,
                                                                        random_time, seed)
        ingre_disease_pd = pd.DataFrame.from_dict(ingre_disease_dict, orient='index',
                                                  columns=['d', 'z', 'm', 's', 'p_val'])
        ingre_disease_pd['herb'] = self.name_trans_herb(herb)
        ingre_disease_pd['herb_id'] = herb
        ingre_disease_pd['ingre_id'] = ingre_disease_pd.index
        ingre_disease_pd['ingre_name'] = ingre_disease_pd['ingre_id'].apply(self.name_find)
        return ingre_disease_dict, ingre_disease_pd


# generate all methods in center,
def get_ingre_dis(herb_distance_obj, herb_info, method_list):
    dis_list = []
    center_dict = defaultdict()
    for m in method_list:
        huang_gan_network = Herb_Pair_network(herb_distance_obj, 'HUANG QI', 'GAN CAO', m, 'center', herb_info)
        huang_gan_dis_pd = huang_gan_network.pd_ingre_pairs_dis
        dis_list.append(huang_gan_dis_pd)
        distance = huang_gan_network.herb_level_distance
        huang_gan_network.center_ingredients.update({'distance': distance})
        center_dict[m] = huang_gan_network.center_ingredients
    return dis_list, center_dict


def prepare_center_distance_list(center_dict):
    center_pd = pd.DataFrame.from_dict(center_dict, orient='index')
    center_pd['methods'] = list(center_pd.index)
    center_pd_2 = center_pd.melt(id_vars=['methods', 'distance'], var_name='herb')
    center_prepared = []
    for index, value in center_pd_2.iterrows():
        dict_col = value['value']
        for k, v in dict_col.items():
            col_S = [value['methods']] + [value['herb']] + [k] + [v] + [value['distance']]
            center_prepared.append(col_S)
    center_prepared_new = pd.DataFrame(center_prepared, columns=['methods',  'herb', 'ingre_id', 'ingre_name', 'distance'])
    return center_prepared_new


def prepare_center_distance_list_2(center_dict, mean_pd_top):
    center_pd = pd.DataFrame.from_dict(center_dict, orient='index')
    center_pd['Ingredient-level distance type'] = list(center_pd.index)
    center_pd['HUANG QI_2'] = center_pd['HUANG QI'].apply(lambda x: list(x.values()))
    center_pd['GAN CAO_2'] = center_pd['GAN CAO'].apply(lambda x: list(x.values()))

    mean_pd_top = mean_pd_top[mean_pd_top['Herb-level distance type'] == 'center']
    pre_pd = pd.merge(mean_pd_top, center_pd, 'inner', left_on='Ingredient-level distance type', right_on='Ingredient-level distance type')
    # prepare to paper format
    center_pd = pre_pd.copy()
    center_pd['center of HUANG QI'] = center_pd['HUANG QI_2'].apply(lambda x: ';'.join(x))
    center_pd['center of GAN CAO'] = center_pd['GAN CAO_2'].apply(lambda x: ';'.join(x))
    center_pd = center_pd.rename(
        columns={'random': 'Distance for random herb pairs', 'distance': 'Distance for HuangQiGanCao', 'top': 'Distance for top herb pairs',
                 'herb_method': 'Herb-level distance type',
                 'ingre_method': 'Ingredient-level distance type'})
    center_pd = center_pd[['center of HUANG QI',
                           'center of GAN CAO',
                           'Ingredient-level distance type',
                           'Distance for HuangQiGanCao',
                           'Distance for top herb pairs',
                           'Distance for random herb pairs']]
    center_pd = center_pd.round(decimals=2)
    return pre_pd, center_pd


# prepare cytoscape
def ingre_target_network(ingre_1_list, ingre_2_list, g_obj, ingredients_obj, how):
    ingre_name_1 = ';'.join([ingredients_obj.ingre_id_name_dict_value[ingre_1] for ingre_1 in ingre_1_list])
    ingre_name_2 = ';'.join([ingredients_obj.ingre_id_name_dict_value[ingre_2] for ingre_2 in ingre_2_list])
    t1 = []
    for ingre_1 in ingre_1_list:
        t1 += ingredients_obj.ingre_tar_dict[ingre_1]

    t2 = []
    for ingre_2 in ingre_2_list:
        t2 += ingredients_obj.ingre_tar_dict[ingre_2]

    t_s = list(set(t1 + t2))
    t_all = []
    for t_1, t_2 in combinations(t_s, 2):
        nodes_related = nx.shortest_path(g_obj.G, t_1, t_2)
        t_all += nodes_related

    t_all = list(set(t_all))
    g_sub_focus = nx.Graph(g_obj.G.subgraph(t_all))

    # add ingre-target
    edge_1 = [(ingre_1, i) for i in list(ingredients_obj.ingre_tar_dict[ingre_1]) for ingre_1 in ingre_1_list ]
    edge_2 = [(ingre_2, i) for i in list(ingredients_obj.ingre_tar_dict[ingre_2]) for ingre_2 in ingre_2_list]
    ingre_list = ingre_1_list + ingre_2_list
    g_sub_focus.add_edges_from(edge_1)
    g_sub_focus.add_edges_from(edge_2)


     # change lable
    node_label = {t:get_symble_from_entries(t[1:]) if t not in ingre_list else ingredients_obj.ingre_id_name_dict_value[t] for t in g_sub_focus.nodes }
    print(node_label.values())


    # generate node color dict by group
    node_color_dict = defaultdict()
    for t in g_sub_focus.nodes:
        if t in t1 and t not in t2:
            node_color_dict[node_label[t]] = 'red'
        elif t in t2 and t not in t1:
            node_color_dict[node_label[t]] = 'green'
        elif t in t2 and t in t1:
            node_color_dict[node_label[t]] = 'yellow'
        elif t in ingre_list:
            node_color_dict[node_label[t]] = 'orange'
        else:
            node_color_dict[node_label[t]] = 'grey'

    H = nx.relabel_nodes(g_sub_focus, node_label)
    node_color = []
    for n in H.nodes:
        node_color.append(node_color_dict[n])
    options = {
        'node_color': node_color,
        'node_size': 2000,
        'width': 2,
        'alpha': 0.9,
        'font_size':12
    }
    ax = plt.figure(figsize=(15, 10))
    nx.draw(H, with_labels= True,
            edge_color = 'grey',
            font_weight='bold',
            **options)
    plt.title('center ingredient nodes {} (yellow node) and {} (blue nodes)'.format(ingre_name_1, ingre_name_1))
    print('center ingredient nodes {} of GAN CAO (red node) and {} of HUANG QI (green nodes)'.format(ingre_name_1, ingre_name_2))
    if how == 'save_figure':
        plt.savefig('figure/Figure 6.png')
    elif how == 'plot_figure':
        plt.show()
    return H


def plot_center_network(center_dict, method, g_obj, ingredients_obj, how):
    ingre_1_list = list(center_dict[method]['GAN CAO'].keys())
    ingre_2_list = list(center_dict[method]['HUANG QI'].keys())
    ingre_target_network(ingre_1_list, ingre_2_list, g_obj, ingredients_obj, how)


def prepare_table_s5(herb_pair_ingre):
    herb_anno = annotaion_herbs_from_mqsql()
    herb_anno .columns = ['TCMID herb ' + c for c in herb_anno.columns]
    col_replace = {
        'herb_name': 'TCMID herb pinyin name',
        'ingredient_name': 'TCMID ingredient_name',
        'Ingredient Website': 'TCMID Ingredient Website',
        'ingredient_id_x': 'TCMID ingredient id',
        'cid': 'PubChem id',
        'name': 'Stitch name'
    }
    herb_pair_ingre = herb_pair_ingre.rename(columns = col_replace)
    selected_columns = [
        'TCMID herb pinyin name',
        'TCMID ingredient_name',
        'TCMID ingredient id',
        'TCMID Ingredient Website',
        'PubChem id',
        'PubChem ID website',
        'Pubchem_inchikey',
        'Pubchem_iupac',
        'Pubchem_molecular_formula',
        'Pubchem_smiles',
        'Stitch_cid_m',
        'Stitch_cid_s',
        'Stitch name'
    ]
    herb_pair_ingre = herb_pair_ingre[selected_columns]
    herb_pair_ingre = pd.merge(herb_anno,
             herb_pair_ingre,
             left_on='TCMID herb Pinyin Name',
             right_on='TCMID herb pinyin name',
             how='right')
    herb_pair_ingre['TCMID Ingredient Website'] = herb_pair_ingre['TCMID Ingredient Website'].apply(lambda x:x.replace('http://www.megabionet.org/',
                                                                                                                  'http://119.3.41.228:8000/'))
    herb_pair_ingre['TCMID herb website'] = herb_pair_ingre['TCMID herb website'].apply(lambda x:x.replace(
        'http://www.megabionet.org/',
        'http://119.3.41.228:8000/'))
    herb_pair_ingre = herb_pair_ingre.drop(columns=['TCMID herb pinyin name'])
    herb_pair_ingre = herb_pair_ingre.rename(columns={'TCMID herb herb-id':'TCMID herb id'})
    return herb_pair_ingre


def main():
    # huangqi example

    # g_obj.get_degree_equivalents('T441531')
    # pairs = g_obj.get_random_equivalents_set(['T441531','T2537'], 100, 45)
    # huang_gan_network = Herb_Pair_network(herb_distance_obj, 'HUANG QI', 'GAN CAO', 'shortest', 'center', herb_info)
    # huang_gan_network.pd_ingre_pairs_dis.to_csv('result/huangqitang/ingre_pairs_dis.csv')
    # herb_pair_ingre = huang_gan_network_closest.pd_herb_ingre

    # table_S5 = prepare_table_s5(herb_pair_ingre)
    # table_S5.to_csv('result/huangqitang/herb_ingre_pairs.txt', spe= '\t)
    # table_S5.to_csv('result/Table_S5.txt', spe='\t)
    # table_S5.to_csv('table/Table_S5.txt', spe='\t)

    # huang_gan_network_closest = Herb_Pair_network(herb_distance_obj, 'HUANG QI', 'GAN CAO', 'closest', 'center',
    #                                              herb_info)

    # huang_gan_network_closest.pd_ingre_pairs_dis.to_csv('result/huangqitang/ingre_pairs_closest.csv')

    # method_list = ['closest', 'shortest', 'center', 'separation', 'kernel']
    # dis_pd_all, center_dict = get_ingre_dis(herb_distance_obj, herb_info, method_list)
    # mean_pd_top = pd.read_csv('result/top_random/mean.csv')

    # # pre_pd, center_pd = prepare_center_distance_list_2(center_dict, mean_pd_top)
    # # center_pd.to_csv('result/huangqi/center_records.csv')

    # pre_pd_back = prepare_center_distance_list(center_dict)
    # pre_pd_back.to_csv('result/huangqi/center_records_back.csv')

    # plot
    ingre_1 = ['I23134']
    ingre_2 = ['I13091']
    ingre_target_network(ingre_1, ingre_2, g_obj, ingredients_obj, 'save_figure')

    #plot_center_network(center_dict, 'shortest', g_obj.G, ingredients_obj)


def disease_liver():
    # g_obj.get_degree_equivadis_pd_all,center_dict = get_ingre_dis(herb_distance_obj, herb_info, method_list)
    get_shortest_path_length_between('T441531')
    # pairs = g_obj.get_random_equivalents_set(['T441531','T2537'], 100, 45)

    disease_file_name = 'source_data/CTD_D008103_genes_20200306081710.csv'
    disease = Disease(disease_file_name, g_obj)

    # huangqi_ingre = huang_gan_network.herb_ingre_id_name['H6919'].keys()
    pass

    # disease.cal_disease_ingre_dis('Liver Cirrhosis', 'I13631',
    #                               'closest', ingredients_obj.ingre_tar_dict)

    # disease.cal_disease_ingre_z_score('Liver Cirrhosis', 'I13631',
    #                                   'closest', ingredients_obj.ingre_tar_dict, 1000, 3333)
    #
    # disease.cal_disease_herb_z_score('Liver Cirrhosis', 'H9606',
    #                                   'closest', herb_obj.herb_ingretargets_dic, 1000, 3333)

    # gancao_disease_ingre_dict, gancao_disease_ingre_pd = huang_gan_network_closest.get_disease_herb_ingre_z(disease, 'Liver Cirrhosis', 'H6801', 'closest', ingredients_obj.ingre_tar_dict, 1000, 3333)
    #
    # huang_disease_ingre_dict, huang_disease_ingre_pd = huang_gan_network_closest.get_disease_herb_ingre_z(disease,
    #                                                                                                      'Liver Cirrhosis',
    #                                                                                                      'H6919',
    #                                                                                                      'closest',
    #                                                                                                      ingredients_obj.ingre_tar_dict,
    #                                                                                                      1000, 3333)
    #
    # gancao_disease_ingre_pd.to_csv('result/huangqitang/gancao_ingre_p_1000.csv')
    #
    # huang_disease_ingre_pd.to_csv('result/huangqitang/huangqi_ingre_p_1000.csv')

# prepare to network
