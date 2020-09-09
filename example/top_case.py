#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 10.6.2019 15.21
# @Author  : YINYIN
# @Site    : 
# @File    : top_case.py
# @Software: PyCharm


from generate_objects import *
from figures_draw import herb_overlap_ingredient
import figures_draw
import permuation_herbs
from permuation_herbs import Result_distance


def change_herb_order(top_1500):
    top_1500['pairs'] = top_1500[['herb1', 'herb2']].apply(list).apply(
        lambda x: x[0] + x[1] if x[0] <= x[1] else x[1] + x[0], axis=1)
    top_1500_ranged = list(top_1500[['herb1','herb1_name', 'herb2', 'herb2_name']].astype(list).apply(lambda x:[x[0], x[1], x[2], x[3]] if x[0] <= x[2] else [x[2], x[3], x[0], x[1]], axis=1))
    top_ranged = pd.DataFrame(top_1500_ranged, columns=['herb1','herb1_name', 'herb2', 'herb2_name'])
    top_ranged['pairs'] = top_ranged['herb1'] + top_ranged['herb2']
    top_1500 = top_1500.drop(columns=['herb1','herb1_name', 'herb2', 'herb2_name'])
    top_1500_real = pd.merge(top_ranged, top_1500, left_on='pairs', right_on='pairs', sort= False)
    top_1500_real = top_1500_real.sort_values(by='frequency', ascending=False)
    top_1500_real = top_1500_real.drop_duplicates(subset=['pairs', 'Ingredient-level distance type'])
    return top_1500_real

# step 1. plot S1--fangji length, figure 2
figures_draw.plot_S1_fangji_length(fangji, 'save_figure')

figures_draw.plot_figure_2(fangji, 10, 'save_figure')

#len(fungji.herbid_frequency_dict) 79609
# !! DO NOT NEED TO RUN UNLESS GENERATE NEW RESULT
# Step 2. choose pairs and generate distance, Table S1, top_200.csv, random_10000
def get_distance_n (fangji, herb_distance_obj, n):
    top_pairs = fangji.choose_common_herbpairs(0, n+500, herb_obj.herb_ingre_dict,herb_info.pinyin_herbid_dic)
    top_result = herb_distance_obj.generator_result(top_pairs)
    top_result = top_result.sort_values(by='frequency', ascending=False)
    top_result['pairs'] = top_result['herb1'] + top_result['herb2']

    top_result_list = list(top_result['pairs'].unique())
    if len(top_result_list) > n:
        top_result_list = top_result_list[0:n]
    top_result_real = top_result[top_result['pairs'].isin(top_result_list)]

    return top_result_real

def get_random(fangji, herb_distance_obj, n):
    random_pairs = fangji.generate_random_pairs(n+500,
                                                herb_obj.herb_ingre_dict.keys(),
                                                herb_info.pinyin_herbid_dic)
    random_result = herb_distance_obj.generator_result(random_pairs)
    random_result = random_result.iloc[0:n, :]
    random_result.to_csv('result/result_random_0_{}.csv'.format(n))

# get frequency, table S1
frequency_top = pd.DataFrame.from_dict(dict(fangji.herb_frequency_dict.most_common()), orient = 'index' ).reset_index()
frequency_top.columns = ['Herb', 'frequency']
frequency_top.to_csv('result/Table S1_new.csv')
frequency_top.to_csv('Table/Table S1_new.csv')

#get the  top 200
n_top_pair = 200
top_result_real = get_distance_n(fangji, herb_distance_obj, n_top_pair)
top_result_real.to_csv('result/top_{}.csv'.format(n_top_pair))
table_S2 = top_result_real['herb1','herb1_name', 'herb2', 'herb2_name', 'frequency'].drop_duplicates()
table_S2.to_csv('result/Table S2.csv')
table_S2.to_csv('Table/Table S2.csv')

# get top 10000
n_tops = 10000
top_result_real_10000 = get_distance_n(fangji, herb_distance_obj, n_tops)
top_result_real_10000.to_csv('result/top_{}.csv'.format(n_tops))


#get_random(fangji, herb_distance_obj, 10000)

# Step 3. get the mean for top 200, generate table 1
filename_top = 'result/top_200.csv'
filename_random = 'result/result_random_0_10000.csv'
result = Result_distance(filename_top, filename_random) #47 time
file_folder = 'top_random_new'
result.save_analysis(file_folder, forbid_include=False)
result.pd_mean.to_csv('result/Table 1.csv')
result.pd_mean.to_csv('Table/Table 1.csv')


# get real table 1 with 10000
filename_top_10000 = 'result/top_10000.csv'
filename_random = 'result/result_random_0_10000.csv'
result_10000 = Result_distance(filename_top_10000, filename_random)
result_10000.get_mean()
mean_10000 = result_10000.pd_mean
mean_10000 = mean_10000.rename(columns = {'Distance for top herb pairs':'Distance for top 10000 herb pairs'})

# get 114 ingreidnets, use step 9
result_no_overlap = Result_distance('result/no_overlap_dis.csv', 'result/result_random_0_10000.csv' )
result_no_overlap.get_mean()
mean_no_overlap = result_no_overlap.pd_mean
mean_no_overlap = mean_no_overlap.rename(columns = {'Distance for top herb pairs':'Distance for top no_overlap 114 herb pairs'})

# merge means
# read
mean_200 = pd.read_csv('result/Table 1.csv')

mean_200['Distance for top 10000 herb pairs'] =  mean_10000['Distance for top 10000 herb pairs']
mean_200['Distance for top no_overlap 114 herb pairs'] =  mean_no_overlap['Distance for top no_overlap 114 herb pairs']
cols_names = ['Herb-level distance type',
              'Ingredient-level distance type',
              'Distance for top herb pairs',
              'Distance for top 10000 herb pairs',
              'Distance for top no_overlap 114 herb pairs',
              'Distance for random herb pairs',
              'p-value',
              'AUROC',
              'AUPRC',
       ]
mean_200 = mean_200[cols_names]
mean_200.to_csv('Table/Table 1.csv')


# plot figure S2
values_200_pd = pd.read_csv('result/Table S1.csv')
figures_draw.plot_S2_fangi_frequency(fangji, values_200_pd, 'save_figure')

# step 4. plot top and random distance
mean_pd = pd.read_csv('result/Table 1.csv')
figures_draw.plot_fig_3(mean_pd, 'save_figure')



# Step 5. plot the AUC PRC curve figure 4
filename_top = 'result/top_200.csv'
filename_random = 'result/result_random_0_10000.csv'
figures_draw.plot_figure_4(filename_top,
                           filename_random,
                           'center',
                           'separation',
                           'save_figure')

# step 6. plot PRC AUC bar of different methods figure 5, use R code
# mean_pd = pd.read_csv('result/Table 1.csv')
# figures_draw.plot_figure_5(mean_pd, 'save_figure')

# Step 7. plot the common ingredient, figure S3
random_pd = pd.read_csv('result/result_random_0_10000.csv')
top_pd = pd.read_csv('result/top_200.csv')
figures_draw.plot_S3_ingredient_overlap(top_pd, random_pd, herb_obj, 'save_figure')

# step 8. keep herbs with no overlap
top_pd = pd.read_csv('result/top_200.csv')
top_pd = top_pd.sort_values(by='frequency', ascending=False)
pairs_top = list(set(list(zip(top_pd['herb1'], top_pd['herb2']))))
overlap_count_top, over_list_dict_top = herb_overlap_ingredient(pairs_top, herb_obj)
no_overlap = []
for k, v in over_list_dict_top.items():
    if v != 0:
        if k[0] < k[1]:
            k_order = k[0] + k[1]
        else:
            k_order = k[1] + k[0]
        no_overlap.append(k_order)
top_pd['pairs'] = top_pd['herb1'] + top_pd['herb2']
pairs_top_no_overlap = top_pd[top_pd['pairs'].isin(no_overlap)]
pairs_top_no_overlap.to_csv('result/no_overlap_dis.csv')

# Step 9. get the mean and value of no overlap herb pairs
result_no_overlap = Result_distance('result/no_overlap_dis.csv', 'result/result_random_0_10000.csv' ) # 97 time
file_folder = 'top_random_no_overlap'
result_no_overlap.save_analysis(file_folder, forbid_include=False)


# Step 9. plot the  figure S5
mean_pd_no_overlap = pd.read_csv('result/top_random_no_overlap/mean.csv')
figures_draw.plot_fig_S5(mean_pd_no_overlap, 'save_figure')

# Step 11. Generate the literature collected pairs, get the distance PRC and UAC, plot figure S5
filename_classic = 'result/classic_distances_recommend.csv'
classic_pd = pd.read_csv(filename_classic)
table_S2 = classic_pd[['herb1','herb1_name', 'herb2', 'herb2_name', 'frequency']].drop_duplicates()
table_S2.to_csv('result/Table S3.csv')
table_S2.to_csv('Table/Table S3.csv')

filename_random = 'result/result_random_0_10000.csv'
result = Result_distance(filename_classic, filename_random)
file_folder = 'recom_random_new'
result.save_analysis(file_folder, forbid_include=False) #repeat 35 time
result.pd_mean.to_csv('result/literature_mean.csv')
result.pd_mean.to_csv('result/Table S3.csv')
result.pd_mean.to_csv('Table/Table S3.csv')

# plot figure S6， use R code
# mean_pd_literature = pd.read_csv('result/literature_mean.csv')
# figures_draw.plot_figure_S6(mean_pd_literature, 'save_figure')

# Step 12, use casestudy2 to plot HUANGQI GANCAO network
from example.casestudy2 import  get_ingre_dis, prepare_center_distance_list_2, ingre_target_network

method_list = ['closest', 'shortest', 'center', 'separation', 'kernel']
dis_pd_all, center_dict = get_ingre_dis(herb_distance_obj, herb_info, method_list)

# generate Table 2
mean_pd_top = pd.read_csv('result/Table 1.csv')
pre_pd, center_pd = prepare_center_distance_list_2(center_dict, mean_pd_top)
center_pd.to_csv('result/Table 2.csv')
center_pd.to_csv('Table/Table 2.csv')

# plot the Figure 6
figures_draw.plot_figure_6(g_obj, ingredients_obj, 'save_figure')
# plot_center_network(center_dict, 'shortest', g_obj.G, ingredients_obj)


# Step 13. check the numbers in paper

def dict_statistic(dictionary:dict):
    n = 0
    n_ob = []
    for k,v in dictionary.items():
        n += len(v)
        n_ob += v
    return len(dictionary), len(set(n_ob)), n


N_herb_formulae = len(fangji.fangji_herbid_dict) # 46929
dict_statistic(fangji.fangji_herbid_dict) # (46929, 1107, 231251)

# the number of formulae less than 20
fang_len = dict(Counter([len(v) for k, v in fangji.fangji_herb_dict.items()]))
sum([v for k, v in fang_len.items() if k <20])/sum(list(fang_len.values()))
# 0.9790
#
average_number = 231251/46929 # 4.927677981631827
fangji.herb_frequency_dict.most_common(10)
'''
(('GAN CAO',), 12518),
(('DANG GUI',), 7417),
 (('REN SHEN',), 7390),
 (('BAI ZHU',), 5259),
 (('HUANG QIN',), 4163),
 (('FANG FENG',), 4074),
 (('CHUAN QIONG',), 4007),
 (('FU LING',), 3666),
 (('CHEN PI',), 3650),
 (('MU XIANG',), 3642)]
'''

N_herb_total = len(fangji.herb_frequency_dict) #16767
N_herb_pairs = len(fangji.herb_pair_frequency_dict) #349197

a = []
for h in list(fangji.herb_pair_frequency_dict.keys()):
    a += h
N_herb_for_pair = len(set(a)) # 12129

N_herb_pairs_with_id = len(fangji.fangji_herbid_dict) # 24552
b = []
for h in list(fangji.herbid_frequency_dict.keys()):
    b += h
N_herb_for_pair_with_id = len(set(b)) # 1076



# herb-ingredient
herb_ingredient_total = len(herb_obj.herb_ingre_dict_all) #8199

# 17753 pairs between 4415 herbs and 3917 ingredient
dict_statistic(herb_obj.herb_ingre_dict) #(4415, 3917, 17753)

# targets related
# 4330 ingredients with 3171 targtes 25050 pairs
dict_statistic(ingredients_obj.ingre_tar_dict)# (4330, 3171, 25050)


# frequency range of top 200
top_pd = pd.read_csv('result/top_200.csv')
top_pd = top_pd.sort_values(by='frequency', ascending=False)
list(top_pd['frequency'])[0] # 3846
list(top_pd['frequency'])[-1] # 358


# fangji
fangji.herbid_frequency_dict.most_common(10)
'''
[(('GAN CAO', 'REN SHEN'), 3846),
 (('DANG GUI', 'GAN CAO'), 2907),
 (('BAI ZHU', 'REN SHEN'), 2599),
 (('BAI ZHU', 'GAN CAO'), 2578),
 (('CHUAN QIONG', 'DANG GUI'), 2296),
 (('GAN CAO', 'HUANG QIN'), 2265),
 (('DANG GUI', 'REN SHEN'), 2097),
 (('CHEN PI', 'GAN CAO'), 2002),
 (('FANG FENG', 'GAN CAO'), 1960),
 (('FU LING', 'GAN CAO'), 1817)]
 '''

#for the 349197 herb pairs, the majority of them (99.9%) occurred in less than 100 herbal formulae
n_100 = 0
n_500 = 0
for k in fangji.herb_pair_frequency_dict.values():
    if k < 100:
        n_100 += 1
    if k> 500:
        n_500 +=1
n_100/len(fangji.herb_pair_frequency_dict) # 0.9944157595855635

len(fangji.herb_pair_frequency_dict) - n_100 # 1950
print(n_500)

# the number of unique herb from top 200 herb pairs
herb_unique = set(list(top_pd['herb1'].unique()) + list(top_pd['herb2'].unique())) # 61
len(herb_unique)
# average number of ingredients
n_i = 0
for h in herb_unique:
    ingredients = dict(herb_obj.herb_ingre_dict)[h]
    n_i += len(ingredients)
n_i/len(herb_unique)

from figures_draw import herb_overlap_ingredient
top_pd = top_pd.sort_values(by='frequency', ascending=False)
pairs_top = list(set(list(zip(top_pd['herb1'],top_pd['herb2']))))
pairs_random = list(set(list(zip(random_pd['herb1'],random_pd['herb2']))))
overlap_count_top, over_list_dict_top = herb_overlap_ingredient(pairs_top, herb_obj)
overlap_count_random, over_list_dict_random = herb_overlap_ingredient(pairs_random, herb_obj)
n_0_common_ingre = overlap_count_top[0]

n_1_more_common_ingre = sum(overlap_count_top.values()) - n_0_common_ingre
n_3_more_common_ingre = sum(overlap_count_top.values()) - \
                        n_0_common_ingre -overlap_count_top[1] -\
                        overlap_count_top[2] - \
                        overlap_count_top[3]
n_1_more_common_ingre / sum(overlap_count_top.values())

n_1_random = sum(overlap_count_random.values()) - overlap_count_random[0]

n_1_random/sum(overlap_count_random.values())


# qiang huo and du huo, ('H4776', 'QIANG HUO', 'H4449', 'DU HUO'
max_ingre = max(over_list_dict_top.items(), key=operator.itemgetter(1))
dict(fangji.herbid_frequency_dict)[('H4449', 'H4776')] #  522


from figures_draw import show_detail_one_pair
print(show_detail_one_pair(herb_info, herb_obj, max_ingre[0]))
'''
'H4449', 'DU HUO', 'H4776', 'QIANG HUO', ['gamma-amin.yri.', 'Camphor', 'columbianetin', 'guaiol', 'guanidinium', 'isoimperatorin', 'isopimpinellin', 'nodakenin', 'scopoletin', 'osthole'])
'''

## show huangqi gancao common ingredient, huangqi H6919,  gancao H6801,
tabel_S1 = pd.read_csv('result/Table S1.csv')
tabel_S1[tabel_S1['pairs']=='H6801H6919']
huangqi_gancao_set =  (( 'H6801', 'H6919'), 1411)
huangqi_gancao_set =  ( 'H6801', 'H6919')
print(show_detail_one_pair(herb_info, herb_obj, huangqi_gancao_set))
'''
('H6801', 'GAN CAO', 'H6919', 'HUANG QI', ['formononetin', 'clionasterol', 'clionasterol'])
'''


# PRC AUC
mean_pd = pd.read_csv('result/Table 1.csv')
mean_pd[mean_pd['p-value'] <0.05].shape[0] # 16
mean_pd['difference'] = mean_pd['Distance for random herb pairs'] - mean_pd['Distance for top herb pairs']

best_difference = mean_pd[mean_pd['difference'] == max(mean_pd['difference'])]
mean_pd[mean_pd['p-value'] == min(mean_pd['p-value'])]


# common ingredients
n_1_more_common_ingre = sum(overlap_count_top.values()) - n_0_common_ingre
n_3_more_common_ingre = sum(overlap_count_top.values()) - n_0_common_ingre -overlap_count_top[1] -overlap_count_top[2] - overlap_count_top[3]

#PRC AUC

# AVERAGE
average_auc = np.mean(mean_pd['AUROC']) # 0.65
average_prc = np.mean(mean_pd['AUPRC']) # 0.72


# center level comparison
mean_pd.groupby('Herb-level distance type')['AUROC'].mean()
'''
center        0.803612
closest        0.792437
kernel        0.504596
separation    0.673833
shortest      0.462200
'''

mean_pd.groupby('Herb-level distance type')['AUPRC'].mean()
'''
center        0.828411
closest        0.797316
kernel        0.634952
separation    0.738881
shortest      0.607763
'''

# The literature collected pairs
mean_pd_literature = pd.read_csv('result/Table S3.csv')
average_auc_liter = np.mean(mean_pd_literature['AUROC']) # 0.62
average_prc_liter = np.mean(mean_pd_literature['AUPRC']) # 0.65
# max of literature
liter_max = mean_pd_literature[mean_pd_literature['AUPRC'] == mean_pd_literature['AUPRC'].max()]
liter_max[['Herb-level distance type',	'Ingredient-level distance type'
]]
liter_max[['AUPRC',	'AUROC']]

# THE common herbs between top 200 and literatures
tabel_S1 = pd.read_csv('result/Table S1.csv')
table_S2 = pd.read_csv('result/Table S2.csv')

tabel_S1['pairs'] = tabel_S1['herb1'] + tabel_S1['herb2']

table_S2['pairs'] = table_S2['herb1'] + table_S2['herb2']

not_common_herb = tabel_S1[~tabel_S1['pairs'].isin(table_S2['pairs'])]

danggui = not_common_herb[(not_common_herb ['herb1']=='H2538')|(not_common_herb ['herb2']=='H2538')]
len([i  for i in list(tabel_S1['pairs']) if i in list(table_S2['pairs'])]) # 32

# huangqi and gancao

table_S4 = pd.read_csv('result/Table S4.csv')
table_S4.groupby('Herb name')['Ingredient id'].count()

'''
GAN CAO     27
HUANG QI    15
'''

#


#

# # result test, skip next steps.
# import permuation_herbs
# from permuation_herbs import *
#filename_top = 'result/top_1_1500.csv'
#top_200 = pd.DataFrame(top_pairs, columns=['herb1','herb1_name', 'herb2','herb2_name', 'frequency'])

# top_200['pairs'] = top_200[['herb1', 'herb2']].astype(list).apply(lambda x:x[0]+x[1] if x[0] <= x[1] else x[1]+x[0], axis=1)
#
# dict_200 = dict(zip(top_200['pairs'],top_200['frequency']))
#top_1500 = pd.read_csv(filename_top)

#top_1500['pairs'] = top_1500[['herb1', 'herb2']].astype(list).apply(lambda x:x[0]+x[1] if x[0] <= x[1] else x[1]+x[0], axis=1)
#top_1500['frequency'] = top_1500['pairs'].apply(lambda x:dict_200[x] if x in dict_200 else None)




# top_1500_real = change_herb_order(top_1500)
# top_200 = list(top_1500_real['pairs'].unique())[0:200]
# top_200_pd = top_1500_real[top_1500_real['pairs'].isin(top_200)]
# top_200_pd = top_200_pd[['separation',
#                          'closest',
#                          'shortest',
#                          'kernel',
#                          'center',
#                          'Ingredient-level distance type',
#                          'herb1',
#                          'herb1_name',
#                          'herb2',
#                          'herb2_name','frequency']]
# top_1500_real.to_csv('result/top_1_1500_new.csv')
# top_200_pd.to_csv('result/top_200.csv')