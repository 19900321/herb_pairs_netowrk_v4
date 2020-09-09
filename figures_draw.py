#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 5.8.2019 16.25
# @Author  : YINYIN
# @Site    : 
# @File    : figures_draw.py
# @Software: PyCharm


# result = pd.read_excel('whole table.xlsx')
# prepare the result to the format of analysis
import pandas as pd
import matplotlib.pyplot as plt
import logging
import operator
import itertools
import numpy as np
from collections import Counter
from functools import reduce
from collections import defaultdict
import seaborn as sns
from generate_objects import *
from sklearn.preprocessing import MultiLabelBinarizer
import matplotlib.gridspec as gridspec
from permuation_herbs import Result_distance, roc_cal, prc_cal
from example.casestudy2 import ingre_target_network
from matplotlib.pyplot import figure, draw
from pathlib import Path

logger = logging.getLogger(__name__)

def stylize_axes(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.xaxis.set_tick_params(top='off', direction='out', width=1)
    ax.yaxis.set_tick_params(right='off', direction='out', width=1)


def flatten_list_yin(df, col):
    b_flat = pd.DataFrame([[i, x]
                           for i, y in df[col].apply(list).iteritems()
                           for x in y], columns=['I', col + '_new'])
    b_flat = b_flat.set_index('I')
    df_new = df.merge(b_flat, left_index=True, right_index=True)
    return df_new


def prepared_data_analysis(p_value_auc):
    p_value_auc = p_value_auc.replace('ingredients_as_set', 'ingredient')
    p_value_auc = p_value_auc.replace('target_set_together', 'target')
    p_value_auc['p_auc_list'] = p_value_auc['p_auc_list'].apply(eval)
    p_value_auc['p_values_list'] = p_value_auc['p_values_list'].apply(eval)
    p_value_auc_new = flatten_list_yin(p_value_auc, 'p_auc_list')
    p_value_auc_new_2 = flatten_list_yin(p_value_auc, 'p_values_list')
    p_value_auc_new_3 = p_value_auc_new.append(p_value_auc_new_2)
    return p_value_auc_new_3


# use seaborn
import seaborn as sns
def plot_box_p_auc_sea(p_value_auc_new_3):
    for col_value in ["p_auc_list_new", "p_value_list_new"]:
        plt.figure(figsize=(27, 15))
        sns.set_context("paper", font_scale=2)
        g = sns.boxplot(x="herb_distance_method",
                        y=col_value,
                        hue="ingredient_distance_method",
                        row='filter_score',
                        data=p_value_auc_new_3,
                        kind="box",
                        height=4,
                        aspect=4,
                        legend=False)
        plt.legend(loc='upper left', ncol=3)
        plt.savefig('five_methods/{}_filter_sea.png'.format(col_value))


# figure 2, top 10 frequency herbs
def plot_cor_top_herb(fangji: classmethod, number_herb, out_type):
    '''out type 1: save_figure,
    type 2:'plot_figure
    type 3: only data'''
    herb_top = [h[0][0] for h in fangji.herb_frequency_dict.most_common(number_herb)]
    herbpair_dict = dict(fangji.herb_pair_frequency_dict)
    kept_herbpair = {i: herbpair_dict[i] for i in list(itertools.permutations(herb_top, 2)) if
                     i in list(herbpair_dict.keys())}

    # sns.heatmap(kept_herbpair)
    pdfre = pd.DataFrame(np.nan, columns=herb_top, index=herb_top)
    for k, v in kept_herbpair.items():
        pdfre.loc[k[0], k[1]] = v
        pdfre.loc[k[1], k[0]] = v

    pdfre = pdfre.fillna(0)
    if out_type in ['save_figure', 'plot_figure']:

        fig, ax = plt.subplots(figsize=(12, 12))
        sns.heatmap(pdfre)

        plt.title('coexist frequency of top {} herb'.format(number_herb), fontsize=20)
        if out_type == 'save_figure':
            plt.savefig('figure/cor_herb_top{}'.format(number_herb))
        elif out_type == 'plot_figure':
            plt.show()
    elif out_type == 'only_data':
        return pdfre


def plot_herb_histom(fangji: classmethod, number_herb, out_type):
    herb_top = [h[0][0] for h in fangji.herb_frequency_dict.most_common(number_herb)]
    frequency_top = [h[1] for h in fangji.herb_frequency_dict.most_common(number_herb)]
    if out_type in ['save_figure', 'plot_figure']:

        plt.subplots(figsize=(12, 9))
        ax1 = sns.barplot(herb_top, frequency_top)
        plt.title('Number of herb formulae'.format(number_herb, fontsize=20))
        plt.ylabel('Frequency of herbs')
        ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45)
        plt.tight_layout()
        if out_type == 'save_figure':
            plt.savefig('figure/bar_herb_top{}'.format(number_herb))
        elif out_type == 'plot_figure':
            plt.show()
    elif out_type == 'only_data':
        return herb_top, frequency_top


def plot_herb_combined(fangji: classmethod, number_herb, out_type):
    herb_top = [h[0][0] for h in fangji.herb_frequency_dict.most_common(number_herb)]
    herb_pairs = list(dict(fangji.herb_pair_frequency_dict).keys())
    herb_pairs = [list(h_p) for h_p in herb_pairs if h_p[0] in herb_top or h_p[1] in herb_top]
    herb_pair_merge = reduce(lambda x, y: x + y, herb_pairs)
    herb_count = dict(Counter(herb_pair_merge))

    herb_uni_combine = [herb_count[h] for h in herb_top]
    logger.info('the unique combined of top {} herb are {}', herb_top, herb_uni_combine)
    if out_type in ['save_figure', 'plot_figure']:

        fig, ax = plt.subplots(figsize=(9, 12))
        sns.barplot(herb_uni_combine, herb_top)
        plt.title('count of unique combined herbs for top {} herbs'.format(number_herb), fontsize=20,  y=-0.1)
        plt.ylabel('count of unique combined herbs')
        plt.tight_layout()
        if out_type == 'save_figure':
            plt.savefig('figure/herb_combined_top{}'.format(number_herb))
        elif out_type == 'plot_figure':
            plt.show()
    elif out_type == 'only_data':
        return herb_top, herb_uni_combine

def show_detail_one_pair(herb_info, herb_obj, max_ingre):
    herb_name_1 = herb_info.pinyin_herbid_dic[max_ingre[0]]
    herb_name_2 = herb_info.pinyin_herbid_dic[max_ingre[1]]
    h_1_i = herb_obj.herb_ingre_dict[max_ingre[0]]
    h_2_i = herb_obj.herb_ingre_dict[max_ingre[1]]
    overlap = [i for i in h_1_i if i in h_2_i]
    overlap_ingre_names = [ingredients_obj.ingre_id_name_dict_value[i] for i in overlap]
    return max_ingre[0],herb_name_1, max_ingre[1], herb_name_2, overlap_ingre_names

def herb_overlap_ingredient(herb_pairs_list, herb_obj: classmethod):
    over_list_dict = defaultdict()
    for h_p in herb_pairs_list:
        if h_p[0] in herb_obj.herb_ingre_dict.keys() and h_p[1] in herb_obj.herb_ingre_dict.keys():
            h_1_i = herb_obj.herb_ingre_dict[h_p[0]]
            h_2_i = herb_obj.herb_ingre_dict[h_p[1]]
            len_overlap = len([i for i in h_1_i if i in h_2_i])
            over_list_dict[h_p] = len_overlap
    over_count = dict(Counter(list(over_list_dict.values())))
    return over_count, over_list_dict



def get_ingre_target(ingredients_obj):
    mlb = MultiLabelBinarizer()

    df = pd.DataFrame(mlb.fit_transform(list(ingredients_obj.ingre_tar_dict.values())), columns=mlb.classes_,    index=list(ingredients_obj.ingre_tar_dict.keys()))
    cor_matrix = pd.DataFrame(np.corrcoef(df), columns= list(ingredients_obj.ingre_tar_dict.keys()), index=list(ingredients_obj.ingre_tar_dict.keys()))
    # keep all corr larger than 0.7
    cor_matrix_two_col = cor_matrix.unstack().reset_index()
    cor_matrix_two_col.columns = ['items_1', 'items_2', 'cor_value']
    cor_matrix_two_col = cor_matrix_two_col[cor_matrix_two_col['items_1'] != cor_matrix_two_col['items_2']]

    cor_matrix_two_col = cor_matrix_two_col.sort_values(by=['cor_value','items_1','items_2'])
    cor_matrix_two_col_selected = cor_matrix_two_col[cor_matrix_two_col['cor_value']>0.7]

    cor_matrix_two_col_selected.to_csv('result/cor_ingre_target.txt', sep = '\t')
    return cor_matrix, cor_matrix_two_col_selected

def show_ingre_target(ingre_1, ingre_2, ingredients_obj):
    ingre_name_1 = ingredients_obj.ingre_id_name_dict_value[ingre_1]
    ingre_name_2 = ingredients_obj.ingre_id_name_dict_value[ingre_2]
    target_1 = [ingredients_obj.entry_ensymbl_dict[t] for t in ingredients_obj.ingre_tar_dict[ingre_1]]
    target_2 = [ingredients_obj.entry_ensymbl_dict[t] for t in ingredients_obj.ingre_tar_dict[ingre_2]]
    target_overlap = [t for t in target_2 if t in target_1]
    return ingre_1, ingre_name_1, target_1, ingre_2, ingre_name_2, target_2, target_overlap

def huangqi_gancao(herb_info, herb_obj, max_ingre, cor_matrix ):
    herb_name_1 = herb_info.pinyin_herbid_dic[max_ingre[0]]
    herb_name_2 = herb_info.pinyin_herbid_dic[max_ingre[1]]
    h_1_i = herb_obj.herb_ingre_dict[max_ingre[0]]
    h_2_i = herb_obj.herb_ingre_dict[max_ingre[1]]
    cor_matrix_2 = cor_matrix.loc[h_1_i, h_2_i]
    cor_matrix_two_col = cor_matrix.unstack().reset_index()
    cor_matrix_two_col.columns = ['items_1', 'items_2', 'cor_value']
    cor_matrix_two_col = cor_matrix_two_col[cor_matrix_two_col['items_1'] != cor_matrix_two_col['items_2']]

    cor_matrix_two_col = cor_matrix_two_col.sort_values(by=['cor_value', 'items_1', 'items_2'])
    cor_matrix_two_col_selected = cor_matrix_two_col[cor_matrix_two_col['cor_value'] > 0.7]
    for i, r in cor_matrix_two_col_selected.iterrows():
        print(show_ingre_target(r['items_1'], r['items_2'], ingredients_obj))
    return cor_matrix.loc[h_1_i, h_2_i]

def plot_figure_2(fangji: classmethod, number_herb, out_type):
    # constrained_layout=True
    fig3 = plt.figure(figsize=(8, 8))
    sns.set_context("paper", font_scale=1.2)

    gs = gridspec.GridSpec(5, 5)
    # top herb cor frequency
    f3_ax1 = fig3.add_subplot(gs[1:, 1:])
    pdfre = plot_cor_top_herb(fangji,
                              number_herb, 'only_data')
    pdfre = pdfre.astype(int)
    mask = np.triu(np.ones_like(pdfre, dtype=np.bool))
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    sns.heatmap(pdfre,
                annot= True,
                mask=mask,
                linewidths=.5,
                yticklabels=False,
                xticklabels=False,
                cmap=cmap,
            square=True,
                cbar_kws=dict(use_gridspec=False,
                              shrink=.7,
                              location='right'),
                cbar=True,
                fmt="d")
    f3_ax1.margins(x=1, y=1)
    # herb frequency
    f3_ax2 = fig3.add_subplot(gs[0, 1:])
    herb_top, frequency_top = plot_herb_histom(fangji, number_herb, 'only_data')
    sns.barplot(herb_top, frequency_top)
    # for i in zip(herb_top, frequency_top):
    #     f3_ax2.text(i[0], i[1], i[1], color='black', ha="center")
    f3_ax2.set_xticklabels(f3_ax2.get_xticklabels(), rotation=45)
    f3_ax2.set_title('Number of herb formulae'.format(number_herb))

    # herb quique combined
    f3_ax3 = fig3.add_subplot(gs[1:, 0])  # herb quique combined

    herb_top, herb_uni_combine = plot_herb_combined(fangji, number_herb, 'only_data')
    sns.barplot(herb_uni_combine, herb_top)

    f3_ax3.set_title('Number of coexistent herbs'.format(number_herb), y=-0.1, fontsize = 12)
    #f3_ax3.margins(x=0.01, y=0.01)
    plt.tight_layout()
    # need change by plot out , change right , left , space
    if out_type == 'save_figure':
        plt.savefig('figure/Figure 2.png'.format(number_herb), dpi = 300)
        plt.savefig('figure/Figure 2.pdf', format = 'pdf')
    elif out_type == 'plot_figure':
        plt.show()


# plot fangji distribution

def plot_S1_fangji_length(fangji, out_type):
    fang_len = dict(Counter([len(v) for k, v in fangji.fangji_herb_dict.items()]))
    x_1 = [str(k) for k in list(fang_len.keys())]
    y_1 = list(fang_len.values())
    fang_total = len(fangji.fangji_herb_dict)
    y_2 = [ y * 100 / fang_total for y in y_1]

    herb_mean = np.mean([len(v) for k, v in fangji.fangji_herbid_dict.items()])
    fang_len_dict = {k:len(v) for k, v in fangji.fangji_herb_dict.items()}
    # max_fangji =  max(fang_len_dict.items(), key=operator.itemgetter(1))[0]
    sns.set_style("white")
    fig, ax1 = plt.subplots(figsize=(8, 6))
    sns.set_context("paper", font_scale=1.2)
    plt.xlabel('Number of herbs')
    ax1.set_ylabel('Count of TCM formulae',
                   color="black")
    ax1 = sns.barplot(x_1, y_1, order=x_1)
    for label in ax1.get_xticklabels():
        label.set_visible(False)
    for label in ax1.get_xticklabels()[::5]:
        label.set_visible(True)
    ax1.tick_params(axis='y')
    ax2 = ax1.twinx()
    ax2.set_ylabel('Percentage %',
                   color="black")
    ax2= sns.lineplot(x_1,
             y_2,
             color="black",
                      marker="o",
                      lw=0.5,markersize=6,
                      sort=False
             )

    #plt.title('The number of herb in one prescription', fontsize=20)
    if out_type == 'save_figure':
        plt.savefig('figure/Figure S1.pdf',   format = 'pdf')
        plt.savefig('figure/Figure S1.png', dpi=300)
    elif out_type == 'plot_figure':
        plt.show()


def plot_S2_fangi_frequency(fangji, values_200_pd, out_type):
    values_200 = values_200_pd.drop_duplicates(subset=['herb1',
                                                       'herb1_name',
                                                       'herb2',
                                                       'herb2_name',
                                                       'frequency'], keep='first')['frequency']
    values = list(dict(fangji.herb_pair_frequency_dict).values())
    figS2 = plt.figure(figsize=(8, 6))
    sns.set_context("paper", font_scale=1.2)
    gs = gridspec.GridSpec(2, 2)
    # top herb cor frequency
    sns.set_style("white")
    f3_ax1 = figS2.add_subplot(gs[0:, 0:])
    plt.hist(values,
                  log=True,
                  bins=2 ** np.arange(12),
                  color='brown',
                  lw=2, alpha=0.8,
                  rwidth=0.90
                  )
    plt.plot([0,1500], [338,338],
             linestyle = '--',
             marker = 'o',
             markersize=3,
             color = 'black')

    plt.xscale('log')
    plt.xlabel('Frequency')
    plt.ylabel('Number of herb pair')
    f3_ax1.spines['top'].set_visible(False)
    f3_ax1.spines['right'].set_visible(False)

    f3_ax1.xaxis.set_tick_params(top='off', direction='out', width=1)
    f3_ax1.yaxis.set_tick_params(right='off', direction='out', width=1)

    f3_ax2 = figS2.add_subplot(gs[0:1, 1:])
    plt.hist(values_200,
             bins=10,
             color='darkgreen',
             lw=2,
             alpha=0.8,
             rwidth=0.90
             )
    f3_ax2.spines['top'].set_visible(False)
    f3_ax2.spines['right'].set_visible(False)
    f3_ax2.xaxis.set_tick_params(top='off', direction='out', width=1)
    f3_ax2.yaxis.set_tick_params(right='off', direction='out', width=1)
    plt.xlabel('Frequency')
    plt.ylabel('Number of herb pairs')
    plt.tight_layout()

    if out_type == 'save_figure':
        plt.savefig('figure/Figure S2.pdf', format='pdf')
        plt.savefig('figure/Figure S2.png', dpi=300)
    elif out_type == 'plot_figure':
        plt.show()

# ingredients oberlap
def plot_S3_ingredient_overlap( top_pd, random_pd, herb_obj: classmethod, out_type):
    pairs_random = list(set(list(zip(random_pd['herb1'],random_pd['herb2']))))
    top_pd = top_pd.sort_values(by='frequency', ascending=False)
    pairs_top = list(set(list(zip(top_pd['herb1'],top_pd['herb2']))))

    overlap_count_top, over_list_dict_top = herb_overlap_ingredient(pairs_top, herb_obj)
    overlap_count_random, over_list_dict_random = herb_overlap_ingredient(pairs_random, herb_obj)
    max_ingre = max(over_list_dict_top.items(), key=operator.itemgetter(1))
    #print(show_detail_one_pair(herb_info, herb_obj, max_ingre[0]))
    #logger.info('the max overlap {}',  max_ingre)
    #logger.info('the number of ingredients overlap {} between {}', overlap_count_top, over_list_dict_top)
    figS2 = plt.figure(figsize=(8, 8))
    gs = gridspec.GridSpec(2, 1)
    sns.set_style("white")
    sns.set_context("paper", font_scale=1.2)
    # top herb cor frequency
    f3_ax1 = figS2.add_subplot(gs[0:1, :])
    top_dict = dict(Counter(list(over_list_dict_top.values())))
    #sns.distplot(list(over_list_dict_top.values()), norm_hist = False)
    pd_data = pd.DataFrame.from_dict(top_dict, orient='index').reset_index()
    pd_data.columns = ['group', 'count']
    #pd_data['group'] = pd_data['group'] .astype(str)
    pd_data = pd_data.sort_values(by='group', ascending=True)
    pd_data['count_per'] = pd_data['count'].apply(lambda x:x/sum(pd_data['count'])*100)
    sns.barplot(x='group',
                y='count_per',
                data=pd_data)

    f3_ax1.set_ylim(0, pd_data['count_per'].max() + 15)
    for p in f3_ax1.patches:
        f3_ax1.annotate(format(p.get_height(), '.2f'),
                        (p.get_x() + p.get_width() / 2.+0.1, p.get_height()),
                        ha='center',
                       va='center',
                        xytext=(0, 10),
                        textcoords='offset points',
                        )
    plt.ylabel('Percentage %')
    plt.xlabel('')
    f3_ax2 = figS2.add_subplot(gs[1:, :])
    random_dict = dict(Counter(list(over_list_dict_random.values())))
    pd_data_r = pd.DataFrame.from_dict(random_dict, orient='index').reset_index()
    pd_data_r.columns = ['group', 'count']
    pd_data_r = pd_data_r.sort_values(by = 'group', ascending=True)
    pd_data_r['count_per'] = pd_data_r['count'].apply(lambda x: x / sum(pd_data_r['count']) * 100)
    sns.barplot(x='group',
                y='count_per',
                data=pd_data_r)
    f3_ax2.set_ylim(0, pd_data_r['count_per'].max() + 15)
    for p in f3_ax2.patches:
        f3_ax2.annotate(format(p.get_height(), '.2f'), (p.get_x() + p.get_width() / 2., p.get_height()), ha='center',
                       va='center', xytext=(0, 10), textcoords='offset points')
    plt.ylabel('Percentage %')
    plt.xlabel('Number of common ingredients')
    if out_type == 'save_figure':
        plt.savefig('figure/Figure S3.png', dpi = 300)
        plt.savefig('figure/Figure S3.pdf', format='pdf')
    elif out_type == 'plot_figure':
        plt.show()

def plot_fig_3(mean_pd, out_type):
    mean_pd['distance_reduce'] = mean_pd['Distance for random herb pairs'] - mean_pd['Distance for top herb pairs']
    mean_pd_matrix = mean_pd.pivot_table(columns=['Herb-level distance type'],index=['Ingredient-level distance type'], values=['distance_reduce'])
    mean_pd_melt = mean_pd.melt(id_vars=['Ingredient-level distance type',	'Herb-level distance type'
], value_vars=['Distance for top herb pairs',	'Distance for random herb pairs'])
    mean_pd_melt = mean_pd_melt.rename(columns = {'value':'Distance', 'variable':'Source'})
    fig3 = plt.figure(figsize=(5, 5))
    sns.set_context("paper", font_scale=1)
    gs = gridspec.GridSpec(3, 3)
    ax_3 = fig3.add_subplot(gs[:, :])
    sns.factorplot(x='Ingredient-level distance type',
                   y='Distance',
                   hue='Source',
                   col='Herb-level distance type',
                       col_wrap=3,
                   data = mean_pd_melt,
                   kind='bar', ci=None, legend=False, palette='Paired' )
    plt.legend(loc='best')
    ax_4 = fig3.add_subplot(gs[1:, 2:])
    sns.heatmap(mean_pd_matrix, annot=True, ax=ax_4)
    plt.tight_layout()
    if out_type == 'save_figure':
        plt.savefig('figure/Figure 3.png')
    elif out_type == 'plot_figure':
        plt.show()

def new_figue_3(out_type):
    filename_top = 'result/top_200.csv'
    filename_random = 'result/result_random_0_10000.csv'
    result = Result_distance(filename_top, filename_random)  # 47 time
    file_folder = 'top_random_new'
    result.get_mean()

    fig = plt.figure(figsize=(10, 8))
    pd_melt = result.pd_melt
    pd_melt = pd_melt.rename(columns={'Herb-level distance type':'H','Ingredient-level distance type':'I','class':'group'})
    sns.set_context('paper', font_scale=2)
    sns.set_style('white')
    g = sns.FacetGrid(pd_melt,
                      row='H',
                      col= 'I',
                      hue_kws={"alpha": [0.25, 0.25]},
                      margin_titles = True,
                      palette="Set2",
                      despine = False)
    g = g.map(sns.violinplot,
                      'group',
              "Distance",
              palette="Set2")
    g.add_legend()
    plt.tight_layout()
    if out_type == 'save_figure':
        plt.savefig('figure/Figure 3.png', dpi = 300)
        plt.savefig('figure/Figure 3.pdf', format='pdf')
    elif out_type == 'plot_figure':
        plt.show()

    if out_type == 'save':
        plt.savefig('result/{}/density.png'.format(file_folder))
    elif out_type == 'show':
        plt.show()


def plot_fig_S5( mean_pd, out_type):
    mean_pd['distance_reduce'] = mean_pd['Distance for random herb pairs'] - mean_pd['Distance for top herb pairs']
    mean_pd_matrix = mean_pd.pivot_table(columns=['Herb-level distance type'], index=['Ingredient-level distance type'],
                                         values=['distance_reduce'])
    mean_pd_melt = mean_pd.melt(id_vars=['Ingredient-level distance type', 'Herb-level distance type'
                                         ],
                                value_vars=['Distance for top herb pairs', 'Distance for random herb pairs'])
    mean_pd_melt = mean_pd_melt.rename(columns={'value': 'Distance', 'variable': 'Source'})
    fig3 = plt.figure(figsize=(5, 5))
    sns.set_context("paper", font_scale=1)
    gs = gridspec.GridSpec(3, 3)
    ax_3 = fig3.add_subplot(gs[:, :])
    sns.factorplot(x='Ingredient-level distance type',
                   y='Distance',
                   hue='Source',
                   col='Herb-level distance type',
                   col_wrap=3,
                   data=mean_pd_melt,
                   kind='bar',
                   ci=None,
                   legend=False,
                   palette='Paired' )
    plt.legend(loc='best')
    ax_4 = fig3.add_subplot(gs[1:, 2:])
    sns.heatmap(mean_pd_matrix, annot=True, ax=ax_4)
    plt.tight_layout()
    if out_type == 'save_figure':
        plt.savefig('figure/Figure S5.png')
    elif out_type == 'plot_figure':
        plt.show()

def plot_figure_4(filename_top, filename_random, herb_m, ingre_m, how):
    result = Result_distance(filename_top, filename_random)
    result.distance_split()
    repeat = len(list(result.new_grouped_list['random'])[1]) - 1
    fig = plt.figure(figsize=(10, 5))
    sns.set_context('paper', font_scale=1)
    ax = fig.add_subplot(1, 2, 1)
    top = list(result.new_grouped_list.loc[(herb_m, ingre_m), 'top'])
    randoms = list(result.new_grouped_list.loc[(herb_m, ingre_m), 'random'])
    auc_list = []; fpr_list=[];tpr_list=[]
    for i, random_list in enumerate(randoms):
        if i < repeat:
            out = roc_cal(top, random_list)
            fpr, tpr, roc_auc = out['fpr'], out['tpr'], out['roc']
            fpr_list.append(fpr); tpr_list.append(tpr); auc_list.append(roc_auc)
            ax.plot(fpr, tpr, lw=2, alpha=0.3)
    ax.plot(np.mean(pd.DataFrame(fpr_list)), np.mean(pd.DataFrame(tpr_list)), '',
            lw=2, alpha=0.8, color='blue',
            label='AUROC = {:.3f}'.format(np.mean(auc_list)))
    ax.plot([0, 1], [0, 1], '--')
    ax.set_ylabel('True Positive Rate'); ax.set_xlabel('False Positive Rate')
    ax.legend(loc='lower right', fontsize='large')

    # PRC
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.set_ylabel('PRC')
    ax2.set_alpha(.6)
    auc_list = []
    precision_list = []
    recall_list = []
    for i, random_list in enumerate(randoms):
        if i < repeat:
            out = prc_cal(top, random_list)
            precision, recall, prc_auc = out['precision'], out['recall'], out['roc']
            precision_list.append(precision)
            recall_list.append(recall)
            auc_list.append(prc_auc)
            ax2.plot(recall, precision, lw=2, alpha=0.3)
    ax2.plot(np.mean(pd.DataFrame(recall_list)), np.mean(pd.DataFrame(precision_list)), '', lw=2, alpha=0.8,
             color='blue',
             label='AUPRC = {:.3f}'.format(np.mean(auc_list)))
    ax2.set_ylabel('Precision')
    ax2.set_xlabel('Recall')
    ax2.legend(loc='lower right', fontsize='large')
    plt.subplots_adjust(hspace=1.2, wspace=0.50)
    plt.tight_layout()

    plt.subplots_adjust(hspace=1.2, wspace=0.50); plt.tight_layout()
    if how == 'save_figure':
        plt.savefig('figure/Figure 4.png')
    elif how == 'plot_figure':
        plt.show()


def plot_figure_5(mean_pd, how):
    mean_pd_melt = mean_pd.melt(id_vars=['Herb-level distance type',
                       'Ingredient-level distance type'], value_vars=['AUROC',
                       'AUPRC'])
    ax = plt.figure(figsize=(5, 5))
    sns.set_context('paper', font_scale=1)
    g = sns.barplot(x= 'Herb-level distance type',
                        y='value',
                        hue="variable",
                        data=mean_pd_melt,
                    palette="Set1",
                    ci=80,
                    errcolor='black',
                    errwidth=1.5
                    )
    plt.legend(title = '')
    plt.ylabel('')
    plt.xlabel('')
    if how == 'save_figure':
        plt.savefig('figure/Figure 5.png')
    elif how == 'plot_figure':
        plt.show()



def plot_figure_S6(mean_pd, how):
    #mean_pd_melt = mean_pd.melt(id_vars=['Herb-level distance type',
                                         #'Ingredient-level distance type'], #value_vars=['AUROC'])
    ax = plt.figure(figsize=((5, 5)))
    sns.set_context('paper', font_scale=1)
    g = sns.barplot(x='Herb-level distance type',
                    y='AUROC',
                    ci=95,
                    data=mean_pd,
                    errcolor='black',
                    errwidth=1.5,
                    palette='Dark2'

                    )
    plt.legend(title='')
    plt.ylabel('AUROC')
    plt.xlabel('')
    if how == 'save_figure':
        plt.savefig('figure/Figure S6.png')
    elif how == 'plot_figure':
        plt.show()

def plot_figure_6(g_obj, ingredients_obj, how):
    ingre_1 = ['I23134']
    ingre_2 = ['I13091']
    ingre_target_network(ingre_1, ingre_2, g_obj, ingredients_obj, how)


def plot_box_diatance_top_no_overlap():
    mean_pd_top = pd.read_csv('result/top_random/mean.csv')
    mean_pd_top['source'] = 'top'
    mean_pd_no_over = pd.read_csv('result/recom_random/mean.csv')
    mean_pd_no_over['source'] = 'top_no_overlap'
    mean_pd = pd.concat([mean_pd_top,mean_pd_no_over])
    ax = plt.figure()
    sns.set_context("paper")
    my_pal = {"top": "r", "top_no_overlap": "g"}
    g = sns.boxplot(x='herb_method',
                        y='top',
                        hue="source",
                        data=mean_pd, whis=200, palette=my_pal)
    g.set_title('The distance of top and top without overlaped ingreidnets herb pairs')
    g.set_xlabel('method')
    g.set_ylabel('distance_mean_value')
    ax.savefig('figure/top_no_overlap_dis_mean.png')

def main():
    # plot_fig_3('plot_figure')
    # plot_figure_4('center', 'shortest', 'save_figure')
    # plot_figure_5()
    # plot_figure_S6('save_figure')
    # plot_figure_2(fangji, 10, 'save_figure')
    # plot_S1_fangji_length(fangji,'save_figure')
    # plot_S2_fangi_frequency(fangji, 'save_figure')
    plot_figure_6(g_obj, ingredients_obj, 'save_figure')


