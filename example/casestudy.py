# good and bad herb pairs from literature
import pandas as pd

# read data

import itertools
def prepare_herbpairs(data_pd, herb_pinyin_dic):
    herb_pairs_list = []
    data_pd['Herb A'] = data_pd['Herb A'].str.upper()
    data_pd['Herb B'] = data_pd['Herb B'].str.upper()
    for i, herb_pair in data_pd.iterrows():
            herb_name1 = herb_pair['Herb A']

            herb_name2 = herb_pair['Herb B']
            frequency = herb_pair['group']
            try:
                herb_id1 = herb_pinyin_dic[herb_name1]
                herb_id2 = herb_pinyin_dic[herb_name2]
                herb_pairs_list += [(herb_id1, herb_name1, herb_id2, herb_name2, frequency)]
            except:
                continue
    return herb_pairs_list

# here to reuse the functions , we treat the known as top, the forbid
def standarize_distance(classic_distances):
    classic_distances['class'] = classic_distances['frequency']
    classic_distances_back = classic_distances[classic_distances.columns]
    classic_distances_back = classic_distances_back.drop(['class'], axis=1)
    classic_distances_back['frequency'] = classic_distances_back['frequency'].str.replace('recommend', '10000')
    classic_distances_back['frequency'] = classic_distances_back['frequency'].str.replace('forbid', '-10000')
    classic_distances_back['frequency'] = classic_distances_back['frequency'].astype(int)
    return classic_distances_back

# read as result objects
# random and known herbs
from permuation_herbs import *

os.chdir('.')
# !! DO NOT NEED TO RUN UNLESS GENERATE NEW RESULT
# classic_pairs = pd.read_csv('source_data/herb_pairs_classic.csv')
# classic_herb_pairs = prepare_herbpairs(classic_pairs, herb_info.herb_pinyin_dic)
# classic_distances = herb_distance_obj.generator_result(classic_herb_pairs)
# classic_distances_2 = standarize_distance(classic_distances)
# classic_distances_2[classic_distances_2['frequency'] == 10000].to_csv('result/classic_distances_recommend.csv')
# classic_distances_2[classic_distances_2['frequency'] == -10000].to_csv('result/classic_distances_forbid.csv')

# here to reuse the functions , we treat the known as top, the forbid
filename_top = 'result/classic_distances_recommend.csv'
filename_forbid = 'result/classic_distances_forbid.csv'
filename_random = 'result/result_random_0_10000.csv'

result_recom_forbid = Result_distance(filename_top, filename_forbid)
set_path = 'result/recom_forbid/'
result_recom_forbid.save_analysis(-10001, 12, 1, set_path, forbid_include=True)

# use recommend and random
result_recom_random = Result_distance(filename_top, filename_random)
set_path = 'result/recom_random/'
result_recom_random.save_analysis(0, 268, 35, set_path, forbid_include=False)


