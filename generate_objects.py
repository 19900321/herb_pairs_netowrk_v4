
from herb_ingre_tar import *
from construct_network import *
from herb_distance_generation import *

# generate object

g_obj = Construct_Network("source_data/toy.sif")
ingredients_obj = Ingredients('source_data/stitch_com_entry_ids.csv', 0)
herb_obj = Herb('source_data/herb_ingredient_pairs.csv')
herb_distance_obj = Herb_Distance(g_obj, ingredients_obj, herb_obj)


from herb_herb_pairs import *
herb_info = Herb_Info('source_data/herb_info.txt')
fangji = FangJi('source_data/prescription.txt')
fangji.herbid_frequency_dic(2, herb_info.herb_pinyin_dic)