from figures_draw import *

from generate_objects import herb_obj, herb_info, fangji



#test the functions
def test_herb_ingre(herb_distance_obj, g_obj):
    herb_distance_obj.Ingredients.ingre_ingre_dis('I1981', 'I33248', g_obj.G, ['separation'])
    herb_distance_obj.Ingredients.ingre_ingre_dis_all('I1981', 'I33248', g_obj.G)
    #
    herb_distance_obj.herb_herb_dis('H10', 'H1010', 'closest',  ['closest'])
    herb_distance_obj.herb_herb_distance_uni('H10', 'H1010', 'closest')
    s = herb_distance_obj.herb_herb_dis_all('H10', 'H1010')


def main():
    test_figure_obj = test_figure(fangji, herb_obj, herb_info)
    test_figure_obj.test_herb_ingredient()




