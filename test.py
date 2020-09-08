from figures_draw import *

from generate_objects import herb_obj, herb_info, fangji

class test_figure:
    def __init__(self, fangji, herb_obj, herb_info):
        self.fangji = fangji
        self.herb_obj = herb_obj
        self.herb_info = herb_info


    def test_herb_cor(self):
        plot_cor_top_herb(self.fangji, 10, 'plot_figure')

    def test_figure2(self):
        plot_figure_2(self.fangji, 10, 'save_figure')

    def test_fangji_length(self):
        plot_fangji_length(self.fangji, 'save_figure')

    def test_herb_ingredient(self):
        plot_ingredient_overlap(self.herb_obj,
                                self.fangji,
                                self.herb_info,
                                200, 'top',
                                'save_figure' )



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
    #test_figure_obj.test_herb_cor()
    # test_figure_obj.test_figure2()
    # test_figure_obj.test_fangji_length()
    test_figure_obj.test_herb_ingredient()




