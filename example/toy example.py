from construct_network import *
from herb_ingre_tar import *
from proximity_key import *
import networkx as nx

# construct network
g_obj = Construct_Network("example/toy3.sif")

import matplotlib.pyplot as plt

options = {
    'node_color': 'grey',
    'node_size': 400,
    'width': 3
}
plt.subplot(121)
nx.draw(g_obj.G, with_labels=True, font_weight='bold',**options)
plt.show()


# example drug_drug [dAB, closest, shortest, kernel, center]

ingre_targets_dict = {'network1': ['Ta1', 'Ta2', 'Ta3', 'Ta4','Ta5','Ta6'],
                      'network2': ['Tb1', 'Tb2', 'Tb3', 'Tb4', 'Tb5'],
                      'network3': ['Ta5', 'Ta6'],
                      'network4': ['Tb2', 'Tb4', 'Tb5'],
                      'network5': ['Tb1', 'Tb3', 'Tb4', 'Ta4']}

def distances(nodes_from, nodes_to):

    length_dict = Sets_Lengths(nodes_from, nodes_to).target_lengths(g_obj.G)

    d = Network_Distance(nodes_from, nodes_to, length_dict)

    method_list = ['separation', 'closest', 'shortest', 'kernel', 'center']

    return [d.network_distance(method) for method in method_list]

nodes_from = ingre_targets_dict['network1']
nodes_to = ingre_targets_dict['network2']

distances(nodes_from, nodes_to)

nodes_from = ingre_targets_dict['network4']
nodes_to = ingre_targets_dict['network3']

distances(nodes_from, nodes_to)

# herb_herb, one thread
herb_ingre_dict = {'X': ['network1', 'network2'],
                   'Y': ['network3', 'network4', 'network5']}

nodes_from = herb_ingre_dict['X']
nodes_to = herb_ingre_dict['Y']

def distance_ingre(nodes_from, nodes_to, method_ingre):

    real_nodes_from = ingre_targets_dict[nodes_from]
    real_nodes_to = ingre_targets_dict[nodes_to]
    length_dict = Sets_Lengths(real_nodes_from, real_nodes_to).target_lengths(g_obj.G)

    d = Network_Distance(real_nodes_from, real_nodes_to, length_dict)

    return d.network_distance(method_ingre)


length_dict = Sets_Lengths(nodes_from, nodes_to).ingre_length(distance_ingre, 'closest')

d = Network_Distance(nodes_from, nodes_to, length_dict)

d.network_distance('closest')


# parallel

def herb_herb_distance_test(nodes_from, nodes_to):
    method_list_ingre = ['separation', 'closest', 'shortest', 'kernel', 'center']
    method_list_herb = ['separation', 'closest', 'shortest', 'kernel', 'center']
    dis_dict = defaultdict()
    for m in method_list_ingre:
        length_dict = Sets_Lengths(nodes_from, nodes_to).ingre_length(distance_ingre, m)
        d = Network_Distance(nodes_from, nodes_to, length_dict)

        dis_dict_ingre = defaultdict()
        for n in method_list_herb:
            dis_dict_ingre[n] = d.network_distance(n)

        dis_dict[m] = dis_dict_ingre
    return dis_dict

nodes_from = herb_ingre_dict['X']
nodes_to = herb_ingre_dict['Y']
herb_herb_distance_test(nodes_from, nodes_to)
