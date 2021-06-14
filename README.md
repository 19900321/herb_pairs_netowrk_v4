
# This is a network-based method written by Python to calculate the herb-herb ingredients in the PPI network. The structure of the whole project is like below: 


# 1 construct_network.py: 
help to set up protein-protein interaction network
For instance, for a given file path with PPI edgelist pairs e.g. example/toy.sif, Construct_Network object will create by "network = Construct_Network("example/toy.sif")".

# 2 proximity_key.py: 
contains the key script that generates five different network distances, such as center, separation, closest, shortest and kernel.

# 3 herb_ingre_tar.py: 
the object named ‘Ingredients’ contains ingredient-related information, especially ingredient-target relationship. The object ‘Herb’ wraps herb-related information, especially ingredient pairs.

# 4 herb_distance_generation.py: 
calculate the herb-herb network distance by "herb_herb_dis(self, herb_from, herb_to, distance_method, distance_method_herb_list)" and "herb_herb_dis_all(self, herb_from, herb_to)" for all the distance. 

# 5 herb_herb_pairs.py:
contain herb informtaion and TCM formulae information

# 6 generate_objects.py: 
generate key objects that will be used in the methods.

# 7 The example folder contains two examples.
