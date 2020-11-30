from proximity_key import *

class Disease:

    def __init__(self, disease_file_name, network):
        self.network = network
        self.disease_file_name = disease_file_name
        self.data = self.read_data()
        self.disease_name = self.data['Disease Name'].unique()
        self.disease_tar_dict_all = self.get_disease_target_dict_all()
        self.disease_tar_dict = self.get_disease_target_dict()

    def read_data(self):
        data = pd.read_csv(self.disease_file_name)
        data['Gene ID'] = 'T' + data['Gene ID'].astype(str)
        return data

    def get_disease_target_dict_all(self):
        dis_tar_dict = dict(self.data.groupby(['Disease Name'])['Gene ID'].apply(list))
        return dis_tar_dict

    def get_disease_target_dict(self):
        dict_filter = {k: [v2 for v2 in v if v2 in self.network.G.nodes]
                      for k, v in self.disease_tar_dict_all.items()}
        dict_filter = {k: v for k, v in dict_filter.items() if len(v) !=0}
        return dict_filter

    def update_disease_target_dict(self, new_dict):
        dict_filter = {k: [v2 for v2 in v if v2 in self.network.G.nodes]
                       for k, v in new_dict.items()}
        dict_filter = {k: v for k, v in dict_filter.items() if len(v) != 0}
        self.disease_tar_dict.update(dict_filter)


    def cal_distance_ob(self, nodes_from, nodes_to):
        length_dict = Sets_Lengths(nodes_from, nodes_to).target_lengths(self.network.G)
        dis_obj = Network_Distance(nodes_from, nodes_to, length_dict)

        return dis_obj

    def cal_disease_ingre_dis(self, disease, ingre, distance_method, ingre_tar_dict):
        if ingre not in ingre_tar_dict.keys():
            print('{} not in ingre_tar_dict dictionary'.format(ingre))
            return None
        elif disease not in self.disease_tar_dict.keys():
            print('{} not in disease_tar_dict dictionary'.format(disease))
            return None
        else:
            nodes_from = self.disease_tar_dict[disease]
            nodes_to = ingre_tar_dict[ingre]
            dis_obj= self.cal_distance_ob(nodes_from, nodes_to)
            distance = dis_obj.network_distance(distance_method)

            return distance

    def cal_disease_ingre_z_score(self, disease, ingre, distance_method, ingre_tar_dict, random_time, seed):
        if ingre not in ingre_tar_dict.keys():
            print('{} not in ingre_tar_dict dictionary'.format(ingre))
            return None
        elif disease not in self.disease_tar_dict.keys():
            print('{} not in disease_tar_dict dictionary'.format(disease))
            return None
        else:
            nodes_from = self.disease_tar_dict[disease]
            nodes_to = ingre_tar_dict[ingre]
            dis_obj = self.cal_distance_ob(nodes_from, nodes_to)
            d, z, (m, s), pval = dis_obj.cal_z_score(distance_method, random_time, self.network, seed)

        return d, z, m, s, pval

    def cal_disease_herb_z_score(self, disease, herb, distance_method, herb_ingretargets_dic, random_time, seed):
        if herb not in herb_ingretargets_dic.keys():
            print('{} not in ingre_tar_dict dictionary'.format(herb))
            return None
        elif disease not in self.disease_tar_dict.keys():
            print('{} not in disease_tar_dict dictionary'.format(disease))
            return None
        else:
            nodes_from = self.disease_tar_dict[disease]
            nodes_to = herb_ingretargets_dic[herb]
            dis_obj = self.cal_distance_ob(nodes_from, nodes_to)
            d, z, (m, s), pval = dis_obj.cal_z_score(distance_method, random_time, self.network, seed)

        return d, z, m, s, pval

    def cal_disease_herb_ingre_z_score(self, disease, herb, distance_method, herb_ingre_dict, ingre_tar_dict, random_time, seed):

        ingres = herb_ingre_dict[herb]
        ingre_z_score_dict = {ingre: self.cal_disease_ingre_z_score(disease,
                                                                    ingre, distance_method, ingre_tar_dict,
                                                                    random_time, seed) for ingre in ingres}
        return ingre_z_score_dict
