import mygene
import pubchempy as pcp
from generate_objects import *
from herb_ingre_tar import *

n = 0
def drug_name_inchikey(drug_name):
    global n
    print(n)
    n +=1
    drug_identider = pcp.get_synonyms(drug_name, 'name')
    if len(drug_identider) != 0:
        drug_identider_dic = drug_identider[0]
        CID = drug_identider_dic['CID']
        c = pcp.Compound.from_cid(CID)
        inchikey = c.inchikey
        return inchikey
    else:
        return 'unknown'

def joint_yin(x):
    if x == None:
        return None
    else:
        return ','.join(list(set(x)))

def transfer_data():
    def expand_list(df, list_column, new_column):
        lens_of_lists = df[list_column].apply(len)
        origin_rows = range(df.shape[0])
        destination_rows = np.repeat(origin_rows, lens_of_lists)
        non_list_cols = (
          [idx for idx, col in enumerate(df.columns)
           if col != list_column]
        )
        expanded_df = df.iloc[destination_rows, non_list_cols].copy()
        expanded_df[new_column] = (
          [item for items in df[list_column] for item in items]
          )
        expanded_df.reset_index(inplace=True, drop=True)
        return expanded_df

    def float_judge(x):
        return isinstance(x, float)

    drug_targte_yin = pd.read_csv('../herb_pairs_materials/case study/drugcom_dss/drug_target_yin.csv')
    drug_targte_yin['symbols'] = drug_targte_yin['symbols'].str.split(',')
    drug_targte_yin = drug_targte_yin[~drug_targte_yin['symbols'].apply(float_judge)]
    drug_targte_yin = expand_list(drug_targte_yin,  'symbols', 'symbol')
    drug_targte_yin.to_csv('../herb_pairs_materials/case study/drugcom_dss/drug_target_yin_expand.csv')


class Drug_Info:

    def __init__(self):
        pass

    def get_drug_info_yin(self, drug_target_info_file):
        self.drug_target_info_file = drug_target_info_file
        self.gene_column = 'symbol'
        self.pd_drug = pd.read_csv(self.drug_target_info_file)

    def get_drug_info_zia(self, drug_info_file, drug_target_file):
        self.gene_column = 'uniprot'
        self.drug_info_file = drug_info_file
        self.drug_info = pd.read_csv(self.drug_info_file, sep='\t')
        self.drug_target_file = drug_target_file
        self.drug_target = pd.read_csv(self.drug_target_file)
        self.drug_target = self.drug_target[self.drug_target['interaction_strength'] > 0.7]
        self.pd_drug = pd.merge(self.drug_info, self.drug_target,
                                how='inner', left_on='inchi_key', right_on='standard_inchi_key')

    def get_gene_name_entryid_dict(self):
        self.mygene = mygene.MyGeneInfo()
        geneSyms = self.mygene.querymany(list(self.pd_drug[self.gene_column].unique()),
                                         scopes=self.gene_column, fields='entrezgene', species='human')
        self.geneSyms = {i['query']: i['entrezgene'] for i in geneSyms if 'entrezgene' in i.keys()}


    def trans_some_gene_name_entryid(self, some_gene_name):
        if str(some_gene_name) in self.geneSyms.keys():
            return self.geneSyms[str(some_gene_name)]
        else:
            return None

    def add_gene_id_column(self):
        self.pd_drug['gene_id'] = self.pd_drug[self.gene_column].astype(str).apply(self.trans_some_gene_name_entryid)

    def get_ingredient_ogj(self):
        drug_save = self.pd_drug.dropna(subset=['gene_id'])
        drug_save = drug_save.groupby(['inchi_key', 'pert_iname'])['gene_id'].apply(list)
        drug_save = drug_save.reset_index()
        drug_save['score'] = '1.0,1.0'
        drug_save['ensymble'] = 'unknown'
        drug_save.columns = ['ingredients_id', 'name', 'targets', 'score', 'ensymble']
        drug_save['targets'] = drug_save['targets'].apply(joint_yin)
        drug_save.to_csv('../herb_pairs_materials/case study/drugcom_dss/ingredients_format.csv', sep='\t')

        self.drug_ingredients_obj = Ingredients('../herb_pairs_materials/case study/drugcom_dss/ingredients_format.csv', 0)

        self.drug_ingredients_obj.ingredients_target_dict(g_obj.G.nodes)



import pandas as pd

drug_info_file = '../herb_pairs_materials/case study/drugcom_dss/GSE70138_compound_list.tsv'
drug_target_file = '../herb_pairs_materials/case study/drugcom_dss/Alberto_compound_targets.csv'
drug_info_ogj = Drug_Info()
drug_info_ogj.get_drug_info_zia(drug_info_file, drug_target_file)
drug_info_ogj.get_gene_name_entryid_dict()
drug_info_ogj.add_gene_id_column()
drug_info_ogj.get_ingredient_ogj()
drug_info_ogj.drug_ingredients_obj.ingre_ingre_dis_all('IALWKGYPQUAPLQC-UHFFFAOYSA-N',
                                                      'IANGKOCUUWGHLCE-UHFFFAOYSA-N', g_obj.G)
# if use the zia excel, we still need

class DrugCom:
    def __init__(self, drug_combin_file):
        self.drug_combin_file = drug_combin_file

    def get_drug_name_inchey_dict(self, drug_inchey_file):
        self.drug_inchey_dict_file = drug_inchey_file
        pd_drug_inchey_dict = pd.read_csv(self.drug_inchey_dict_file)
        pd_drug_inchey_dict = pd_drug_inchey_dict[['drug_name', 'drug_inchey']].drop_duplicates(keep='first')
        self.drug_inchey_dict = dict(zip(pd_drug_inchey_dict['drug_name'], pd_drug_inchey_dict['drug_inchey']))

    def get_targets_drugcom(self, drug_chemble):
        self.drug_combin = pd.read_csv(self.drug_combin_file)

    def trans_drug_name_inchey(self, drug_name):
        if drug_name in self.drug_inchey_dict.keys():
            return self.drug_inchey_dict[drug_name]
        else:
            return None

    def add_inchkey_column(self, data_pd):
        data_pd['drug_row_inchey'] = data_pd['drug_row'].apply(self.trans_drug_name_inchey)
        data_pd['drug_col_inchey'] = data_pd['drug_col'].apply(self.trans_drug_name_inchey)
        return data_pd

    def get_drug_pairs(self):
        self.drug_pairs = list(set(zip(self.drug_combin['drug_row_inchey'],
                                       self.drug_combin['drug_col_inchey'])))

    def cal_drugs_net_distances(self, drug_ingredients_obj, g_obj, drug_1, drug_2):
        return drug_ingredients_obj.ingre_ingre_dis_all(drug_1, drug_2, g_obj.G)

    def cal_drugs_net_distances_dict(self, drug_ingredients_obj, g_obj):
        self.drugs_net_distances_dict = {drugs:self.cal_drugs_net_distances(drug_ingredients_obj, g_obj,
                                                                      list(drugs)[0], list(drugs)[1]) for drugs in self.drug_pairs}

    def add_drugs_net_distances_all(self):
        self.drug_combin['drug_pairs'] = list(zip(self.drug_combin['drug_row_inchey'],
                                                  self.drug_combin['drug_col_inchey']))
        drugs_net_distances = pd.DataFrame.from_dict(self.drugs_net_distances_dict,
                                                     orient='index')
        drugs_net_distances['drugs'] = drugs_net_distances.index
        self.drug_combin = pd.merge(self.drug_combin, drugs_net_distances, how='inner',
                 left_on='drug_pairs', right_on='drugs')

    def get_drug_inckey_from_start(self, ingre_tar_dict):
        self.drug_combin = pd.read_csv(self.drug_combin_file)
        drug_list = list(self.drug_combin['drug_row'].append(self.drug_combin['drug_col']).unique())[:-1]
        pd_drug = pd.DataFrame(drug_list, columns=['drugs'])
        pd_drug.to_csv('../herb_pairs_materials/case study/drugcom_dss/drug_com_list.csv')
        drug_inchey_pubchem_id_dict = {drug: drug_name_inchikey(drug) for drug in drug_list if drug is not np.nan}
        drug_inchey_dict = {k:'I' + v[0]for k,v in drug_inchey_pubchem_id_dict.items()}
        self.drug_pubchem_id = {k: v[1]for k,v in drug_inchey_pubchem_id_dict.items()}
        self.drug_inchey_dict = {k: v for k, v in drug_inchey_dict.items()
                            if v in ingre_tar_dict.keys() and v != 'Iunknown'}
        pd_drug_inchey_dict = pd.DataFrame.from_dict(self.drug_inchey_dict, orient='index').reset_index()
        pd_drug_inchey_dict.columns = ['drug_name', 'drug_inchey']
        pd_drug_inchey_dict.to_csv('../herb_pairs_materials/case study/drugcom_dss/pd_drug_inchey_dict.csv')
        pd_drug_pubcem_dict = pd.DataFrame.from_dict(self.drug_pubchem_id, orient='index').reset_index()
        pd_drug_pubcem_dict.columns = ['drug_name', 'drug_pubchem']
        pd_drug_pubcem_dict.to_csv('../herb_pairs_materials/case study/drugcom_dss/pd_drug_pubchem_dict.csv')

    def get_drug_targte_drucomm(self, drug_target_file):
        pass
        # need to prepare in local computer by drug_inckey dict. see codes in case study 3
    def get_drug_target_stitch(self):
        pass

    def read_filter(self, chrunk_pd):
        chrunk_pd_left = chrunk_pd[(chrunk_pd['drug_row'].isin(self.drug_inchey_dict.keys())) &
                         (chrunk_pd['drug_col'].isin(self.drug_inchey_dict.keys()))]
        chrunk_pd_left = self.add_inchkey_column(chrunk_pd_left)
        return chrunk_pd_left

    def get_simple_drug_combin_from_start(self):
        self.drug_combin = pd.concat([self.read_filter(chrunk_pd)
                                      for chrunk_pd in pd.read_csv(self.drug_combin_file, chunksize=10 ** 3)])
        self.drug_combin.to_csv('../herb_pairs_materials/case study/drugcom_dss/drug_com_simple.csv')


ingre_tar_dict = drug_ingredients_obj.ingre_tar_dict
drug_combin_file = '../herb_pairs_materials/case study/drugcom_dss/DataTable5e60faa831385690299.csv'

# from  start both of drug_inchey_dict and drugcom
drug_com_1 = DrugCom(drug_combin_file)
drug_com_1.get_drug_inckey_from_start(ingre_tar_dict) #need long tim eto map drug name to inchey
drug_com_1.get_simple_drug_combin_from_start()
drug_com_1.cal_drugs_net_distances_dict(drug_info_ogj.drug_ingredients_obj, g_obj)
drug_com_1.add_drugs_net_distances_all()


drug_target_info = '../herb_pairs_materials/case study/drugcom_dss/drug_target_yin_expand.csv'
drug_inco_v2 = Drug_Info()
drug_inco_v2.get_drug_info_yin(drug_target_info)

drug_inco_v2.get_gene_name_entryid_dict()
drug_inco_v2.add_gene_id_column()
drug_inco_v2.get_ingredient_ogj()
drug_inco_v2.drug_ingredients_obj.ingre_ingre_dis_all('IACTOXUHEUCPTEW-JMRHEKERSA-N', 'IAYJRTVVIBJSSKN-UHFFFAOYSA-N', g_obj.G)


drug_com_2 = DrugCom('../herb_pairs_materials/case study/drugcom_dss/DataTable5e60faa831385690299.csv')
drug_com_2.get_drug_name_inchey_dict('../herb_pairs_materials/case study/drugcom_dss/pd_drug_inchey_dict.csv')
drug_com_2.get_simple_drug_combin_from_start()
drug_com_2.get_drug_pairs()
drug_com_2.cal_drugs_net_distances_dict(drug_inco_v2.drug_ingredients_obj, g_obj.G.nodes)