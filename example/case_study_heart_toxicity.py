from process.drug_annotation_pipline import *
import pandas as pd


# 1.get the drug information, including stitch target, drugcommon target
#in fact, this is done in local computer

def cid_info_dict(test1, info_list, cid):
    drug_info_dict = {'key_name': test1.drug}
    drug_info_dict.update(test1.properities_dict[cid])

    if 'drug_common_target' in info_list:
        drug_info_dict.update(
        {'drug_common_entrezgene': [i for i in test1.drug_common_targets_dict[cid]['entrezgene'] if i != None]})
    if 'drug_stitch_target' in info_list:
        drug_info_dict.update({'stitch_CIDm': test1.stitch_id_dict[cid]})
        drug_stitch_pd = pd.DataFrame(test1.stitch_targets_dict[cid])
        drug_stitch_pd.dropna(how='any')
        drug_info_dict.update({col: list(drug_stitch_pd[col]) for col in drug_stitch_pd.columns})
    return drug_info_dict

def get_drug_info(drug_name,info_list):
    drug_dict_list = []
    test1 = Mapping_Drugs(drug_name)
    test1.get_cid()
    test1.get_properties_cids()
    test1.get_unichem_ids_from_cids()
    if 'drug_common_target' in info_list:
        test1.get_drug_common_target()
    if 'drug_stitch_target' in info_list:
        test1.get_stitch_id_by_search()
        test1.get_stitch_target_cids()

    return [cid_info_dict(test1, info_list, cid) for cid in test1.cids]

def get_drug_info_pd(drug_name_list,info_list):

    for drug in drug_name_list:
        drug_dict_list = []
        if drug == np.nan:
            continue
        else:
            try:
                a = get_drug_info(drug, info_list)
            except:
                continue
            else:
                drug_dict_list += a

    drug_info_pd = pd.DataFrame(drug_dict_list)
    return drug_info_pd

# done
drug = pd.read_csv('../herb_pairs_materials/case study/heart toxicity/Drug_gene_list_YY.csv')

drug_info_pd = get_drug_info_pd(list(drug['CHEMBLID']))
drug_info_pd.to_csv('../herb_pairs_materials/case study/heart toxicity/drug_info_pd.csv', sep = '\t')

# prepareed to ingreidnets object
import ast
def string_list(x):
    key_list = [i for i in ast.literal_eval(x) if i !=None]
    return ','.join(key_list)

def drug_info_format(drug_info_pd):
    drug_info_pd['ingredients_id'] = drug_info_pd['cid']
    drug_info_pd['name'] = drug_info_pd['key_name']
    drug_info_pd['targets'] = drug_info_pd['entrezgene'].apply(string_list)
    drug_info_pd['ensymble'] = drug_info_pd['symbol'].apply(string_list)
    drug_info_pd['score'] = '10,10'

    drug_info_pd_ingre = drug_info_pd[['ingredients_id', 'name', 'targets', 'ensymble', 'score']]
    return drug_info_pd_ingre

#for drug positive
drug_info_pd = pd.read_csv('../herb_pairs_materials/case study/heart toxicity/drug_info_pd.csv', sep = '\t')
drug_info_pd_ingre = drug_info_format(drug_info_pd )
drug_info_pd_ingre.to_csv('../herb_pairs_materials/case study/heart toxicity/ingredients_format.csv', sep='\t')

drug_info_pd_negative = pd.read_csv('../herb_pairs_materials/case study/heart toxicity/drug_info_pd_2.csv', sep = '\t')
drug_info_pd_ingre_negative = drug_info_format(drug_info_pd_negative )
drug_info_pd_ingre_negative = drug_info_pd_ingre_negative.drop_duplicates(subset=['ingredients_id'])
drug_info_pd_ingre_negative.to_csv('../herb_pairs_materials/case study/heart toxicity/ingredients_negativeformat.csv', sep='\t')

from generate_objects import *
drug_ingredients_obj = Ingredients('../herb_pairs_materials/case study/heart toxicity/ingredients_format.csv', 0)
drug_ingredients_obj.ingredients_target_dict(g_obj.G.nodes)

drug_ingredients_negative_obj = Ingredients('../herb_pairs_materials/case study/heart toxicity/ingredients_negativeformat.csv', 0)

drug_ingredients_negative_obj.ingredients_target_dict(g_obj.G.nodes)

# transfer gene symbol to entrze name
import mygene
mg = mygene.MyGeneInfo()

def get_entries_from_symbol(symbol_id):
    result_list = [i['entrezgene'] for i in mg.querymany([symbol_id], scopes='symbol',
                                                  fields='entrezgene', species='human') if 'entrezgene' in i.keys()]
    if len(result_list) ==0:
        return None
    else:
        return result_list[0]

#drug = pd.read_csv('../herb_pairs_materials/case study/heart toxicity/Drug_gene_list_YY.csv')
# proxessed already, drug['targets'] = drug['Gene list'].apply(get_entries_from_symbol)
# drug.to_csv('../herb_pairs_materials/case study/heart toxicity/Drug_gene_list_YY_enrzgene.csv')

drug = pd.read_csv('../herb_pairs_materials/case study/heart toxicity/Drug_gene_list_YY_enrzgene.csv')

from example.casestudy2 import *
from disease import *
drug = drug.dropna(subset=['targets'])
new_disease_dict = {'heart_toxicity':['T' + str(int(i)) for i in drug['targets'] if i != None]}
g_obj.get_degree_binning(1001)
disease_file_name = '../herb_pairs_materials/case study/huangqitang/CTD_D008103_genes_20200306081710.csv'
disease = Disease(disease_file_name, g_obj)
disease.update_disease_target_dict(new_disease_dict)

from collections import defaultdict
dict_z = defaultdict()
for ingre in list(drug_ingredients_obj.ingre_tar_dict.keys()):
    d, z, m, s, pval = disease.cal_disease_ingre_z_score('heart_toxicity',
                                  ingre, 'closest', drug_ingredients_obj.ingre_tar_dict, 1000, 3333)

    dict_z[ingre] = {'d': d, 'z': z,'m': m, 'pval': list(pval)[0]}


from collections import defaultdict
dict_z = defaultdict()
for ingre in list(drug_ingredients_negative_obj.ingre_tar_dict.keys())[0:50]:
    d, z, m, s, pval = disease.cal_disease_ingre_z_score('heart_toxicity',
                                  ingre, 'closest', drug_ingredients_negative_obj.ingre_tar_dict, 1000, 3333)

    dict_z[ingre] = {'d': d, 'z': z,'m': m, 'pval': list(pval)[0]}
z_no_toxicity = pd.DataFrame.from_dict(dict_z, orient= 'index')
z_no_toxicity['cid'] = list(z_no_toxicity.index)
z_no_toxicity.to_csv('../herb_pairs_materials/case study/heart toxicity/no_toxicity.csv')

def get_name(x):
    return drug_ingredients_obj.ingredients_info(x)['name']

z_toxicity['chemble'] = z_toxicity['cid'].apply(get_name)
z_toxicity_pd = pd.merge(z_toxicity, drug, how='inner', left_on='chemble', right_on='CHEMBLID')
z_toxicity_pd.to_csv('../herb_pairs_materials/case study/heart toxicity/toxicity_2.csv')
