#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 13.3.2019 16.51
# @Author  : YINYIN
# @Site    : 
# @File    : drug_annotation_pipline.py
# @Software: PyCharm

import pubchempy as pcp

#from chembl_webresource_client.new_client import new_client
#molecule = new_client.molecule
#drug_indication = new_client.drug_indication
import urllib
import json
#import html5lib
#import lxml.html as lh

import mygene
mg = mygene.MyGeneInfo()

from process.drug_commom_target import *

def delete_suffix(x):
    return x.split('.')[1].strip()

def get_entries_from_ensembl(ensembl_protein):
    result_list = [i['entrezgene'] for i in mg.querymany([ensembl_protein], scopes='ensembl.protein',
                                                  fields='entrezgene', species='human') if 'entrezgene' in i.keys()]
    if len(result_list) ==0:
        return None
    else:
        return result_list[0]

def get_symble_from_entries(entry):
    result_list = [i['symbol'] for i in mg.querymany([entry], scopes='entrezgene', fields='symbol') if 'symbol' in i.keys()]
    if len(result_list) == 0:
        return None
    else:
        return result_list[0]


def get_entries_from_uniprot(uniprot_id):
    result_list = [i['entrezgene'] for i in mg.querymany([uniprot_id], scopes='uniprot',
                                                  fields='entrezgene', species='human') if 'entrezgene' in i.keys()][0]
    if len(result_list) ==0:
        return None
    else:
        return result_list[0]

def unichem_trans(original_value, original_id_number, to_id_number):
    url_2 = 'https://www.ebi.ac.uk/unichem/rest/src_compound_id/{}/{}/{}'.format(original_value, original_id_number,
                                                                                 to_id_number)
    try:
        html_2 = requests.get(url_2).content
        responseStr_2 = html_2.decode("utf-8")
    except:
        print('error:no html')
        return None
    else:
        json_acceptable_string_2 = responseStr_2.replace("'", "\"")
        try:
            d_2 = json.loads(json_acceptable_string_2)
        except:
            print('error2:no html read')
            return None
        else:
            if len(d_2) == 0:
                print('no information')
                return None
            else:
                id_get_2 = list(d_2[0].values())[0]
                return id_get_2

def from_inchikey_to_stitch_id(inchikeys):
    import MySQLdb
    # connect mysql
    db_2 = MySQLdb.connect(host="127.0.0.1", user="yin", passwd="Mqxs320321wyy", db="stitch")
    c = db_2.cursor()
    format_strings = ','.join(['%s'] * len(inchikeys))
    c.execute("SELECT * from inchey where inchekey in (%s)" % format_strings, tuple(inchikeys))

    inchey_used_2 = c.fetchall()
    pd_result = pd.DataFrame(list(inchey_used_2), columns=['stitch_CIDm', 'stitch_CIDs',
                                                           'stitch_cid', 'inchikey'])
    inchikey_stitch_id_dict = dict(zip(pd_result['inchikey'], pd_result['stitch_CIDm']))
    return inchikey_stitch_id_dict


def get_stitch_target(stitch_id):

    url = 'http://stitch.embl.de/api/psi-mi-tab/interactions?identifier={}&species=9606'.format(stitch_id)
    try:
        responseStr = urllib.request.urlopen(url).read()
        responseStr = responseStr.decode("utf-8")
    except:
        print('no cid infor')
        return None
    else:
        pd_network = pd.DataFrame([j.split('\t') for j in responseStr.split('\n')]).iloc[:,0:4]
        pd_network.columns = ['nodes1', 'nodes2','name1','name2']
        pd_network = pd_network.dropna()

        pd_network['nodes1'] = pd_network['nodes1'].apply(delete_suffix)
        pd_network['nodes2'] = pd_network['nodes2'].apply(delete_suffix)

        tuple_dict = set(list(zip(pd_network['nodes1'].append(pd_network['nodes2']),
                                  pd_network['name1'].append(pd_network['name2']))))

        return [{'ensembl.protein': i[0], 'symbol': i[1], 'entrezgene': get_entries_from_ensembl(i[0])} for i in tuple_dict if not i[0].startswith('CID')]

#1. Get cid by CAS or name , then get  Pubchem_inchikey,Pubchem_molecular_formula,Pubchem_smiles,
#Pubchem_synonyms,
class Mapping_Drugs:
    def __init__(self, drug):
        self.drug = drug

    def get_cid(self):
        cid_records = pcp.get_synonyms(self.drug, 'name')
        self.cids = [cid['CID'] for cid in cid_records]


    def get_properties(self, drug_cid):

        c = pcp.Compound.from_cid(drug_cid)
        property_dict = c.to_dict(['cid', 'canonical_smiles',
                                   'molecular_formula', 'inchikey', 'iupac_name'])
        return property_dict

    def get_properties_cids(self):
        if len(self.cids) != 0:
            self.properities_dict = {cid: self.get_properties(cid) for cid in self.cids}
        else:
            self.properities_dict = None

    def get_unichem_ids_from_cid(self, cid):
        unichem_dict = {'chembl_id':unichem_trans(cid, 22, 1),
        'drugbank_id': unichem_trans(cid, 22, 2),
        'kegg_id':unichem_trans(cid, 22, 6),
        'chebi_i':unichem_trans(cid, 22, 7)}
        return unichem_dict

    # 2.Get mutiple other id from one id unichem,for example Unichem_chembl_id,
    # DrugBank ID,Unichem_kegg_compoundid,Unichem_chebi_id

    def get_unichem_ids_from_cids(self):
        if len(self.cids) != 0:
            self.unichem_dict = {cid: self.get_unichem_ids_from_cid(cid) for cid in self.cids}
        else:
            self.unichem_dict = None

    def get_stitch_id_by_dict(self,inchikey_stitch_id_dict):
        if len(self.cids) != 0:
            self.stitch_id_dict = {cid: inchikey_stitch_id_dict[self.unichem_dict[cid]['inchikey']] for cid in self.cids}
        else:
            self.stitch_id_dict = None

    def get_stitch_id_by_search(self):
        self.stitch_id_dict = {cid:list(from_inchikey_to_stitch_id([self.properities_dict[cid]['inchikey']]).values())[0] for cid in self.cids}


    #4. Get stitch compound target and targets interactions and score
    def get_stitch_target_cids(self):

        self.stitch_targets_dict = {cid: get_stitch_target(self.stitch_id_dict[cid]) for cid in self.cids}

    def get_drug_common_target(self):
        drug_common_targets_dict = {cid: get_drugcomm_target(self.unichem_dict[cid]['chembl_id'])
                                                            for cid in self.cids}

        self.drug_common_targets_dict ={cid:{'uniprot':proteins.split(','),'entrezgene':[get_entries_from_uniprot(protein) for protein in proteins.split(',')]} for cid, proteins in drug_common_targets_dict.items()}

    def get_drug_info(self,info_list):
        self.get_cid()
        self.get_properties_cids()
        self.get_unichem_ids_from_cids()
        if 'drug_common_target' in info_list:
            self.get_drug_common_target()
        if 'drug_stitch_target' in info_list:
            self.get_stitch_id_by_search()
            self.get_stitch_target_cids()


# #4.get drug_target_interactions from pubchem
#
# def get_drug_targets_interaction(compound_df,cid_col_name):
#     for id in compound_df.index:
#         print(id)
#         compound_item = compound_df.loc[id]
#         cid = int(compound_df.loc[id, cid_col_name])
#
#         url_2 = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/' + str(
#             cid) + '/JSON?heading=DrugBank+Interactions'
#         try:
#             responseStr = urllib.request.urlopen(url_2).read()
#             responseStr = responseStr.decode("utf-8")
#         except:
#             print('no drug target')
#         else:
#             d_2 = json.loads(responseStr)
#             d_3 = d_2['Record']["Section"][0]["Section"][0]["Information"]
#             pd_try = pd.DataFrame(d_3)
#             pd_try_item = [[i['Row'][0]["Cell"][0]["StringValue"]] + [
#                 i['Row'][0]["Cell"][1]["StringValue"].replace('</a>', '').split('">')[-1]] + [
#                                i['Row'][1]["Cell"][1]["StringValue"]] for i in pd_try.iloc[:, 2]]
#             pd_try_item_name = [[i['Row'][0]["Cell"][1]["StringValue"].replace('</a>', '').split('">')[-1]] for i in
#                                 pd_try.iloc[:, 2]]
#             compound_df.loc[id, 'Pubchem_drug_target_all'] = str(pd_try_item_name)
#             compound_df.loc[id, 'Pubchem_drug_target_all_real'] = str(pd_try_item)
#             list_unipro = []
#             list_symbol = []
#             for i in pd_try.iloc[:, 2]:
#                 row_infor = i['Row']
#                 for k in row_infor:
#                     if k["Cell"][0]["StringValue"] == "PubChem Protein Target":
#                         unipro_target = k["Cell"][1]["StringValue"]
#                         list_unipro.append(unipro_target)
#                     if k["Cell"][0]["StringValue"] == "PubChem Gene Target":
#                         unipro_gene = k["Cell"][1]["StringValue"]
#                         list_symbol.append(unipro_gene)
#             compound_df.loc[id, 'Pubchem_drug_target_uniprot_id'] = str(list_unipro)
#             compound_df.loc[id, 'Pubchem_drug_target_symbol'] = str(list_symbol)
#     return compound_df
#
# #5. get drugbank annotation for compounds, fill more blank drugbankid,keggid
#
# ##5.1 get other ids by CAS map from links of drugbank
# def get_link_from_drugbank(compound_df,CAS_col_name):
#     data_drug_map = pd.read_csv(
#         'C:\\Users\\yinyin\\Desktop\\Project\\drugcom_1135\\target\\drugbank_small_molecule_drug_links.csv\\drug links.csv')
#     data_drug_map.columns = ['drugbank_link_' + i for i in data_drug_map.columns]
#     compound_df_3 = pd.concat([compound_df, pd.DataFrame(columns=data_drug_map.columns)], sort=False)
#
#     for id in compound_df_3.index:
#         print(id)
#         compound_item = compound_df_3.loc[id]
#         if pd.isna(compound_df_3.loc[id,CAS_col_name]) == False:
#             cas_no = compound_df_3.loc[id,CAS_col_name].strip()
#             if cas_no in data_drug_map.iloc[:, 2].tolist():
#                 drug_id_map_item = data_drug_map[data_drug_map.iloc[:, 2] == cas_no]
#                 compound_df_3.loc[id, data_drug_map.columns[0]:data_drug_map.columns[-1]] = drug_id_map_item.iloc[0, :].tolist()
#     return compound_
#
# #5.2 combined drugbank from unichem and drugbank link
# def drugbank_id_combined(compound_df,unichem_drugbankid_col_name,drugbank_link_drugbankid):
#     compound_df['combined_drugbank_id_combined'] = compound_df[drugbank_link_drugbankid].fillna(compound_df[unichem_drugbankid_col_name])
#     return compound_df
#
# ##5.3 get indication and description from drugbank
#
# def get_drugbank_indcation(compound_df_4,drugbank_id_col_name):
#     import xml.etree.ElementTree as ET
#     e = ET.parse('drugbank_all_full_database.xml\\full database.xml').getroot()
#     DBs = compound_df_4[drugbank_id_col_name].tolist()
#     for i in e:
#         db_id = i[0].text
#         if db_id in DBs:
#             index_use_2 = compound_df_4.index[compound_df_4[drugbank_id_col_name] == db_id].tolist()
#             indication_drug = i.find('{http://www.drugbank.ca}indication').text
#             compound_df_4.loc[index_use_2, 'drugbank_indication'] = indication_drug
#     return compound_df_4
#
# def get_drugbank_description(compound_df_4,drugbank_id_col_name):
#     import xml.etree.ElementTree as ET
#     e = ET.parse('drugbank_all_full_database.xml\\full database.xml').getroot()
#     DBs = compound_df_4[drugbank_id_col_name].tolist()
#     for i in e:
#         db_id = i[0].text
#         if db_id in DBs:
#             print('in_', db_id)
#             index_use_2 = compound_df_4.index[compound_df_4[drugbank_id_col_name] == db_id].tolist()
#             description_drug = i.find('{http://www.drugbank.ca}description').text
#             compound_df_4.loc[index_use_2, 'drugbank_description'] = description_drug
#     return compound_df_4
#
#
#
# #6. Get fill, get hembl_mechanism,Chembl_mechanism_type,chemble indication
# ##6.1 get from pubchem
# def chemble_id_by_pubchem(compound_df,inchey_col_name,cid_col_name,synomy_col_name):
#     for id in compound_df.index:
#         print(id)
#         compound_item = compound_df.loc[id]
#         cid = int(compound_df.loc[id, cid_col_name])
#         inchikey = compound_df.loc[id, inchey_col_name]
#         synomy = compound_df.loc[id, synomy_col_name]
#     try:
#         chemble_infor = molecule.get(inchikey)
#     except:
#         try:
#             chemble_infor = molecule.get(smiles)
#         except:
#             chembl_id = 'NA'
#         else:
#             chembl_id = chemble_infor['molecule_chembl_id']
#     else:
#         chembl_id = chemble_infor['molecule_chembl_id']
#
#     compound_df.loc[id, 'chemble_id_from_pubchem'] = str(chembl_id)
#
#     return compound_df
#
# ##6.2 combine chemble id from two way, one from pubchem one from unichem
# def chemble_id_combined(compound_df,inchey_col_name,cid_col_name,synomy_col_name,id_number, col_name_used, original_id_number,id_name_add):
#     compound_df_1 = unichem_trans(id_number, compound_df, col_name_used, original_id_number,id_name_add)
#     compound_df_2 = chemble_id_by_pubchem(compound_df_1,inchey_col_name,cid_col_name,synomy_col_name)
#     compound_df_2 = compound_df_2[id_name_add].fillna(compound_df_2['chemble_id_from_pubchem'])
#     return compound_df_2
#
# ##6.3 get indication from chemble
#
# def get_mechanism_from_chemble_yin(compound_df,id_name_add):
#     for id in compound_df.index:
#         print(id)
#         compound_item = compound_df.loc[id]
#         cheml_id = compound_item.loc[id,id_name_add]
#         if cheml_id != 'NA':
#             cheml_link = 'https://www.ebi.ac.uk/chembl/api/data/mechanism.json?molecule_chembl_id__exact=' + str(
#                 cheml_id)
#             responseStr = urllib.request.urlopen(cheml_link).read()
#             responseStr = responseStr.decode("utf-8")
#
#             pattern = re.compile(r'mechanism_of_action":.*?,')
#             result = pattern.findall(responseStr)
#             cheml_mechaisms = [i[23:][:-2].strip() for i in result]
#
#             pattern_2 = re.compile(r'action_type": .*?,')
#             result_2 = pattern_2.findall(responseStr)
#             cheml_mechaisms_type = [i[15:][:-2].strip() for i in result_2]
#             compound_df.loc[id, 'chembl_mechanism'] = str(cheml_mechaisms)
#             compound_df.loc[id, 'chembl_mechanism_type'] = str(cheml_mechaisms_type)
#     return compound_df
#
# ##6.4 get indication from chemble
# from chembl_webresource_client.new_client import new_client
# molecule = new_client.molecule
# drug_indication = new_client.drug_indication
#
# def get_indication_from_chemble(compound_df,id_name_add):
#     indic_all = pd.read_csv('C:\\Users\\yinyin\\Desktop\\herbpair\\14drugbank\\chemble\\chemble_indication.txt',
#                             sep='\t')
#     indic_all2 = pd.DataFrame(indic_all,id_name_add)
#     for id in compound_df.index:
#         print(id)
#         compound_item = compound_df.loc[id]
#         cheml_id = compound_item.loc[id,id_name_add]
#         if chembl_id in indic_all2['MOLECULE_CHEMBL_ID'].tolist():
#             mesh_id_list = indic_all2.loc[indic_all2['MOLECULE_CHEMBL_ID'] == chembl_id]
#             mesh_id = mesh_id_list['MESH_ID'].tolist()
#             mesh_head_id = mesh_id_list['MESH_HEADING'].tolist()
#             compound_df.loc[id, 'MESH_ID'] = str(mesh_id)
#             compound_df.loc[id, 'MESH_HEADING'] = str(mesh_head_id)
#     return compound_df
#
# ## 6.5 get predicted targets from chemble
# # url https://www.ebi.ac.uk/chembl/api/data/target_prediction.json?molecule_chembl_id=CHEMBL939&value=10&limit=200
# def get_predict_target_from_chemble(compound_df,level_of_concentration,numer_of_target,id_name_add):
#     for id in compound_df.index:
#         print(id)
#         compound_item = compound_df.loc[id]
#         cheml_id = compound_item.loc[id,id_name_add].strip()
#         url_2 = 'https://www.ebi.ac.uk/chembl/api/data/target_prediction.json?molecule_chembl_id='+str(cheml_id)+'&value='+str(level_of_concentration)\
#                 +'&limit='+str(numer_of_target)
#     try:
#         responseStr = urllib.request.urlopen(url_2).read()
#         responseStr = responseStr.decode("utf-8")
#     except:
#         print('no web')
#     else:
#         d_2 = json.loads(responseStr)["target_predictions"]
#         targets = [i['target_accession'] for i in d_2 ]
#         geneSyms = mg.querymany(targets, scopes='uniprot', fields='symbol', species='human')
#         geneSyms_sym = [i['symbol'] for i in geneSyms if len(i) != 2]
#
#         compound_df.loc[id, 'chemble_pretarget_uni'] = str(targets)
#         compound_df.loc[id, 'chemble_pretarget_symbol'] = str(geneSyms_sym)
#     return compound_df
#
# #7.get Pharmacological Classes (Pubchem_MoA(Mechanisms of Action), Pubchem_EPC(Established Pharmacologic Class) ,
# # Pubchem_mesh_name, Pubchem_mesh_id, MeSH Pharmacological Classification
#
# ##7.1 get Pubchem_MoA(Mechanisms of Action), Pubchem_EPC(Established Pharmacologic Class)
#
# def get_pharm_mechism_from_pubchem(compound_df,cid_col_name):
#     for id in compound_df.index:
#         print(id)
#         compound_item = compound_df.loc[id]
#         puchem_id = compound_item.loc[id,cid_col_name]
#         url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/' + str(
#             puchem_id) + '/JSON?heading=Pharmacological+Classes'
#         try:
#             responseStr = urllib.request.urlopen(url).read()
#             responseStr = responseStr.decode("utf-8")
#         except:
#             print('2_situ1:no web')
#         else:
#             responseStr = responseStr.replace('\n', '').replace('  ', '').replace('[', '').replace(']', '').replace(
#                 '\"', '').replace('}', '')
#             pattern_1 = re.compile(r'MoA,StringValue:.*?,')
#             result_1 = pattern_1.findall(responseStr)
#             pattern_2 = re.compile(r'EPC,StringValue:.*?,')
#             result_2 = pattern_2.findall(responseStr)
#             if len(result_1) == 0:
#                 print('result_1 map_nothing')
#             else:
#                 list_create = []
#                 for i in result_1:
#                     k = i.split('</a>')
#                     if len(k) >= 2:
#                         list_create.append(k[-1])
#                     else:
#                         try:
#                             list_create.append(i[17:-1])
#                         except:
#                             print('no enought length')
#                         else:
#                             list_create.append(i[17:-1])
#                     compound_df.loc[id, 'MoA'] = str(list_create)
#             if len(result_2) == 0:
#                 print('result_2 map_nothing')
#             else:
#                 list_create = []
#                 for i in result_2:
#                     k = i.split('</a>')
#                     if len(k) >= 2:
#                         list_create.append(k[-1])
#                     else:
#                         try:
#                             list_create.append(i[17:-1])
#                         except:
#                             print('no enought length')
#                         else:
#                             list_create.append(i[17:-1])
#                     compound_df.loc[id, 'EPC'] = str(list_create)
#     return compound_df
#
# ##7.2 get Pubchem_mesh_id, MeSH Pharmacological Classification
# def get_indication_from_pubchem(compound_df,cid_col_name):
#     for id in compound_df.index:
#         print(id)
#         compound_item = compound_df.loc[id]
#         puchem_id = compound_df.loc[id,cid_col_name]
#         url_2 = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/' + str(
#             puchem_id) + '/JSON?heading=MeSH+Pharmacological+Classification'
#         try:
#             html_2 = requests.get(url_2).content
#             responseStr_2 = html_2.decode("utf-8")
#             json_acceptable_string_2 = responseStr_2.replace("'", "\"")
#             d_2 = json.loads(json_acceptable_string_2)
#         except:
#             print('2_situ1:no mesh')
#         else:
#             json_acceptable_string_2 = responseStr_2.replace("'", "\"")
#             d_2 = json.loads(json_acceptable_string_2)
#             if len(d_2) != 0:
#                 try:
#                     chembi_id_get_2 = d_2['Record']['Reference']
#                     mesh_name_list = [i['Name'] for i in chembi_id_get_2]
#                     mesh_id_list = [i['SourceID'] for i in chembi_id_get_2]
#                 except:
#                     print('no record')
#                 else:
#                     compound_df.loc[id, 'pubchem_mesh_name'] = str(mesh_name_list)
#                     compound_df.loc[id, 'pubcjem_mesh_id'] = str(mesh_id_list)
#     return compound_df
#
# #8 get other inforation from broad, need to add cid then use it directly next time
# def get_infor_from_broad(compound_df,cid_col_name):
#     broad_pd = pd.read_csv('C:\\Users\\yinyin\\Desktop\\Project\\drugcom_1135\\drug_repurseing\\broad_modify.txt', sep='\t')
#     columns_changed = ['broad_'+i for i in broad_pd.columns]
#     compound_df_4 = pd.concat([compound_df, pd.DataFrame(columns=columns_changed)], sort=False)
#
#     for id in compound_df_4.index:
#         print(id)
#         CID =  compound_df_4.loc[id, cid_col_name]
#         if CID in broad_pd['cid'].tolist() and CID != 'NA':
#             add_part = broad_pd[broad_pd['cid']==CID].tolist()
#             compound_df_4.loc[id, columns_changed[0]:columns_changed[-1]] = add_part
#     return compound_df_4
#
#
# import json
# import MySQLdb
#
# connect = MySQLdb.connect(host="127.0.0.1", user="yin", passwd="Mqxs320321wyy", db="TCM_infor")
# c = connect.cursor()
# c.execute('SELECT * from herb_compound')
# herb_compounds = c.fetchall()
# inchey_used_list = [list(i) for i in herb_compounds]
# herb_compounds = pd.DataFrame(inchey_used_list[1:],columns=inchey_used_list[0])
# herb_compounds['cid\r'] = herb_compounds['cid\r'].str.strip()
# herb_compounds.columns.tolist()[5] = herb_compounds.columns.tolist()[5].strip()
#
#
# ##turn inchikey to stitch cid from shuyu
# data_new_drugs = pd.read_csv('C:\\Users\\yinyin\\Desktop\\Project\\drugcom_1135\\drugstitchid\\drug.csv')
#
#
# def form_cid_to_stitch_id(example_2, inchikey_col_name):
#     import json
#     import MySQLdb
#
#     list_drug = example_2[inchikey_col_name].tolist()
#     list_drug_2 = [x for x in list_drug if str(x) != 'nan']
#     # connect mysql
#     db_2 = MySQLdb.connect(host="127.0.0.1", user="yin", passwd="Mqxs320321wyy", db="stitch")
#     c = db_2.cursor()
#     format_strings = ','.join(['%s'] * len(list_drug_2))
#     c.execute("SELECT * from inchey where inchekey in (%s)" % format_strings, tuple(list_drug_2))
#     inchey_used_2 = c.fetchall()
#     inchey_used_list = [list(i) for i in inchey_used_2]
#     pd_used = pd.DataFrame(inchey_used_list)
#     index_with_example_2 = example_2.index[example_2[inchikey_col_name].notna()]
#     for j in index_with_example_2:
#         inchikey_value = example_2.loc[j, inchikey_col_name]
#         if inchikey_value in pd_used[3].tolist():
#             index_used = pd_used.index[pd_used[3] == inchikey_value].tolist()
#             stitch_cid_m = pd_used.loc[index_used[0]][0]
#             stitch_cid_s = pd_used.loc[index_used[0]][1]
#             example_2.loc[j, 'Stitch_cid_m'] = stitch_cid_m
#             example_2.loc[j, 'Stitch_cid_s'] = stitch_cid_s
#     return example_2
#
# maped_stitch_id = form_cid_to_stitch_id(data_new_drugs,'inchikey')
#
#
# def get_stich_name_yin(compound_df,stitch_id_col_name):
#     use_index = compound_df.index[compound_df[stitch_id_col_name].notna()]
#     for id in use_index:
#         print(id)
#         stich = compound_df.loc[id, stitch_id_col_name]
#         url = 'http://stitch.embl.de/api/json/resolve?identifier='+str(stich) +'&species=9606'
#         try:
#             responseStr = urllib.request.urlopen(url).read().decode("utf-8")
#         except:
#             print('no cid infor',id)
#             compound_df.loc[id,'stitch_name'] = np.nan
#         else:
#             try:
#                 html_2 = ast.literal_eval(responseStr)
#                 stitch_name = html_2[0]['preferredName']
#             except:
#                 compound_df.loc[id, 'stitch_name'] = np.nan
#             else:
#                 compound_df.loc[id, 'stitch_name'] = stitch_name
#     return compound_df