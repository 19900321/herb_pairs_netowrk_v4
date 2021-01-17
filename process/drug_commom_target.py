import requests
import pandas as pd

from pandas.io.json import json_normalize


def prepare_filter(bioassay_dict):
    data = json_normalize(bioassay_dict)
    data = data[data['endpoint_standard_value'].astype(float) <= 1000]
    data.replace('', pd.np.nan, inplace=True)
    data.fillna(value=pd.np.nan, inplace=True)
    data = data.dropna(subset=['uniprot_id', 'gene_name'], how='all')
    data = data[data['activity_comment'] != 'inactive']
    return data



def joint_list(x):
    if x == None:
        return None
    else:
        return ','.join(list(set(x)))

from functools import reduce
def get_drugcomm_target(chemble_id):
    bioassay_dict=[]
    def get_targtes(chemble_id, url):
        url_2_start = 'https://drugtargetcommons.fimm.fi/api/data/bioactivity/?chembl_id={}&&format=json&&endpoint_standard_units=NM&&limit=1000&&offset=0'.format(
            chemble_id)
        if url == None:
            html_2 = requests.get(url_2_start).json()
        else:
            html_2 = requests.get(url).json()

        bioassay_dict.append(html_2['bioactivities'])

        if html_2['meta']['next'] != None:
            url_new = 'https://drugtargetcommons.fimm.fi' + html_2['meta']['next']
            get_targtes(chemble_id, url_new)

    get_targtes(chemble_id, None)
    bioassay_dict = reduce(lambda x,y: x+y, bioassay_dict )
    data = prepare_filter(bioassay_dict)
    data_left = data.groupby(['chembl_id', 'compound_name'])['uniprot_id'].apply(list).reset_index()
    data_left.columns = ['chembl_id', 'pert_iname', 'uniprot']
    data_left['uniprot'] = data_left['uniprot'].apply(joint_list)
    return data_left.loc[:, 'uniprot'][0]


#a2 = get_drugcomm_target('CHEMBL185')
#a3 = get_drugcomm_target('CHEMBL1976040')



