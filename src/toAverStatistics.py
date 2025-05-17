import sys
import pandas as pd
import os

input_list_file=sys.argv[1]

res_list=['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
residue_three_to_one_letter = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

#prepare pro-list list
input_list_list=[]
with open(f'{input_list_file}','r')as input:
    for line in input:
           input_list_list.append(line.strip())
#prepare all protein_toAS dir
if not os.path.exists('./All_protein_toAS'):
    os.mkdir('./All_protein_toAS')
#calculate all protein_toAS data
for model in ('6A_ball','4.5A_ball','4.5A_to_6A_shell'):
    result_sum = {}
    result_sum_count = {}
    result_sum_pro = {}
    result_sum_2 = {}
    result_sum_pro_2 = {}
    for list_file in input_list_list:
        if os.path.exists(list_file):
            with open(f'{list_file}','r') as file:
                for line in file:
                    pdb_file=line.strip()
                    if os.path.exists(f'./FEATURE.feat/{list_file}/{model}/{pdb_file}_{model}.csv'):
                        if model == '6A_ball' or model == '4.5A_ball':
                            feat = pd.read_csv(
                                f'./FEATURE.feat/{list_file}/{model}/{pdb_file}_{model}.csv',
                                header=None)
                            Feature_res_name = feat.iloc[:, -1].tolist()
                            res_name = []
                            res_id = []
                            Res = []
                            id = 1
                            for i in Feature_res_name:
                                i.split(':')
                                Res.append(i[:3])
                                res_name.append(residue_three_to_one_letter[i[:3]])
                                res_id.append(id)
                                id += 1
                            for ID in range(1, len(res_id) + 1):
                                res_element = Res[ID - 1]
                                if res_element in result_sum_count:
                                    result_sum_count[res_element] += 1
                                else:
                                    result_sum_count[res_element] = 1
                                if res_element in result_sum.keys():
                                    result_sum[res_element] += feat.iloc[ID - 1, 1:55].astype('float64')
                                else:
                                    result_sum[res_element] = feat.iloc[ID - 1, 1:55].astype('float64')
                                if res_element in result_sum_2.keys():
                                    result_sum_2[res_element] += feat.iloc[ID - 1, 56:81].astype('float64')
                                else:
                                    result_sum_2[res_element] = feat.iloc[ID - 1, 56:81].astype('float64')

                        if model == '4.5A_to_6A_shell':
                            feat = pd.read_csv(
                                f'./FEATURE.feat/{list_file}/{model}/{pdb_file}_{model}.csv',
                                header=None)
                            Feature_res_name = feat.iloc[:, -1].tolist()
                            res_name = []
                            res_id = []
                            Res = []
                            id = 1
                            for i in Feature_res_name:
                                i.split(':')
                                Res.append(i[:3])
                                res_name.append(residue_three_to_one_letter[i[:3]])
                                res_id.append(id)
                                id += 1
                            for ID in range(1, len(res_id) + 1):
                                res_element = Res[ID - 1]
                                if res_element in result_sum_count:
                                    result_sum_count[res_element] += 1
                                else:
                                    result_sum_count[res_element] = 1
                                if res_element in result_sum.keys():
                                    result_sum[res_element] += feat.iloc[ID - 1, 241:295].astype('float64')
                                else:
                                    result_sum[res_element] = feat.iloc[ID - 1, 241:295].astype('float64')
                                if res_element in result_sum_2.keys():
                                    result_sum_2[res_element] += feat.iloc[ID - 1, 296:321].astype('float64')
                                else:
                                    result_sum_2[res_element] = feat.iloc[ID - 1, 296:321].astype('float64')
        else:
            with open('Not_exist_list.log','a')as wrong_list:
                wrong_list.write(f'{list_file}+\n')
    for residue in res_list:
        if residue in result_sum.keys():
            result_sum_pro[residue] = result_sum[residue] / result_sum_count[residue]
            result_sum_pro_2[residue] = result_sum_2[residue] / result_sum_count[residue]
        else:
            result_sum_pro[residue] = pd.Series([0]*54, dtype = 'float64')
            result_sum_pro_2[residue] = pd.Series([0]*25, dtype='float64')

    content1 = ''
    for residue in res_list:
        feat_1 = ','.join([f"{num:.8f}" for num in result_sum_pro[residue]])
        content1 += feat_1 + '\n'
    with open(f'./All_protein_toAS/data_{model}_toAS1.csv', "w") as f:
        f.write(content1)
    content2 = ''
    for residue in res_list:
        feat_2 = ','.join([f"{num:.8f}" for num in result_sum_pro_2[residue]])
        content2 += feat_2 + '\n'
    with open(f'./All_protein_toAS/data_{model}_toAS2.csv', "w") as f:
        f.write(content2)


    df1 = pd.read_csv(f'./All_protein_toAS/data_{model}_toAS1.csv',header=None)
    df1.insert(0, 'new_column', res_list)
    df1.to_csv(f'./All_protein_toAS/data_{model}_toAS1.csv', index=False,header=False)
    df2 = pd.read_csv(f'./All_protein_toAS/data_{model}_toAS2.csv',header=None)
    df2.insert(0, 'new_column', res_list)
    df2.to_csv(f'./All_protein_toAS/data_{model}_toAS2.csv', index=False,header=False)
    print(f'\n\ntoAverStatistics.py--{model}:The mutagenesis characteristics of the protein database have been calculated!')
