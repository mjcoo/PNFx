import pandas as pd
import os
import sys
import numpy as np

feature_dir = os.environ.get('FEATURE_DIR')
input_file_name = sys.argv[1]
model = sys.argv[2]

res_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO',
            'SER', 'THR', 'TRP', 'TYR', 'VAL']

residue_three_to_one_letter = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

feature_list = pd.read_csv(f'{feature_dir}/Feature_feature_{model}_toAS_list.csv', header=None)
pro = ','.join(feature_list.iloc[:, 0].astype(str))


def read_csv_and_assign(csv_file):
    data = pd.read_csv(csv_file, header=None)
    result = {}
    for index, row in data.iterrows():
        key = row.iloc[0]
        value = row.iloc[1:]
        result[key] = value
    return result


with open(f'{input_file_name}', 'r') as file:
    for line in file:
        pdb_file = line.strip()
        if os.path.exists(f'./FEATURE.feat/{input_file_name}/{model}/{pdb_file}_{model}.csv'):
            result = {}
            result_count = {}
            result_pro = {}
            result_2 = {}
            result_pro_2 = {}
            if model == '6A_ball' or model == '4.5A_ball':
                feat = pd.read_csv(f'./FEATURE.feat/{input_file_name}/{model}/{pdb_file}_{model}.csv', header=None)
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
                with open(f'./toAS.feat/{input_file_name}/{model}/{pdb_file}_{model}_toAS.feat', "a") as f:
                    f.write(f'{pro}')
                    f.write('\n')

                    for ID in range(1, len(res_id) + 1):
                        res_element = Res[ID - 1]
                        if res_element in result_count:
                            result_count[res_element] += 1
                        else:
                            result_count[res_element] = 1
                        if res_element in result.keys():
                            result[res_element] += feat.iloc[ID - 1, 1:55].astype('float64')
                        else:
                            result[res_element] = feat.iloc[ID - 1, 1:55].astype('float64')
                        if res_element in result_2.keys():
                            result_2[res_element] += feat.iloc[ID - 1, 56:81].astype('float64')
                        else:
                            result_2[res_element] = feat.iloc[ID - 1, 56:81].astype('float64')
                    for item in res_list:
                        if item in result.keys():
                            result_pro[item] = result[item] / result_count[item]
                            result_pro_2[item] = result_2[item] / result_count[item]

                    for ID in range(1, len(res_id) + 1):
                        for Remaining_residues in res_list:
                            if Res[ID - 1] != Remaining_residues:
                                if Remaining_residues in result.keys():

                                    numpy_arr1 = np.array(result_pro[Remaining_residues])
                                    numpy_arr2 = np.array(feat.iloc[ID - 1, 1:55]).astype('float64')
                                    feat_tmp1 = np.subtract(numpy_arr1, numpy_arr2).astype('float64')
                                    feat_1 = ','.join([f"{num:.3f}" for num in feat_tmp1])

                                    numpy_arr3 = np.array(result_pro_2[Remaining_residues])
                                    numpy_arr4 = np.array(feat.iloc[ID - 1, 56:81]).astype('float64')
                                    feat_tmp2 = np.subtract(numpy_arr3, numpy_arr4).astype('float64')
                                    feat_2 = ','.join([f"{num:.3f}" for num in feat_tmp2])

                                    f.write(
                                        f'{feat_1},{feat_2}')  # {res_name[ID-1]}{ID}{residue_three_to_one_letter[Remaining_residues]}
                                    f.write("\n")
                                else:
                                    data_feat1 = read_csv_and_assign(f'./All_protein_toAS/data_{model}_toAS1.csv')
                                    data_feat2 = read_csv_and_assign(f'./All_protein_toAS/data_{model}_toAS2.csv')

                                    numpy_arr1 = np.array(data_feat1[Remaining_residues])
                                    numpy_arr2 = np.array(feat.iloc[ID - 1, 1:55]).astype('float64')
                                    feat_tmp1 = np.subtract(numpy_arr1, numpy_arr2).astype('float64')
                                    feat_1 = ','.join([f"{num:.3f}" for num in feat_tmp1])

                                    numpy_arr3 = np.array(data_feat2[Remaining_residues])
                                    numpy_arr4 = np.array(feat.iloc[ID - 1, 56:81]).astype('float64')
                                    feat_tmp2 = np.subtract(numpy_arr3, numpy_arr4).astype('float64')
                                    feat_2 = ','.join([f"{num:.3f}" for num in feat_tmp2])

                                    f.write(
                                        f'{feat_1},{feat_2}')  # {res_name[ID-1]}{ID}{residue_three_to_one_letter[Remaining_residues]}
                                    f.write("\n")

            if model == '4.5A_to_6A_shell':
                feat = pd.read_csv(
                    f'./FEATURE.feat/{input_file_name}/{model}/{pdb_file}_{model}.csv',
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
                with open(
                        f'./toAS.feat/{input_file_name}/{model}/{pdb_file}_{model}_toAS.feat',
                        "a") as f:
                    f.write(f'{pro}')
                    f.write('\n')

                    for ID in range(1, len(res_id) + 1):
                        res_element = Res[ID - 1]
                        if res_element in result_count:
                            result_count[res_element] += 1
                        else:
                            result_count[res_element] = 1
                        if res_element in result.keys():
                            result[res_element] += feat.iloc[ID - 1, 241:295].astype('float64')
                        else:
                            result[res_element] = feat.iloc[ID - 1, 241:295].astype('float64')
                        if res_element in result_2.keys():
                            result_2[res_element] += feat.iloc[ID - 1, 296:321].astype('float64')
                        else:
                            result_2[res_element] = feat.iloc[ID - 1, 296:321].astype('float64')
                    for item in res_list:
                        if item in result.keys():
                            result_pro[item] = result[item] / result_count[item]
                            result_pro_2[item] = result_2[item] / result_count[item]

                    for ID in range(1, len(res_id) + 1):
                        for Remaining_residues in res_list:
                            if Res[ID - 1] != Remaining_residues:
                                if Remaining_residues in result.keys():
                                    numpy_arr1 = np.array(result_pro[Remaining_residues])
                                    numpy_arr2 = np.array(feat.iloc[ID - 1, 241:295]).astype('float64')
                                    feat_tmp1 = np.subtract(numpy_arr1, numpy_arr2).astype('float64')
                                    feat_1 = ','.join([f"{num:.3f}" for num in feat_tmp1])

                                    numpy_arr3 = np.array(result_pro_2[Remaining_residues])
                                    numpy_arr4 = np.array(feat.iloc[ID - 1, 296:321]).astype('float64')
                                    feat_tmp2 = np.subtract(numpy_arr3, numpy_arr4).astype('float64')
                                    feat_2 = ','.join([f"{num:.3f}" for num in feat_tmp2])

                                    f.write(
                                        f'{feat_1},{feat_2}')  # {res_name[ID-1]}{ID}{residue_three_to_one_letter[Remaining_residues]}
                                    f.write("\n")
                                else:
                                    data_feat1 = read_csv_and_assign(
                                        f'./All_protein_toAS/data_{model}_toAS1.csv')
                                    data_feat2 = read_csv_and_assign(
                                        f'./All_protein_toAS/data_{model}_toAS2.csv')
                                    numpy_arr1 = np.array(data_feat1[Remaining_residues])
                                    numpy_arr2 = np.array(feat.iloc[ID - 1, 241:295]).astype('float64')
                                    feat_tmp1 = np.subtract(numpy_arr1, numpy_arr2).astype('float64')
                                    feat_1 = ','.join([f"{num:.3f}" for num in feat_tmp1])

                                    numpy_arr3 = np.array(data_feat2[Remaining_residues])
                                    numpy_arr4 = np.array(feat.iloc[ID - 1, 296:321]).astype('float64')
                                    feat_tmp2 = np.subtract(numpy_arr3, numpy_arr4).astype('float64')
                                    feat_2 = ','.join([f"{num:.3f}" for num in feat_tmp2])

                                    f.write(
                                        f'{feat_1},{feat_2}')  # {res_name[ID-1]}{ID}{residue_three_to_one_letter[Remaining_residues]}
                                    f.write("\n")
            print(f'{model}:The mutagenesis characteristics of the protein of {input_file_name} have been calculated!')

        else:
            print(f'     ./FEATURE.feat/{input_file_name}/{model}/{pdb_file}_{model}.csv doesn\'t exist,\n'
                '     This method takes the feature file of the FEATURE and places it under Feature.feat/${your list}/${model}\n'
                '     And you are supposed to rename the file\'s name to ${pdb}_${model}.csv\n\n')