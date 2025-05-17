import pandas as pd
import os
import sys
from Bio.PDB import PDBParser
import  warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

feature_dir = os.environ.get('FEATURE_DIR')
input_file_name = sys.argv[1]
model = sys.argv[2]
pdb_dir = sys.argv[3]

res_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO',
            'SER', 'THR', 'TRP', 'TYR', 'VAL']

residue_three_to_one_letter = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}
'''19ci
with open(f'{input_file_name}', 'r') as file:
    for line in file:
        pdb_file = line.strip()
        if os.path.exists(f'./FEATURE.feat/{input_file_name}/{model}/{pdb_file}_ca_{model}.ff'):
            with open(f'./FEATURE.feat/{input_file_name}/{model}/{pdb_file}_ca_{model}.ff', 'r') as infile, open(
                    f'./FEATURE.feat/{input_file_name}/{model}/{pdb_file}_{model}.csv', 'w') as outfile:
                for i in range(36):
                    next(infile)
                for line in infile:
                    outfile.write(line.replace('\t', ','))
        else:
            continue

feature_list = pd.read_csv(f'{feature_dir}/Feature_feature_{model}_list.csv', header=None)
pro = ','.join(feature_list.iloc[:, 0].astype(str))

with open(f'{input_file_name}', 'r') as file:
    for line in file:
        try:
            pdb_file = line.strip()
            if os.path.exists(f'./FEATURE.feat/{input_file_name}/{model}/{pdb_file}_{model}.csv'):
                p = PDBParser()
                file_path = f"{pdb_dir}/{pdb_file}.pdb"
                structure = p.get_structure("protein", file_path)
                residues = structure.get_residues()

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
                if model == '6A_ball' or model == '4.5A_ball':
                    with open(f'./FEATURE.feat/{input_file_name}/{model}/{pdb_file}_feature_{model}.feat', "a") as f:
                        f.write(f"{pro}")
                        f.write("\n")
                        for ID in range(1, len(res_id) + 1):
                            for Remaining_residues in res_list:
                                if Res[ID - 1] != Remaining_residues:
                                    feat_1 = ','.join(feat.iloc[ID - 1, 1:55].astype(str))
                                    feat_2 = ','.join(feat.iloc[ID - 1, 56:81].astype(str))
                                    f.write(
                                        f'{feat_1},{feat_2}')  # {res_name[ID-1]}{ID}{residue_three_to_one_letter[Remaining_residues]}
                                    f.write("\n")
                    os.remove(f'./FEATURE.feat/{input_file_name}/{model}/{pdb_file}_ca_{model}.ptf')
                    os.remove(f'./FEATURE.feat/{input_file_name}/{model}/{pdb_file}_ca_{model}.ff')

                if model == '4.5A_to_6A_shell':
                    with open(f'./FEATURE.feat/{input_file_name}/{model}/{pdb_file}_feature_{model}.feat', "a") as f:
                        f.write(f"{pro}")
                        f.write("\n")
                        for ID in range(1, len(res_id) + 1):
                            for Remaining_residues in res_list:
                                if Res[ID - 1] != Remaining_residues:
                                    feat_1 = ','.join(feat.iloc[ID - 1, 241:295].astype(str))
                                    feat_2 = ','.join(feat.iloc[ID - 1, 296:321].astype(str))
                                    f.write(
                                        f'{feat_1},{feat_2}')  # {res_name[ID-1]}{ID}{residue_three_to_one_letter[Remaining_residues]}
                                    f.write("\n")
                    os.remove(f'./FEATURE.feat/{input_file_name}/{model}/{pdb_file}_ca_{model}.ptf')
                    os.remove(f'./FEATURE.feat/{input_file_name}/{model}/{pdb_file}_ca_{model}.ff')
'''
with open(f'{input_file_name}', 'r') as file:
    for line in file:
        pdb_file = line.strip()
        if os.path.exists(f'./FEATURE.feat/{input_file_name}/{model}/{pdb_file}_ca_{model}.ff'):
            with open(f'./FEATURE.feat/{input_file_name}/{model}/{pdb_file}_ca_{model}.ff', 'r') as infile, open(
                    f'./FEATURE.feat/{input_file_name}/{model}/{pdb_file}_{model}.csv', 'w') as outfile:
                for i in range(36):
                    next(infile)
                for line in infile:
                    outfile.write(line.replace('\t', ','))
        else:
            continue

feature_list = pd.read_csv(f'{feature_dir}/Feature_feature_{model}_list.csv', header=None)
pro = ','.join(feature_list.iloc[:, 0].astype(str))

with open(f'{input_file_name}', 'r') as file:
    for line in file:
        try:
            pdb_file = line.strip()
            if os.path.exists(f'./FEATURE.feat/{input_file_name}/{model}/{pdb_file}_{model}.csv'):
                p = PDBParser()
                file_path = f"{pdb_dir}/{pdb_file}.pdb"
                structure = p.get_structure("protein", file_path)
                residues = structure.get_residues()

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
                if model == '6A_ball' or model == '4.5A_ball':
                    with open(f'./FEATURE.feat/{input_file_name}/{model}/{pdb_file}_feature_{model}.feat', "a") as f:
                        f.write(f"{pro}")
                        f.write("\n")
                        for ID in range(1, len(res_id) + 1):
                            feat_1 = ','.join(feat.iloc[ID - 1, 1:55].astype(str))
                            feat_2 = ','.join(feat.iloc[ID - 1, 56:81].astype(str))
                            f.write(
                                f'{feat_1},{feat_2}')  # {res_name[ID-1]}{ID}{residue_three_to_one_letter[Remaining_residues]}
                            f.write("\n")
                    
                if model == '4.5A_to_6A_shell':
                    with open(f'./FEATURE.feat/{input_file_name}/{model}/{pdb_file}_feature_{model}.feat', "a") as f:
                        f.write(f"{pro}")
                        f.write("\n")
                        for ID in range(1, len(res_id) + 1):
                            feat_1 = ','.join(feat.iloc[ID - 1, 241:295].astype(str))
                            feat_2 = ','.join(feat.iloc[ID - 1, 296:321].astype(str))
                            f.write(
                                f'{feat_1},{feat_2}')  # {res_name[ID-1]}{ID}{residue_three_to_one_letter[Remaining_residues]}
                            f.write("\n")
                    
        except Exception as e:
            continue
