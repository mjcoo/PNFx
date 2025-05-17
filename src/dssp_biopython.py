from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import os
import sys
import shutil
import  warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

res_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO',
            'SER', 'THR', 'TRP', 'TYR', 'VAL']

residue_three_to_one_letter = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def DSSP_run(pdb_file):
    p = PDBParser()
    structure = p.get_structure("protein", pdb_file)
    model = structure[0]
    dssp = DSSP(model, f'{pdb_file}', dssp='mkdssp')
    RES_chain = list(t[0] for t in dssp.keys())
    RES_num = list(f'{t[1][1]}{t[1][2]}' for t in dssp.keys())
    RES = [str(x) + str(y) for x, y in zip(RES_chain, RES_num)]
    Secondary_structure = [res[2] for res in dssp]
    RASA = [res[3] for res in dssp]
    RES_name=[residue.get_resname() for residue in structure.get_residues()]
    RES_dssp_dict = {item1: [item2, item3,item4] for item1, item2, item3,item4 in zip(RES, Secondary_structure, RASA,RES_name)}
    return RES_dssp_dict

def feat_process(RES_dssp_dict,output_file):
    with open(output_file, 'a') as f:
        f.write(
            "Secondary_structure_is_H,Secondary_structure_is_B,Secondary_structure_is_I,Secondary_structure_is_G,Secondary_structure_is_T,Secondary_structure_is_S,Secondary_structure_is_E,Secondary_structure_is_C,res_RASA")
        f.write('\n')
        for key in RES_dssp_dict.keys():
            for Remaining_residues in res_list:
                if Remaining_residues != RES_dssp_dict[key][2]:
                    H, B, I, G, T, S, E, C = 0, 0, 0, 0, 0, 0, 0, 0
                    if RES_dssp_dict[key][0] == 'H':
                        H = 1
                    if RES_dssp_dict[key][0] == 'B':
                        B = 1
                    if RES_dssp_dict[key][0] == 'I':
                        I = 1
                    if RES_dssp_dict[key][0] == 'G':
                        G = 1
                    if RES_dssp_dict[key][0] == 'T':
                        T = 1
                    if RES_dssp_dict[key][0] == 'S':
                        S = 1
                    if RES_dssp_dict[key][0] == 'E':
                        E = 1
                    if RES_dssp_dict[key][0] == '-':
                        C = 1
                    f.write(f'{H},{B},{I},{G},{T},{S},{E},{C},{RES_dssp_dict[key][1]}')
                    f.write("\n")



def main():
    input_file_name = sys.argv[1]
    pdb_dir = sys.argv[2]

    if not os.path.exists('dssp_biopython.feat'):
        os.mkdir('dssp_biopython.feat')
    if not os.path.exists(f'dssp_biopython.feat/{input_file_name}'):
        os.mkdir(f'dssp_biopython.feat/{input_file_name}')
    else:
        shutil.rmtree(f'dssp_biopython.feat/{input_file_name}')
        os.mkdir(f'dssp_biopython.feat/{input_file_name}')
    print(f'\ndssp_diopython.py:Start to calculate {input_file_name}\'s RSA and the One-Hot of secondary structure')
    with open(f"{input_file_name}", 'r') as file:
        for line in file:
            pdb_file = line.strip()
            file_path = f"{pdb_dir}/{pdb_file}.pdb"
            if os.path.exists(file_path):
                DSSP_data=DSSP_run(file_path)
                feat_process(DSSP_data,f'./dssp_biopython.feat/{input_file_name}/{pdb_file}.feat')

if __name__ == "__main__":
    main()