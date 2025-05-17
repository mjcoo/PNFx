import os
import sys
import shutil
import pandas as pd

res_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO',
            'SER', 'THR', 'TRP', 'TYR', 'VAL']

residue_three_to_one_letter = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

''' 19ci
def get_position(input_file_name, pdb_dir):
    if os.path.exists(f'{input_file_name}'):
        with open(f'{input_file_name}', 'r') as file:
            for line in file:
                pdb_file = line.strip()
                file_path = f"{pdb_dir}/{pdb_file}.pdb"
                if os.path.exists(file_path):
                    new_line = []
                    chain_res_seq = {}
                    with open(file_path, 'r') as f:
                        lines = f.readlines()
                        for line in lines:
                            if line.startswith('ATOM'):
                                resName = line[17:20].strip()
                                resseq = line[22:26].strip()
                                chainID = line[21]
                                iCode = line[26]
                                id_iCode = (resseq + iCode).strip()
                                if chainID not in chain_res_seq.keys():
                                    chain_res_seq[chainID] = []
                                if id_iCode not in chain_res_seq[chainID]:
                                    chain_res_seq[chainID].append(id_iCode)
                                    for Remaining_residues in res_list:
                                        if Remaining_residues != resName:
                                            new_line.append(
                                                f"{residue_three_to_one_letter[resName]},{id_iCode},{residue_three_to_one_letter[Remaining_residues]},{chainID}\n")
                    with open(f'./position/{input_file_name}/{pdb_file}.pos', 'w') as out:
                        out.writelines(new_line)
                    with open(f'./position/{input_file_name}/{pdb_file}.pos', 'r') as f:
                        original_lines = f.readlines()
                    original_lines.insert(0, "Residue_name,position,Alt_Residue_name,chain\n")
                    with open(f'./position/{input_file_name}/{pdb_file}.pos', 'w') as out:
                        out.writelines(original_lines)
    else:
        print(f'{input_file_name} doesn\'t exist')


def merge(pos_dir, hse_dir, dssp_dir, FEATURE_dir, input_list_name, pdb_dir):
    if os.path.exists(f'{input_list_name}'):

        with open(f'{input_list_name}', 'r') as file:
            PDBs_processed = 0
            for line in file:
                pdb_file = line.strip()
                if os.path.exists(f'{pdb_dir}/{pdb_file}.pdb'):
                    pos_file = f'{pos_dir}/{input_list_name}/{pdb_file}.pos'  # position
                    hse_file = f'{hse_dir}/{input_list_name}/{pdb_file}.feat'  # hse_feat
                    dssp_file = f'{dssp_dir}/{input_list_name}/{pdb_file}.feat'  # dssp_feat
                    FEATURE_file1 = f'{FEATURE_dir}/{input_list_name}/6A_ball/{pdb_file}_feature_6A_ball.feat'  # Feature_6A_ball
                    FEATURE_file2 = f'{FEATURE_dir}/{input_list_name}/4.5A_ball/{pdb_file}_feature_4.5A_ball.feat'  # Feature_4.5A_ball
                    FEATURE_file3 = f'{FEATURE_dir}/{input_list_name}/4.5A_to_6A_shell/{pdb_file}_feature_4.5A_to_6A_shell.feat'  # Feature_4.5A_to_6A_shell
                    if os.path.exists(pos_file) and os.path.exists(hse_file) and os.path.exists(
                            dssp_file) and os.path.exists(FEATURE_file1) and os.path.exists(
                        FEATURE_file2) and os.path.exists(FEATURE_file3):
                        with (open(pos_file, 'r') as f1, open(hse_file, 'r') as f2, open(dssp_file, 'r') as f3, open(
                                FEATURE_file1, 'r') as f4,
                              open(FEATURE_file2, 'r') as f5, open(FEATURE_file3, 'r') as f6, open(
                            f'./feat/{input_list_name}/{pdb_file}.feat',
                            'w') as out):
                            lines1 = f1.readlines()
                            lines2 = f2.readlines()
                            lines3 = f3.readlines()
                            lines4 = f4.readlines()
                            lines5 = f5.readlines()
                            lines6 = f6.readlines()
                            for l1, l2, l3, l4, l5, l6 in zip(lines1, lines2, lines3, lines4, lines5, lines6):
                                out.write(
                                    l1.strip() + ',' + l2.strip() + ',' + l3.strip() + ',' + l4.strip() + ',' + l5.strip() + ',' + l6.strip() + '\n')
                            PDBs_processed += 1
                    else:
                        if os.path.exists(pos_file) == False:
                            print(pos_file + ' doesn\'t exist')
                        if os.path.exists(hse_file) == False:
                            print(hse_file + ' doesn\'t exist')
                        if os.path.exists(dssp_file) == False:
                            print(dssp_file + ' doesn\'t exist')
                        if os.path.exists(FEATURE_file1) == False:
                            print(FEATURE_file1 + ' doesn\'t exist')
                        if os.path.exists(FEATURE_file2) == False:
                            print(FEATURE_file2 + ' doesn\'t exist')
                        if os.path.exists(FEATURE_file3) == False:
                            print(FEATURE_file3 + ' doesn\'t exist')
                else:
                    print(f'\nmerge.py:{pdb_dir}/{pdb_file}.pdb doesn\'t exist\n')

        print(f'\nmerge.py:PDBs_processed:{PDBs_processed}!!!!\n')
    else:
        print(f'\nmerge.py:{input_list_name} doesn\'t exist')
'''


def get_position(input_file_name, pdb_dir):
    if os.path.exists(f'{input_file_name}'):
        with open(f'{input_file_name}', 'r') as file:
            for line in file:
                pdb_file = line.strip()
                file_path = f"{pdb_dir}/{pdb_file}.pdb"
                if os.path.exists(file_path):
                    new_line = []
                    chain_res_seq = {}
                    with open(file_path, 'r') as f:
                        lines = f.readlines()
                        for line in lines:
                            if line.startswith('ATOM'):
                                resName = line[17:20].strip()
                                resseq = line[22:26].strip()
                                chainID = line[21]
                                iCode = line[26]
                                id_iCode = (resseq + iCode).strip()
                                if chainID not in chain_res_seq.keys():
                                    chain_res_seq[chainID] = []
                                if id_iCode not in chain_res_seq[chainID]:
                                    chain_res_seq[chainID].append(id_iCode)
                                    new_line.append(f"{residue_three_to_one_letter[resName]},{id_iCode},{chainID}\n")
                    with open(f'./position/{input_file_name}/{pdb_file}.pos', 'w') as out:
                        out.writelines(new_line)
                    with open(f'./position/{input_file_name}/{pdb_file}.pos', 'r') as f:
                        original_lines = f.readlines()
                    original_lines.insert(0, "Residue_name,position,chain\n")
                    with open(f'./position/{input_file_name}/{pdb_file}.pos', 'w') as out:
                        out.writelines(original_lines)
    else:
        print(f'{input_file_name} doesn\'t exist')


def merge(pos_dir, hse_dir, dssp_dir, FEATURE_dir, input_list_name, pdb_dir):
    if os.path.exists(f'{input_list_name}'):

        with open(f'{input_list_name}', 'r') as file:
            PDBs_processed = 0
            for line in file:
                pdb_file = line.strip()
                if os.path.exists(f'{pdb_dir}/{pdb_file}.pdb'):
                    pos_file = f'{pos_dir}/{input_list_name}/{pdb_file}.pos'  # position

                    hse_file = f'{hse_dir}/{input_list_name}/{pdb_file}.feat'  # hse_feat
                    hse = pd.read_csv(hse_file, header=None, low_memory = False)
                    hse_selected_lines = pd.concat([hse.iloc[[0]], hse.iloc[1::19]])
                    hse_selected_lines.to_csv(f'{hse_dir}/{input_list_name}/{pdb_file}_2.feat',index=False,header=False)
                    hse_file_2=f'{hse_dir}/{input_list_name}/{pdb_file}_2.feat'

                    dssp_file = f'{dssp_dir}/{input_list_name}/{pdb_file}.feat'  # dssp_feat
                    dssp = pd.read_csv(dssp_file, header=None, low_memory = False)
                    dssp_selected_lines = pd.concat([dssp.iloc[[0]], dssp.iloc[1::19]])
                    dssp_selected_lines.to_csv(f'{dssp_dir}/{input_list_name}/{pdb_file}_2.feat',index=False,header=False)
                    dssp_file_2 = f'{dssp_dir}/{input_list_name}/{pdb_file}_2.feat'

                    FEATURE_file1 = f'{FEATURE_dir}/{input_list_name}/6A_ball/{pdb_file}_feature_6A_ball.feat'  # Feature_6A_ball
                    FEATURE_file2 = f'{FEATURE_dir}/{input_list_name}/4.5A_ball/{pdb_file}_feature_4.5A_ball.feat'  # Feature_4.5A_ball
                    FEATURE_file3 = f'{FEATURE_dir}/{input_list_name}/4.5A_to_6A_shell/{pdb_file}_feature_4.5A_to_6A_shell.feat'  # Feature_4.5A_to_6A_shell
                    if os.path.exists(pos_file) and os.path.exists(hse_file_2) and os.path.exists(dssp_file_2)and os.path.exists(FEATURE_file1) and os.path.exists(FEATURE_file2) and os.path.exists(FEATURE_file3):
                        with (open(pos_file, 'r') as f1,open(hse_file_2,'r')as f2,open(dssp_file_2,'r')as f3, open(f'./feat/{input_list_name}/{pdb_file}.feat','w') as out,open(FEATURE_file1, 'r') as f4, open(FEATURE_file2, 'r') as f5, open(FEATURE_file3, 'r') as f6):
                            lines1 = f1.readlines()
                            lines2 = f2.readlines()
                            lines3 = f3.readlines()
                            lines4 = f4.readlines()
                            lines5 = f5.readlines()
                            lines6 = f6.readlines()
                            for l1, l2, l3, l4, l5, l6 in zip(lines1, lines2, lines3, lines4, lines5, lines6):
                                out.write(
                                    l1.strip() + ',' + l2.strip() + ',' + l3.strip() + ',' + l4.strip() + ',' + l5.strip() + ',' + l6.strip() + '\n')
                            PDBs_processed += 1
                    else:
                        if os.path.exists(pos_file) == False:
                            print(pos_file + ' doesn\'t exist')
                        if os.path.exists(hse_file_2) == False:
                            print(hse_file + ' doesn\'t exist')
                        if os.path.exists(dssp_file_2) == False:
                            print(dssp_file + ' doesn\'t exist')
                        if os.path.exists(FEATURE_file1) == False:
                            print(FEATURE_file1 + ' doesn\'t exist')
                        if os.path.exists(FEATURE_file2) == False:
                            print(FEATURE_file2 + ' doesn\'t exist')
                        if os.path.exists(FEATURE_file3) == False:
                            print(FEATURE_file3 + ' doesn\'t exist')
                else:
                    print(f'\nmerge.py:{pdb_dir}/{pdb_file}.pdb doesn\'t exist\n')

        print(f'\nmerge.py:PDBs_processed:{PDBs_processed}!!!!\n')
    else:
        print(f'\nmerge.py:{input_list_name} doesn\'t exist')


def main():
    input_file_name = sys.argv[1]
    hse_dir = sys.argv[2]
    dssp_dir = sys.argv[3]
    FEATURE_dir = sys.argv[4]
    pdb_dir = sys.argv[5]
    pos_dir = './position'

    # prepare position dir
    if not os.path.exists('./position'):
        os.mkdir('./position')
    if os.path.exists(f'./position/{input_file_name}'):
        shutil.rmtree(f'./position/{input_file_name}')
        os.mkdir(f'./position/{input_file_name}')
    else:
        os.mkdir(f'./position/{input_file_name}')

    if not os.path.exists('./feat'):
        os.mkdir('./feat')
    if os.path.exists(f'./feat/{input_file_name}'):
        shutil.rmtree(f'./feat/{input_file_name}')
        os.mkdir(f'./feat/{input_file_name}')
    else:
        os.mkdir(f'./feat/{input_file_name}')

    get_position(input_file_name, pdb_dir)
    merge(pos_dir, hse_dir, dssp_dir, FEATURE_dir, input_file_name, pdb_dir)


if __name__ == '__main__':
    main()
