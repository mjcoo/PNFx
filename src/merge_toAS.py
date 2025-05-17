import os
import sys
import shutil
def merge_toAS(merge_dir,FEATURE_dir_toAS,input_list_name,pdb_dir):
    if os.path.exists(f'{input_list_name}'):


        with open(f'{input_list_name}', 'r') as file:
            PDBs_processed = 0
            for line in file:
                pdb_file = line.strip()
                if os.path.exists(f'{pdb_dir}/{input_list_name}/{pdb_file}.pdb'):
                    merge_file= f'{merge_dir}/{input_list_name}/{pdb_file}.feat'
                    FEATURE_toAS_file1 = f'{FEATURE_dir_toAS}/{input_list_name}/6A_ball/{pdb_file}_6A_ball_toAS.feat'  # Feature_6A_ball
                    FEATURE_toAS_file2 = f'{FEATURE_dir_toAS}/{input_list_name}/4.5A_ball/{pdb_file}_4.5A_ball_toAS.feat'  # Feature_4.5A_ball
                    FEATURE_toAS_file3 = f'{FEATURE_dir_toAS}/{input_list_name}/4.5A_to_6A_shell/{pdb_file}_4.5A_to_6A_shell_toAS.feat'  # Feature_4.5A_to_6A_shell
                    if os.path.exists(merge_file) and os.path.exists(FEATURE_toAS_file1) and os.path.exists(
                            FEATURE_toAS_file2) and os.path.exists(FEATURE_toAS_file3):
                        with (open(merge_file, 'r') as f1, open(FEATURE_toAS_file1, 'r') as f2, open(FEATURE_toAS_file2, 'r') as f3, open(FEATURE_toAS_file3,'r') as f4,
                              open(f'./FEAT/{input_list_name}/{pdb_file}.feat',
                                                             'w') as out):
                            lines1 = f1.readlines()
                            lines2 = f2.readlines()
                            lines3 = f3.readlines()
                            lines4 = f4.readlines()
                            for l1, l2, l3, l4 in zip(lines1, lines2, lines3, lines4):
                                out.write(l1.strip() + ',' + l2.strip() + ',' + l3.strip() + ',' + l4.strip()+'\n')
                            PDBs_processed += 1
                    else:
                        if os.path.exists(merge_file) == False:
                            print(merge_file+' doesn\'t exist')
                        if os.path.exists(FEATURE_toAS_file1) == False:
                            print(FEATURE_toAS_file1+' doesn\'t exist')
                        if os.path.exists(FEATURE_toAS_file2) == False:
                            print(FEATURE_toAS_file2+' doesn\'t exist')
                        if os.path.exists(FEATURE_toAS_file3) == False:
                            print(FEATURE_toAS_file3+' doesn\'t exist')
                else:
                    print(f'\nmerge_toAS.py:{pdb_dir}/{input_list_name}/{pdb_file}.pdb doesn\'t exist')

        print(f'\nmerge_toAS.py:PDBs_processed:{PDBs_processed}!!\n\n')
    else:
        print(f'\nmerge_toAS.py:{input_list_name} doesn\'t exist')


def main():
    input_file_name = sys.argv[1]
    merge_dir=sys.argv[2]
    FEATURE_dir_toAS= sys.argv[3]
    pdb_dir=sys.argv[4]
    # prepare
    if not os.path.exists('./FEAT'):
        os.mkdir('./FEAT')
    if os.path.exists(f'./FEAT/{input_file_name}'):
        shutil.rmtree(f'./FEAT/{input_file_name}')
        os.mkdir(f'./FEAT/{input_file_name}')
    else:
        os.mkdir(f'./FEAT/{input_file_name}')
    merge_toAS(merge_dir,FEATURE_dir_toAS,input_file_name,pdb_dir)

if __name__ == '__main__':
    main()

