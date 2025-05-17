import os
import shutil
import sys

def convert_pdb_to_FEATURE_input_pdb(input_file_name,pdb_dir):
    alphabet = 'abcdefghijklmnopqrstuvwxyz'
    count=0
    with open(f'{input_file_name}', 'r') as file:
        In_lines=[]
        Notin_lines=[]
        for line in file:
            pdb_file=line.strip()
            file_path=f"{pdb_dir}/{pdb_file}.pdb"
            if os.path.exists(file_path):
                index = count + 1
                letter4 = ((index - 1) // 26 // 26 // 26 % 26)
                letter1 = alphabet[(index - 1) // 26 // 26 % 26]
                letter2 = alphabet[(index - 1) // 26 % 26]
                letter3 = alphabet[(index - 1) % 26]
                new_name = f'{letter4}{letter1}{letter2}{letter3}'

                print(f"{pdb_dir}/{pdb_file}.pdb,{new_name}")
                shutil.copy(file_path, f"./FEATURE_data/pdb_FEATURE/{input_file_name}/")
                os.rename(f"./FEATURE_data/pdb_FEATURE/{input_file_name}/{pdb_file}.pdb", f"./FEATURE_data/pdb_FEATURE/{input_file_name}/{new_name}.pdb")
                In_lines.append(f'{pdb_file},{new_name}\n')
            else:
                Notin_lines.append(f'{pdb_file}\n')
            with open(f'./rename_list/{input_file_name}/In_list.list','w')as f:
                f.writelines(In_lines)
            with open(f'./rename_list/{input_file_name}/Notin_list.list','w') as f:
                f.writelines(Notin_lines)
            count += 1
def rename_Feature_list(input_file_name):
    alphabet = 'abcdefghijklmnopqrstuvwxyz'
    with open(input_file_name, 'r') as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            index = i + 1
            letter4 = ((index - 1) // 26// 26 // 26 % 26)
            letter1 = alphabet[(index - 1) // 26 // 26 % 26]
            letter2 = alphabet[(index - 1) // 26 % 26]
            letter3 = alphabet[(index - 1) % 26]
            new_name = f'{letter4}{letter1}{letter2}{letter3}'
            lines[i] = new_name + '\n'
    with open(f'./rename_list/{input_file_name}/Feature_PDBIDS.list', 'a') as file:
        file.writelines(lines)

def main():
    input_file_name=sys.argv[1]
    pdb_dir = sys.argv[2]

    # create rename_list-dir
    if not os.path.exists('./rename_list'):
        os.mkdir('./rename_list')
    if os.path.exists(f'./rename_list/{input_file_name}'):
        shutil.rmtree(f'./rename_list/{input_file_name}')
    os.mkdir(f'./rename_list/{input_file_name}')

    if os.path.exists(f'./FEATURE_data/pdb_FEATURE/{input_file_name}'):
        shutil.rmtree(f'./FEATURE_data/pdb_FEATURE/{input_file_name}')
    os.mkdir(f'./FEATURE_data/pdb_FEATURE/{input_file_name}')

    if os.path.exists(f'./FEATURE_data/dssp_FEATURE/{input_file_name}'):
        shutil.rmtree(f'./FEATURE_data/dssp_FEATURE/{input_file_name}')
    os.mkdir(f'./FEATURE_data/dssp_FEATURE/{input_file_name}')

    #RUN
    convert_pdb_to_FEATURE_input_pdb(input_file_name,pdb_dir)
    rename_Feature_list(input_file_name)

if __name__ == '__main__':
    main()
