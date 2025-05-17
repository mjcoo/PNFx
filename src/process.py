import shutil
import sys
import os
import traceback


res_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO',
            'SER', 'THR', 'TRP', 'TYR', 'VAL']

def process_pdb(pdb_file, out_file):
    '''
    FIRST STEP:rename the seq of the res to standard form.
    :param pdb_file:input pdb file path
    :param out_file:output pdb file path
    :return:
    '''
    with open(f"{pdb_file}", 'r') as f:
        lines = f.readlines()
    chain_NUM_dict = {}
    resseq_list = {}
    new_lines = []
    chain_list=[]
    for line in lines:
        if line.startswith('ATOM'):
            chainID = line[21]
            number = line[22:26]
            if chainID not in resseq_list:
                resseq_list[chainID] = []
                chain_list.append(chainID)
            if number not in resseq_list[chainID]:
                resseq_list[chainID].append(int(number))
            chain_NUM_dict[chainID] = min(resseq_list[chainID])

    for line in lines:
        if line.startswith('ATOM'):
            chainID = line[21]
            if chain_NUM_dict[chainID] > 0:
                ATOM=line[0:4]
                serial=line[6:11].rjust(5)
                name=line[12:16].ljust(4)
                altLoc=line[16]
                resName=line[17:20]
                chainID=line[21]
                resSeq=line[22:26].rjust(4)
                iCode=line[26]
                x = "{:>8.3f}".format(float(line[30:38]))
                y = "{:>8.3f}".format(float(line[38:46]))
                z = "{:>8.3f}".format(float(line[46:54]))
                occupancy=line[54:60].rjust(6)
                tempFactor=line[60:66].rjust(6)
                segID=line[72:76].ljust(4)
                element=line[76:78].rjust(2)
                charge=line[78:80]
                result = "".join(
                    [ATOM,"  ",serial," ",name, altLoc, resName," ",chainID, resSeq, iCode,"   ", x, y, z, occupancy, tempFactor,"      ",segID,
                     element, charge])
                new_lines.append(result + '\n')

            if chain_NUM_dict[chainID] <= 0:
                add = abs(chain_NUM_dict[chainID]) + 1
                ATOM=line[0:4]
                serial=line[6:11].rjust(5)
                name=line[12:16].ljust(4)
                altLoc=line[16]
                resName=line[17:20]
                chainID=line[21]
                iCode=line[26]
                x = "{:>8.3f}".format(float(line[30:38]))
                y = "{:>8.3f}".format(float(line[38:46]))
                z = "{:>8.3f}".format(float(line[46:54]))
                occupancy=line[54:60].rjust(6)
                tempFactor=line[60:66].rjust(6)
                segID=line[72:76].ljust(4)
                element=line[76:78].rjust(2)
                charge=line[78:80]
                resSeq=line[22:26]
                number_refactor = str(int(resSeq) + add).rjust(4)
                result = "".join(
                    [ATOM,"  ",serial," ",name, altLoc, resName," ",chainID, number_refactor, iCode,"   ", x, y, z, occupancy, tempFactor,"      ",segID,
                     element, charge])
                new_lines.append(result + '\n')


    '''
    SECOND STEP:process pdb:add TER and END
    '''
    new_lines2 = []
    chain = {}
    chain_for_write = {}
    chain_num = 0
    for line in new_lines:
        if line.startswith('ATOM'):
            if line[21] not in chain.keys():
                chain[line[21]] = 1
            else:
                chain[line[21]] += 1
    for line in new_lines:
        if line.startswith('ATOM'):
            new_lines2.append(line)
            if line[21] not in chain_for_write:
                chain_for_write[line[21]] = 1
            else:
                chain_for_write[line[21]] += 1
            if line[21] in chain_for_write and chain_for_write[line[21]] == chain[line[21]]:
                new_lines2.append('TER\n')
                chain_num += 1

    new_lines2.append('END\n')
    with open(out_file, 'w') as f:
        f.writelines(new_lines2)

def fix(pdb_file):
    from Pras_Server.RunType import InitRunType

    rotamer = ''
    mutation = ''
    pdb_faspr = ''
    keep_ligand = ''
    chain_no = ''
    ofname = ''
    his_p = ''
    fixing = InitRunType(rotamer, mutation, pdb_faspr, keep_ligand, chain_no,
                         addh=False, ss=False, raman=False, ofname=False, pdbid=False, his_p=False)
    fixing.fname = pdb_file
    fixing.ProcessWithDefault()
    os.rename('out.pdb', f'fix_{fixing.fname}')

def main():

    input_file_name = sys.argv[1]
    pdb_dir = sys.argv[2]


    if os.path.exists(f'{input_file_name}'):

        # create fix_pdb
        if not os.path.exists('./fix_pdb'):
            os.mkdir('./fix_pdb')
        if not os.path.exists(f'./fix_pdb/{input_file_name}'):
            os.mkdir(f'./fix_pdb/{input_file_name}')
        else:
            shutil.rmtree(f'./fix_pdb/{input_file_name}')
            os.mkdir(f'./fix_pdb/{input_file_name}')
        #RUN proccess and fix
        with open(f'{input_file_name}', 'r') as file:
            for line in file:
                try:
                    pdb_file = line.strip()
                    file_path = f"{pdb_dir}/{pdb_file}.pdb"
                    if os.path.exists(file_path):
                        shutil.copy(file_path,f'./fix_pdb/{input_file_name}/{pdb_file}.pdb')
                        os.chdir(f'./fix_pdb/{input_file_name}')
                        process_pdb(f'{pdb_file}.pdb', f'{pdb_file}.pdb')
                        os.rename(f'{pdb_file}.pdb',f'tmp.pdb')
                        fix('tmp.pdb')
                        os.remove('tmp.pdb')
                        os.rename(f'fix_tmp.pdb', f'{pdb_file}.pdb')
                        os.remove('log.txt')
                        os.chdir(f'../../')
                    else:
                        print(f'{pdb_dir}/{pdb_file}.pdb dosen\'t exist ')
                        with open(f'Not_exist_pdb.list', 'a') as out:
                            out.write(f'{pdb_dir}/{pdb_file}.pdb dosen\'t exist ')
                except Exception as e:
                    with open(f'../../fix_wrong.log',
                              'a') as wrong_file:
                        wrong_file.write(pdb_file + '\n')
                        wrong_file.write(str(traceback.format_exc()))
                        os.chdir(f'../../')
                    continue
    else:
        print(f'{input_file_name} doesn\'t exists')


if __name__ == '__main__':
    main()