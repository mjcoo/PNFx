import argparse
import os
import sys

__author__ = "Junwen Huang"
__version__ = "v1.0"


def get_args():
    parser = argparse.ArgumentParser(
        prog='PNFx',
        description='Compute different structure features in protein stored in pdb files or pdb list.')
    parser.add_argument(
        '-l',
        type=str,
        help='Specify the protein list file')
    parser.add_argument(
        '-p',
        type=str,
        help='Specify the protein')
    parser.add_argument(
        '-toAS',
        type=str,
        choices=['0', '1'],
        default='0',
        help='Whether to use the toAS feature for the current list, 0 or 1, and 1 represent the calculation of the toAS feature')
    parser.add_argument(
        '-pdbdir',
        type=str,
        help='Specify the PDB data directory')
    args = parser.parse_args()
    error_messages = []
    if not args.l and not args.p and not args.toAS == 0 and not args.pdbdir:
        error_messages.append('Error: -l or -p must be specified,-toAS and -pdbdir must be both specified!\n')
    if not args.pdbdir:
        error_messages.append('Error: -pdbdir must be specified!\n')
    else:
        if not os.path.isdir(args.pdbdir):
            error_messages.append('Error: -pdbdir must be an existing directory!\n')
    if args.l and args.p:
        error_messages.append('Error: -l and -p cannot be specified at the same time!\n')
    if args.p:
        if os.path.splitext(args.p)[1] != '.pdb':
            error_messages.append(
                'Error: Please select a protein pdb file with the suffix ".pdb" for the -p parameter!\n')
    print('......................................................................')
    print('.                                                                    .')
    print('.                              PNFX                                  .')
    print('.      Compute different structure features in protein stored in pdb .')
    print('.      files or pdb list                                             .')
    print('.                                                                    .')
    print('. -l:                                                                .')
    print('.    Specify the protein list file                                   .')
    print('. -p:                                                                .')
    print('.    Specify the protein                                             .')
    print('. -toAS:                                                             .')
    print('.    (1)Whether to use the toAS feature for the current list, 0 or 1 .')
    print('.     and 1 represent the calculation of the toAS feature(default:0) .')
    print('.    (2)toAS refers to the microenvironment features obtained by     .')
    print('.     FEATURE minus the average value of each residue in a single    .')
    print('.     protein or protein list                                        .')
    print('. -pdbdir:                                                           .')
    print('.    Specify the PDB data directory                                  .')
    print('.                                                                    .')
    print('......................................................................')
    if len(error_messages) != 0:
        for error in error_messages:
            print(error)
            sys.exit(1)
    return args.l, args.p, args.toAS, args.pdbdir


if __name__ == "__main__":
    l, p, toAS, pdbdir = get_args()
    try:
        python2bin = os.environ.get('python2bin')
        if p is None and l is not None:
            cmd = f"bash workflow.sh {l} {pdbdir} {python2bin}"
            print(cmd)
            os.system(cmd)
            if toAS == '1':
                cmd = f"bash workflow_toAS.sh {l} ./fix_pdb "
                print(cmd)
                os.system(cmd)
        if p is not None and l is None:
            pdb_base_name = os.path.splitext(p)[0]
            with open(f'{pdb_base_name}.list', 'w') as input_pro:
                input_pro.write(f'{pdb_base_name}')
            cmd = f"bash workflow.sh {pdb_base_name}.list {pdbdir} {python2bin}"
            print(cmd)
            os.system(cmd)
            if toAS == '1':
                cmd = f"bash workflow_toAS.sh {pdb_base_name}.list ./fix_pdb "
                print(cmd)
                os.system(cmd)
    except Exception as e:
        print('Make sure your python2bin is already an environment variable')


