from Bio.PDB import PDBParser
import os
import sys
import traceback
from math import pi
from Bio.PDB.AbstractPropertyMap import AbstractPropertyMap
from Bio.PDB.Polypeptide import is_aa
import biographs as bg
import networkx as nx
from math import acos
import shutil
import warnings
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)


class _AbstractHSExposure2(AbstractPropertyMap):

    def __init__(self, model, radius, offset, hse_up_key, hse_down_key, angle_key=None):

        assert offset >= 0
        self.ca_cb_list = []
        chains = model.get_chains()
        hse_map = {}
        hse_list = []
        hse_keys = []
        for chain1 in chains:
            chain1_seqlist = []
            for residue in chain1:
                if residue.get_id()[2] != ' ':
                    chain1_seqlist.append(f'{residue.get_id()[1]}{residue.get_id()[2]}')
                else:
                    chain1_seqlist.append(residue.get_id()[1])
            for i in chain1_seqlist:
                position1 = ""
                add_code1 = ""
                for c in str(i):
                    if c.isdigit():
                        position1 += c
                    elif c.isalpha():
                        add_code1 += c
                if add_code1 == "":
                    add_code1 += " "
                r1 = chain1[" ", int(position1), add_code1]
                pcb = self._get_cb(r1)
                hse_u = 0
                hse_d = 0
                ca1 = r1["CA"].get_vector()
                ca1_coord = r1["CA"]
                chains = model.get_chains()
                for chain2 in chains:
                    chain2_seqlist = []
                    for residue in chain2:
                        if residue.get_id()[2] != ' ':
                            chain2_seqlist.append(f'{residue.get_id()[1]}{residue.get_id()[2]}')
                        else:
                            chain2_seqlist.append(residue.get_id()[1])
                    for j in chain2_seqlist:
                        position2 = ""
                        add_code2 = ""
                        for c in str(j):
                            if c.isdigit():
                                position2 += c
                            elif c.isalpha():
                                add_code2 += c
                        if add_code2 == "":
                            add_code2 += " "
                        r2 = chain2[" ", int(position2), add_code2]
                        if not is_aa(r2) or not r2.has_id("CA"):
                            continue
                        ca2 = r2["CA"].get_vector()
                        ca2_coord = r2["CA"]
                        distance = ca1_coord - ca2_coord
                        if chain1 is chain2 and distance <= offset:
                            continue
                        d = ca1 - ca2
                        if d.norm() < radius:
                            if d.angle(pcb) < (pi / 2):
                                hse_u += 1
                            else:
                                hse_d += 1
                res_id = r1.get_id()
                chain_id = r1.get_parent().get_id()
                hse_map[(chain_id, res_id)] = (hse_u, hse_d)
                hse_list.append((r1, (hse_u, hse_d)))
                hse_keys.append((chain_id, res_id))
                r1.xtra[hse_up_key] = hse_u
                r1.xtra[hse_down_key] = hse_d
        AbstractPropertyMap.__init__(self, hse_map, hse_keys, hse_list)

    def _get_cb(self, r2):
        return NotImplemented

    def _get_gly_cb_vector(self, residue):
        """
        把甘氨酸的结构近似为成四面体，但是侧链是H原子，和中心氢原子各占据一个角，无法准确定义哪个H原子为侧链方向，则取他们的角平分线的方向向量
        """
        try:
            n_v = residue["N"].get_vector()
            c_v = residue["C"].get_vector()
            ca_v = residue["CA"].get_vector()
        except Exception:
            return None
        n_v = n_v - ca_v
        c_v = c_v - ca_v
        n_v_norm = n_v.norm()
        c_v_norm = c_v.norm()
        n_unit = n_v / n_v_norm if n_v_norm != 0 else n_v  # 单位向量
        c_unit = c_v / c_v_norm if c_v_norm != 0 else c_v  # 单位向量
        bisector = n_unit + c_unit
        bisector_norm = bisector.norm()
        bisector_unit = bisector / bisector_norm if bisector_norm != 0 else bisector  # 计算角平分线向量方向向量
        return bisector_unit


class HSExposureCB(_AbstractHSExposure2):
    """Class to calculate HSE based on the real CA-CB vectors."""

    def __init__(self, model, hse_up_key, hse_down_key, radius, offset):
        _AbstractHSExposure2.__init__(
            self, model, radius, offset, hse_up_key, hse_down_key
        )

    def _get_cb(self, residue):
        if residue.get_resname() == "GLY":
            return self._get_gly_cb_vector(residue)
        else:
            if residue.has_id("CB") and residue.has_id("CA"):
                c_b = residue["CB"].get_vector()
                c_a = residue["CA"].get_vector()
                return (c_b - c_a)

    def theta(self, r1, r2, r3):
        try:
            r1_v = r1["CA"].get_vector()
            r2_v = r2["CA"].get_vector()
            r3_v = r3["CA"].get_vector()
        except Exception:
            return None
        r1_r2_v = r2_v - r1_v
        r3_r2_v = r2_v - r3_v
        cos_theta = (r1_r2_v * r3_r2_v) / (r1_r2_v.norm() + r3_r2_v.norm())
        return acos(cos_theta)


# CA-CA
class ExposureCN(AbstractPropertyMap):
    """Residue exposure as number of CA atoms around its CA atom."""

    def __init__(self, model, CN_key, radius, offset):
        assert offset >= 0
        fs_map = {}
        fs_list = []
        fs_keys = []
        chains = model.get_chains()
        for chain1 in chains:
            chain1_seqlist = []
            for residue in chain1:
                if residue.get_id()[2] != ' ':
                    chain1_seqlist.append(f'{residue.get_id()[1]}{residue.get_id()[2]}')
                else:
                    chain1_seqlist.append(residue.get_id()[1])
            for i in chain1_seqlist:
                fs = 0
                position1 = ""
                add_code1 = ""
                for c in str(i):
                    if c.isdigit():
                        position1 += c
                    elif c.isalpha():
                        add_code1 += c
                if add_code1 == "":
                    add_code1 += " "
                r1 = chain1[" ", int(position1), add_code1]
                if not is_aa(r1) or not r1.has_id("CA"):
                    continue
                ca1 = r1["CA"]
                chains = model.get_chains()
                for chain2 in chains:
                    chain2_seqlist = []
                    for residue in chain2:
                        if residue.get_id()[2] != ' ':
                            chain2_seqlist.append(f'{residue.get_id()[1]}{residue.get_id()[2]}')
                        else:
                            chain2_seqlist.append(residue.get_id()[1])
                    for j in chain2_seqlist:
                        position2 = ""
                        add_code2 = ""
                        for c in str(j):
                            if c.isdigit():
                                position2 += c
                            elif c.isalpha():
                                add_code2 += c
                        if add_code2 == "":
                            add_code2 = " "
                        r2 = chain2[" ", int(position2), add_code2]
                        if not is_aa(r2) or not r2.has_id("CA"):
                            continue
                        ca2 = r2["CA"]
                        distance = ca1 - ca2
                        if chain1 is chain2 and distance <= offset:
                            continue
                        if distance <= radius:
                            fs += 1
                res_id = r1.get_id()
                chain_id = r1.get_parent().get_id()
                fs_map[(chain_id, res_id)] = fs
                fs_list.append((r1, fs))
                fs_keys.append((chain_id, res_id))
                r1.xtra[CN_key] = fs  # Add to xtra
        AbstractPropertyMap.__init__(self, fs_map, fs_keys, fs_list)


# CB-CB
class ExposureCB(AbstractPropertyMap):
    """Residue exposure as number of CB atoms around its CB atom."""

    def __init__(self, model, CB_key, radius, offset):

        assert offset >= 0
        fs_map = {}
        fs_list = []
        fs_keys = []
        special_cases = {}
        chains = model.get_chains()
        for chain1 in chains:
            chain1_seqlist = []
            for residue in chain1:
                if residue.get_id()[2] != ' ':
                    chain1_seqlist.append(f'{residue.get_id()[1]}{residue.get_id()[2]}')
                else:
                    chain1_seqlist.append(residue.get_id()[1])
            for i in chain1_seqlist:
                fs = 0
                position1 = ""
                add_code1 = ""
                for c in str(i):
                    if c.isdigit():
                        position1 += c
                    elif c.isalpha():
                        add_code1 += c
                if add_code1 == "":
                    add_code1 += " "
                r1 = chain1[" ", int(position1), add_code1]
                if not is_aa(r1):
                    continue
                if not r1.has_id("CB"):
                    ca1 = r1["CA"]
                    chains = model.get_chains()
                    for chain2 in chains:
                        chain2_seqlist = []
                        for residue in chain2:
                            if residue.get_id()[2] != ' ':
                                chain2_seqlist.append(f'{residue.get_id()[1]}{residue.get_id()[2]}')
                            else:
                                chain2_seqlist.append(residue.get_id()[1])
                        for j in chain2_seqlist:
                            position2 = ""
                            add_code2 = ""
                            for c in str(j):
                                if c.isdigit():
                                    position2 += c
                                elif c.isalpha():
                                    add_code2 += c
                            if add_code2 == "":
                                add_code2 += " "
                            r2 = chain2[" ", int(position2), add_code2]
                            if not is_aa(r2) or not r2.has_id("CB"):
                                continue
                            ca2 = r2["CB"]
                            distance = ca1 - ca2
                            if chain1 is chain2 and distance <= offset:
                                continue
                            if distance < radius:
                                fs += 1
                    res_id = r1.get_id()
                    chain_id = r1.get_parent().get_id()
                    fs_map[(chain_id, res_id)] = fs
                    fs_list.append((r1, fs))
                    fs_keys.append((chain_id, res_id))
                    r1.xtra[CB_key] = fs
                else:
                    ca1 = r1["CB"]
                    chains = model.get_chains()
                    for chain2 in chains:
                        chain2_seqlist = []
                        for residue in chain2:
                            if residue.get_id()[2] != ' ':
                                chain2_seqlist.append(f'{residue.get_id()[1]}{residue.get_id()[2]}')
                            else:
                                chain2_seqlist.append(residue.get_id()[1])
                        id_groups = {}
                        for j in chain2_seqlist:
                            position2 = ""
                            add_code2 = ""
                            for c in str(j):
                                if c.isdigit():
                                    position2 += c
                                elif c.isalpha():
                                    add_code2 += c
                            if add_code2 == "":
                                add_code2 += " "
                            r2 = chain2[" ", int(position2), add_code2]
                            if chain1 is chain2 and abs(int(position1) - int(position2)) <= offset:
                                continue
                            if not is_aa(r2) or not r2.has_id("CB"):
                                continue
                            ca2 = r2["CB"]
                            distance = ca1 - ca2
                            if distance < radius:
                                fs += 1
                    res_id = r1.get_id()
                    chain_id = r1.get_parent().get_id()
                    fs_map[(chain_id, res_id)] = fs
                    fs_list.append((r1, fs))
                    fs_keys.append((chain_id, res_id))
                    # Add to xtra
                    r1.xtra[CB_key] = fs
            AbstractPropertyMap.__init__(self, fs_map, fs_keys, fs_list)


def net_feature(pdb_file):
    molecule = bg.Pmolecule(pdb_file)
    network = molecule.network()
    # 全残基聚类系数
    Clustering_coefficient = nx.clustering(network)
    # 全残基介数中心性
    Betweenness_centrality = nx.betweenness_centrality(network)
    return Clustering_coefficient, Betweenness_centrality


res_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO',
            'SER', 'THR', 'TRP', 'TYR', 'VAL']

residue_three_to_one_letter = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}


def main():
    input_file_name = sys.argv[1]
    pdb_dir = sys.argv[2]
    
    '''  340 346
    if not os.path.exists('biopython_networkx.feat'):
        os.mkdir('biopython_networkx.feat')
    
    if not os.path.exists(f'biopython_networkx.feat/{input_file_name}'):
        os.mkdir(f'biopython_networkx.feat/{input_file_name}')
    else:
        shutil.rmtree(f'biopython_networkx.feat/{input_file_name}')
        os.mkdir(f'biopython_networkx.feat/{input_file_name}')'''

    print(f'\nhse2.py:Start to calculate {input_file_name}\'s network characteristics')
    with open(f'{input_file_name}', 'r') as file:
        for line in file:
            try:
                pdb_file = line.strip()
                file_path = f"{pdb_dir}/{pdb_file}.pdb"
#340 346
                if os.path.exists(file_path) and not os.path.exists(f"./biopython_networkx.feat/{input_file_name}/{pdb_file}.feat"):

                    # 计算网络属性
                    Clustering_coefficient, Betweenness_centrality = net_feature(file_path)

                    p = PDBParser()
                    structure = p.get_structure("protein", file_path)
                    model = structure[0]
                    total_res_num = 0
                    chain_res_seq = {}
                    for chain in model:
                        Chain_ID = chain.id
                        if Chain_ID not in chain_res_seq:
                            chain_res_seq[Chain_ID] = []
                        for residue in chain:
                            chain_res_seq[Chain_ID].append(residue.id[1])
                    for num_list in chain_res_seq.values():
                        total_res_num += len(num_list)

                    CN_calculator_4A = ExposureCN(model, CN_key='CA_CA_4A', radius=4.0, offset=0)
                    hse_calculator_4A = HSExposureCB(model, hse_up_key='HSE_UP_4A', hse_down_key='HSE_DOWN_4A',
                                                     radius=4.0, offset=0)

                    CN_calculator_5A = ExposureCN(model, CN_key='CA_CA_5A', radius=5.0, offset=0)
                    hse_calculator_5A = HSExposureCB(model, hse_up_key='HSE_UP_5A', hse_down_key='HSE_DOWN_5A',
                                                     radius=5.0, offset=0)

                    CN_calculator_6A = ExposureCN(model, CN_key='CA_CA_6A', radius=6.0, offset=0)
                    hse_calculator_6A = HSExposureCB(model, hse_up_key='HSE_UP_6A', hse_down_key='HSE_DOWN_6A',
                                                     radius=6.0, offset=0)

                    CN_calculator_7A = ExposureCN(model, CN_key='CA_CA_7A', radius=7.0, offset=0)
                    hse_calculator_7A = HSExposureCB(model, hse_up_key='HSE_UP_7A', hse_down_key='HSE_DOWN_7A',
                                                     radius=7.0, offset=0)

                    CN_calculator_8A = ExposureCN(model, CN_key='CA_CA_8A', radius=8.0, offset=0)
                    hse_calculator_8A = HSExposureCB(model, hse_up_key='HSE_UP_8A', hse_down_key='HSE_DOWN_8A',
                                                     radius=8.0, offset=0)

                    CN_calculator_9A = ExposureCN(model, CN_key='CA_CA_9A', radius=9.0, offset=0)
                    hse_calculator_9A = HSExposureCB(model, hse_up_key='HSE_UP_9A', hse_down_key='HSE_DOWN_9A',
                                                     radius=9.0, offset=0)

                    CN_calculator_10A = ExposureCN(model, CN_key='CA_CA_10A', radius=10.0, offset=0)
                    hse_calculator_10A = HSExposureCB(model, hse_up_key='HSE_UP_10A', hse_down_key='HSE_DOWN_10A',
                                                      radius=10.0, offset=0)

                    CN_calculator_11A = ExposureCN(model, CN_key='CA_CA_11A', radius=11.0, offset=0)
                    hse_calculator_11A = HSExposureCB(model, hse_up_key='HSE_UP_11A', hse_down_key='HSE_DOWN_11A',
                                                      radius=11.0, offset=0)

                    CN_calculator_CB_4A = ExposureCB(model, CB_key='CB_CB_4A', radius=4.0, offset=0)
                    CN_calculator_CB_5A = ExposureCB(model, CB_key='CB_CB_5A', radius=5.0, offset=0)
                    CN_calculator_CB_6A = ExposureCB(model, CB_key='CB_CB_6A', radius=6.0, offset=0)
                    CN_calculator_CB_7A = ExposureCB(model, CB_key='CB_CB_7A', radius=7.0, offset=0)
                    CN_calculator_CB_8A = ExposureCB(model, CB_key='CB_CB_8A', radius=8.0, offset=0)
                    CN_calculator_CB_9A = ExposureCB(model, CB_key='CB_CB_9A', radius=9.0, offset=0)
                    CN_calculator_CB_10A = ExposureCB(model, CB_key='CB_CB_10A', radius=10.0, offset=0)
                    CN_calculator_CB_11A = ExposureCB(model, CB_key='CB_CB_11A', radius=11.0, offset=0)

                    with open(f"./biopython_networkx.feat/{input_file_name}/{pdb_file}.feat", "a") as f:
                        f.write(
                            "Residue_contact_density_CACA_4A,hse_u_proportion_4A,hse_d_proportion_4A,"
                            "Residue_contact_density_CACA_5A,hse_u_proportion_5A,hse_d_proportion_5A,"
                            "Residue_contact_density_CACA_6A,hse_u_proportion_6A,hse_d_proportion_6A,"
                            "Residue_contact_density_CACA_7A,hse_u_proportion_7A,hse_d_proportion_7A,"
                            "Residue_contact_density_CACA_8A,hse_u_proportion_8A,hse_d_proportion_8A,"
                            "Residue_contact_density_CACA_9A,hse_u_proportion_9A,hse_d_proportion_9A,"
                            "Residue_contact_density_CACA_10A,hse_u_proportion_10A,hse_d_proportion_10A,"
                            "Residue_contact_density_CACA_11A,hse_u_proportion_11A,hse_d_proportion_11A,"
                            "Residue_contact_density_CBCB_4A,"
                            "Residue_contact_density_CBCB_5A,"
                            "Residue_contact_density_CBCB_6A,"
                            "Residue_contact_density_CBCB_7A,"
                            "Residue_contact_density_CBCB_8A,"
                            "Residue_contact_density_CBCB_9A,"
                            "Residue_contact_density_CBCB_10A,"
                            "Residue_contact_density_CBCB_11A,"
                            "betweenness_centrality,clustering_coefficient")
                        f.write("\n")
                        for chain in model:
                            Chain_ID = chain.id
                            for residue in chain:
                                for Remaining_residues in res_list:
                                    if residue.get_resname() != Remaining_residues:
                                        Residue_Chain_Number = f'{Chain_ID}{residue.id[1]}{residue.id[2]}'
                                        # 4A_CA-CA contact density
                                        Residue_contact_density_CACA_4A = residue.xtra['CA_CA_4A'] / total_res_num
                                        # Proportion of upper hemisphere exposure_4A
                                        if residue.xtra['CA_CA_4A'] != 0:
                                            RES_hse_u_proportion_CACA_4A = residue.xtra['HSE_UP_4A'] / residue.xtra[
                                                'CA_CA_4A']
                                            RES_hse_d_proportion_CACA_4A = residue.xtra['HSE_DOWN_4A'] / residue.xtra[
                                                'CA_CA_4A']
                                        else:
                                            RES_hse_u_proportion_CACA_4A = 0
                                            RES_hse_d_proportion_CACA_4A = 0

                                        # 5A_CA-CA contact density
                                        Residue_contact_density_CACA_5A = residue.xtra['CA_CA_5A'] / total_res_num
                                        # Proportion of upper hemisphere exposure_5A
                                        if residue.xtra['CA_CA_5A'] != 0:
                                            RES_hse_u_proportion_CACA_5A = residue.xtra['HSE_UP_5A'] / residue.xtra[
                                                'CA_CA_5A']
                                            RES_hse_d_proportion_CACA_5A = residue.xtra['HSE_DOWN_5A'] / residue.xtra[
                                                'CA_CA_5A']
                                        else:
                                            RES_hse_u_proportion_CACA_5A = 0
                                            RES_hse_d_proportion_CACA_5A = 0

                                        # 6A_CA-CA contact density
                                        Residue_contact_density_CACA_6A = residue.xtra['CA_CA_6A'] / total_res_num
                                        # Proportion of upper hemisphere exposure_4A
                                        if residue.xtra['CA_CA_6A'] != 0:
                                            RES_hse_u_proportion_CACA_6A = residue.xtra['HSE_UP_6A'] / residue.xtra[
                                                'CA_CA_6A']
                                            RES_hse_d_proportion_CACA_6A = residue.xtra['HSE_DOWN_6A'] / residue.xtra[
                                                'CA_CA_6A']
                                        else:
                                            RES_hse_u_proportion_CACA_6A = 0
                                            RES_hse_d_proportion_CACA_6A = 0

                                        # 7A_CA-CA contact density
                                        Residue_contact_density_CACA_7A = residue.xtra['CA_CA_7A'] / total_res_num
                                        # Proportion of upper hemisphere exposure_7A
                                        if residue.xtra['CA_CA_7A'] != 0:
                                            RES_hse_u_proportion_CACA_7A = residue.xtra['HSE_UP_7A'] / residue.xtra[
                                                'CA_CA_7A']
                                            RES_hse_d_proportion_CACA_7A = residue.xtra['HSE_DOWN_7A'] / residue.xtra[
                                                'CA_CA_7A']
                                        else:
                                            RES_hse_u_proportion = 0
                                            RES_hse_d_proportion = 0

                                        # 8A_CA-CA contact density
                                        Residue_contact_density_CACA_8A = residue.xtra['CA_CA_8A'] / total_res_num
                                        # Proportion of upper hemisphere exposure_8A
                                        if residue.xtra['CA_CA_8A'] != 0:
                                            RES_hse_u_proportion_CACA_8A = residue.xtra['HSE_UP_8A'] / residue.xtra[
                                                'CA_CA_8A']
                                            RES_hse_d_proportion_CACA_8A = residue.xtra['HSE_DOWN_8A'] / residue.xtra[
                                                'CA_CA_8A']
                                        else:
                                            RES_hse_u_proportion_CACA_8A = 0
                                            RES_hse_d_proportion_CACA_8A = 0

                                        # 9A_CA-CA contact density
                                        Residue_contact_density_CACA_9A = residue.xtra['CA_CA_9A'] / total_res_num
                                        # Proportion of upper hemisphere exposure_9A
                                        if residue.xtra['CA_CA_9A'] != 0:
                                            RES_hse_u_proportion_CACA_9A = residue.xtra['HSE_UP_9A'] / residue.xtra[
                                                'CA_CA_9A']
                                            RES_hse_d_proportion_CACA_9A = residue.xtra['HSE_DOWN_9A'] / residue.xtra[
                                                'CA_CA_9A']
                                        else:
                                            RES_hse_u_proportion_CACA_9A = 0
                                            RES_hse_d_proportion_CACA_9A = 0

                                        # 10A_CA-CA contact density
                                        Residue_contact_density_CACA_10A = residue.xtra['CA_CA_10A'] / total_res_num
                                        # Proportion of upper hemisphere exposure_10A
                                        if residue.xtra['CA_CA_10A'] != 0:
                                            RES_hse_u_proportion_CACA_10A = residue.xtra['HSE_UP_10A'] / residue.xtra[
                                                'CA_CA_10A']
                                            RES_hse_d_proportion_CACA_10A = residue.xtra['HSE_DOWN_10A'] / residue.xtra[
                                                'CA_CA_10A']
                                        else:
                                            RES_hse_u_proportion_CACA_10A = 0
                                            RES_hse_d_proportion_CACA_10A = 0

                                        # 11A_CA-CA contact density
                                        Residue_contact_density_CACA_11A = residue.xtra['CA_CA_11A'] / total_res_num
                                        # Proportion of upper hemisphere exposure_11A
                                        if residue.xtra['CA_CA_11A'] != 0:
                                            RES_hse_u_proportion_CACA_11A = residue.xtra['HSE_UP_11A'] / residue.xtra[
                                                'CA_CA_11A']
                                            RES_hse_d_proportion_CACA_11A = residue.xtra['HSE_DOWN_11A'] / residue.xtra[
                                                'CA_CA_11A']
                                        else:
                                            RES_hse_u_proportion_CACA_11A = 0
                                            RES_hse_d_proportion_CACA_11A = 0

                                        # CB-CB contact density
                                        Residue_contact_density_CBCB_4A = residue.xtra['CB_CB_4A'] / total_res_num
                                        Residue_contact_density_CBCB_5A = residue.xtra['CB_CB_5A'] / total_res_num
                                        Residue_contact_density_CBCB_6A = residue.xtra['CB_CB_6A'] / total_res_num
                                        Residue_contact_density_CBCB_7A = residue.xtra['CB_CB_7A'] / total_res_num
                                        Residue_contact_density_CBCB_8A = residue.xtra['CB_CB_8A'] / total_res_num
                                        Residue_contact_density_CBCB_9A = residue.xtra['CB_CB_9A'] / total_res_num
                                        Residue_contact_density_CBCB_10A = residue.xtra['CB_CB_10A'] / total_res_num
                                        Residue_contact_density_CBCB_11A = residue.xtra['CB_CB_11A'] / total_res_num

                                        f.write(
                                            f"{round(Residue_contact_density_CACA_4A, 5)},{round(RES_hse_u_proportion_CACA_4A, 5)},{round(RES_hse_d_proportion_CACA_4A, 5)},"
                                            f"{round(Residue_contact_density_CACA_5A, 5)},{round(RES_hse_u_proportion_CACA_5A, 5)},{round(RES_hse_d_proportion_CACA_5A, 5)},"
                                            f"{round(Residue_contact_density_CACA_6A, 5)},{round(RES_hse_u_proportion_CACA_6A, 5)},{round(RES_hse_d_proportion_CACA_6A, 5)},"
                                            f"{round(Residue_contact_density_CACA_7A, 5)},{round(RES_hse_u_proportion_CACA_7A, 5)},{round(RES_hse_d_proportion_CACA_7A, 5)},"
                                            f"{round(Residue_contact_density_CACA_8A, 5)},{round(RES_hse_u_proportion_CACA_8A, 5)},{round(RES_hse_d_proportion_CACA_8A, 5)},"
                                            f"{round(Residue_contact_density_CACA_9A, 5)},{round(RES_hse_u_proportion_CACA_9A, 5)},{round(RES_hse_d_proportion_CACA_9A, 5)},"
                                            f"{round(Residue_contact_density_CACA_10A, 5)},{round(RES_hse_u_proportion_CACA_10A, 5)},{round(RES_hse_d_proportion_CACA_10A, 5)},"
                                            f"{round(Residue_contact_density_CACA_11A, 5)},{round(RES_hse_u_proportion_CACA_11A, 5)},{round(RES_hse_d_proportion_CACA_11A, 5)},"
                                            f"{round(Residue_contact_density_CBCB_4A, 5)},{round(Residue_contact_density_CBCB_5A, 5)},{round(Residue_contact_density_CBCB_6A, 5)},{round(Residue_contact_density_CBCB_7A, 5)},{round(Residue_contact_density_CBCB_8A, 5)},{round(Residue_contact_density_CBCB_9A, 5)},{round(Residue_contact_density_CBCB_10A, 5)},{round(Residue_contact_density_CBCB_11A, 5)},"
                                            f"{round(Clustering_coefficient[Residue_Chain_Number], 5)},{round(Betweenness_centrality[Residue_Chain_Number], 5)}")
                                        f.write("\n")
                else:
                    with open(f'./biopython_networkx.feat/{input_file_name}/Not_exist_pdb.list', 'a') as out:
                        print(f'\nhse2.py:{pdb_dir}/{pdb_file}.pdb dosen\'t exist ')
                        out.write(f'{pdb_dir}/{pdb_file}.pdb dosen\'t exist\n ')
                    continue
            except Exception as e:
                with open(f'./biopython_networkx.feat/{input_file_name}/wrong_{input_file_name}.log',
                          'a') as wrong_file:
                    wrong_file.write(pdb_file)
                    wrong_file.write(str(traceback.format_exc()))
                    continue
        print(f'   \nhse2.py:{input_file_name}\'s network characteristics are calculated successfully!')


if __name__ == '__main__':
    main()

