import sys
import os
import re
import csv
import pandas as pd
import json
import multiprocessing
from multiprocessing import Pool

def flow_para(genelist):
    cmd = f"bash workpara.sh {genelist}.list /data/jwhuang/data >>  /data/jwhuang/log/{genelist}FEATURE.log 2>&1"
    print(cmd)
    os.system(cmd)

def main():
    prot_list = ["Uniprot" + str(i) for i in range(1, 340)]
    prot_list.append("Uniprot341")
    prot_list.append("Uniprot342")
    prot_list.append("Uniprot343")
    prot_list.append("Uniprot344")
    prot_list.append("Uniprot345")
    prot_list2 = ["Uniprot" + str(i) for i in range(347, 378)]
    prot_list.extend(prot_list2)
    pool = Pool(60)
    pool.map(flow_para, prot_list)



if __name__ == '__main__':
    main()
