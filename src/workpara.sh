#!/bin/sh
input_list_name=$1
pdb_dir=$2
python PNFx.py -l ${input_list_name} -pdbdir ${pdb_dir}