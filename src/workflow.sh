#!/bin/sh
#input_data
input_list_name=$1
pdb_dir=$2
python2bin=$3
if [ -f  ${input_list_name} ];then
#python>3,Bio,dssp
#proccess pdb
#python process.py $input_list_name $pdb_dir
#biopython+networkx
#python hse2.py $input_list_name ./fix_pdb/$input_list_name
#biopython+dssp
#python dssp_biopython.py $input_list_name ./fix_pdb/$input_list_name
#Feature
bash FEATURE_flow.sh ${input_list_name} a ./fix_pdb/$input_list_name ${python2bin}
bash FEATURE_flow.sh ${input_list_name} b ./fix_pdb/$input_list_name ${python2bin}
bash FEATURE_flow.sh ${input_list_name} c ./fix_pdb/$input_list_name ${python2bin}
#merge
python merge.py $input_list_name ./biopython_networkx.feat ./dssp_biopython.feat ./FEATURE.feat  $pdb_dir
else
	echo "${input_list_name} does not exist"
fi