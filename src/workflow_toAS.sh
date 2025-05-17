#!/bin/sh
#input_list and input_pdb_dir
input_list_name=$1
pdb_dir=$2
if [ -f  ${input_list_name} ];then
echo ${input_list_name} > tmp_list_list
python   toAverStatistics.py tmp_list_list
bash Feature_toAS.sh ${input_list_name} a ${pdb_dir}
bash Feature_toAS.sh ${input_list_name} b ${pdb_dir}
bash Feature_toAS.sh ${input_list_name} c ${pdb_dir}
python merge_toAS.py ${input_list_name}  ./feat ./toAS.feat ${pdb_dir}
echo -e "\nPNFx:All the structure features(include toAS) from pdb are extracted successfully!!"
else
	echo "${input_list_name} does not exist"
fi