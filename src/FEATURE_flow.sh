#!/bin/sh
#input_list
input_list_name=$1
Chosen_model=$2
pdb_dir=$3
python2bin=$4
model=''
if [ "$Chosen_model" = "a" ]; then
    model='4.5A_ball'
elif [ "$Chosen_model" = "b" ]; then
    model='6A_ball'
elif [ "$Chosen_model" = "c" ]; then
    model='4.5A_to_6A_shell'
else
    echo -e "\n\nPlease select the correct model from [a,b,c]"
    exit 1
fi
echo -e "\n\nFEATURE_flow.sh:What you have chosen is:${Chosen_model}.${model}"
if [ -f ${input_list_name} ];then
#create FEATURE_DATA and export DATA_PATH
if [ ! -d ./FEATURE_data ]; then
    mkdir ./FEATURE_data
    mkdir ./FEATURE_data/pdb_FEATURE
    mkdir ./FEATURE_data/dssp_FEATURE
fi
#convert pdb name to FEATURE needed name
python Feature_PDBIDS.py $input_list_name $pdb_dir
export PDB_DIR=./FEATURE_data/pdb_FEATURE/$input_list_name
export DSSP_DIR=./FEATURE_data/dssp_FEATURE/$input_list_name
#prepare dssp file
while IFS= read -r line;
do
	pro_name=$(echo "$line" | tr -d ' ' )
	if [ -f ./FEATURE_data/pdb_FEATURE/${input_list_name}/${pro_name}.pdb ]; then
		mkdssp -i ./FEATURE_data/pdb_FEATURE/${input_list_name}/${pro_name}.pdb -o ./FEATURE_data/dssp_FEATURE/${input_list_name}/${pro_name}.dssp
	fi
done < ./rename_list/${input_list_name}/Feature_PDBIDS.list
if [ ! -d ./FEATURE.feat ];then
mkdir ./FEATURE.feat
fi
if [ ! -d ./FEATURE.feat/${input_list_name} ];then
mkdir ./FEATURE.feat/${input_list_name}
fi
if [ -d ./FEATURE.feat/${input_list_name}/$model ];then
rm -rf ./FEATURE.feat/${input_list_name}/$model
fi
#RUN FEATURE
mkdir ./FEATURE.feat/${input_list_name}/$model
bash FEATURE.sh ${input_list_name} ${model} ${python2bin}
python Feature.feat.py ${input_list_name} ${model} ${pdb_dir}
else
	echo -e "\n\n   FEATURE_flow.sh:${input_list_name}.list does not exist\n\n"
fi
