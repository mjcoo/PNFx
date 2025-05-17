#!/bin/sh

#input_list
input_list_name=$1
Chosen_model=$2
#echo "Please select an output Feature_model"
#echo "a: ball with a radius of 4.5A    b:ball with a radius of 6A   c:4.5A to 6A shells"
#read Chosen_model

model=''

if [ "$Chosen_model" = "a" ]; then
    model='4.5A_ball'
elif [ "$Chosen_model" = "b" ]; then
    model='6A_ball'
elif [ "$Chosen_model" = "c" ]; then
    model='4.5A_to_6A_shell'
else
    echo "Please select the correct model from [a,b,c]"
    exit 1
fi
echo -e "\n\nWhat you have chosen is:${Chosen_model}.${model}"

if [ -f ${input_list_name} ];then
if [ ! -d toAS.feat ];then
mkdir toAS.feat
fi
if [ ! -d toAS.feat/${input_list_name} ];then
mkdir toAS.feat/${input_list_name}
fi
if [ -d toAS.feat/${input_list_name}/${model} ];then
rm -rf toAS.feat/${input_list_name}/${model}
fi
if [ ! -d FEAT ];then
mkdir FEAT
fi
#RUN Feature_toAS
mkdir toAS.feat/${input_list_name}/${model}
python Feature.feat_toAS.py ${input_list_name} ${model}

else
	echo -e "\n\n   Feature_toAS.sh:${input_list_name}.list does not exist"
fi
