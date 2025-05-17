#!/bin/sh
list_name=$1
model=$2
python2bin=$3
while IFS= read -r line;
do
	pro_name=$(echo "$line" | tr -d ' ' )
if [ "$model" = "6A_ball" ]; then
    if [ -f ./FEATURE_data/pdb_FEATURE/${list_name}/${pro_name}.pdb ]; then
        $python2bin $atomselect -a ca $pro_name > ./FEATURE.feat/${list_name}/${model}/${pro_name}_ca_${model}.ptf
        featurize -P ./FEATURE.feat/${list_name}/${model}/${pro_name}_ca_${model}.ptf -n 1 -w 6 ${pro_name} > ./FEATURE.feat/${list_name}/${model}/${pro_name}_ca_${model}.ff
    fi
elif [ "$model" = "4.5A_ball" ]; then
    if [ -f ./FEATURE_data/pdb_FEATURE/${list_name}/${pro_name}.pdb ]; then
        $python2bin $atomselect -a ca $pro_name > ./FEATURE.feat/${list_name}/${model}/${pro_name}_ca_${model}.ptf
        featurize -P ./FEATURE.feat/${list_name}/${model}/${pro_name}_ca_${model}.ptf -n 1 -w 4.5 ${pro_name} > ./FEATURE.feat/${list_name}/${model}/${pro_name}_ca_${model}.ff
    fi
elif [ "$model" = "4.5A_to_6A_shell" ]; then
    if [ -f ./FEATURE_data/pdb_FEATURE/${list_name}/${pro_name}.pdb ]; then
        $python2bin $atomselect -a ca $pro_name > ./FEATURE.feat/${list_name}/${model}/${pro_name}_ca_${model}.ptf
        featurize -P ./FEATURE.feat/${list_name}/${model}/${pro_name}_ca_${model}.ptf -n 4 -w 1.5 ${pro_name} > ./FEATURE.feat/${list_name}/${model}/${pro_name}_ca_${model}.ff
    fi
fi
done < ./rename_list/${list_name}/Feature_PDBIDS.list
line_number=1
while IFS= read -r line;
do
	pro_name=$(echo "$line" | tr -d ' ' )
	pro_origin=$(sed -n "${line_number}p" ${list_name})
	pro_origin_name=$(echo "$pro_origin" | tr -d ' ')
	if [ -f ./FEATURE.feat/${list_name}/${model}/${pro_name}_ca_${model}.ptf ] && [ -f ./FEATURE.feat/${list_name}/${model}/${pro_name}_ca_${model}.ff  ]; then
		mv ./FEATURE.feat/${list_name}/${model}/${pro_name}_ca_${model}.ptf ./FEATURE.feat/${list_name}/${model}/${pro_origin_name}_ca_${model}.ptf
		mv ./FEATURE.feat/${list_name}/${model}/${pro_name}_ca_${model}.ff ./FEATURE.feat/${list_name}/${model}/${pro_origin_name}_ca_${model}.ff
	fi
	((line_number++))
done <  ./rename_list/${list_name}/Feature_PDBIDS.list
for file in ./FEATURE.feat/${list_name}/${model}/*; do
    new_file=$(echo "$file" | sed "s/^'//;s/'$//;s/\\r//g")
    if [ "$file" != "$new_file" ]; then
        mv -i -- "$file" "$new_file"
    fi
done
