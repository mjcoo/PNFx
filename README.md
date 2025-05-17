# PNFx
Environment variables to be set
export PATH=$PATH:${dir}/PNFX/bin
export FEATURE_DIR=${dir}/PNFX/FEATURE_feature_list
export PYTHONPATH=$PYTHONPATH:${dir}/PNFX/lib/FEATURE/
export atomselect=${dir}/PNFX/bin/atomselector.py
export PATH=“/data/jwhuang/anaconda3/envs/PNFx/bin/python2.7:$PATH”
export python2bin=/data/jwhuang/anaconda3/envs/PNFx/bin/python2.7

Description of Files Stored in the Created Directories
feat: Contains network features, DSSP features, FEATURE features, and toAS features.

biopython_networkx.feat: Network features.
dssp_biopython.feat: DSSP features.
FEATURE.feat: FEATURE features.
toAS.feat: toAS features.

fix_pdb: Stores fixed PDB files.

FEAT: Contains network features, DSSP features, FEATURE features, and toAS features.

All_protein_toAS: toAS data calculated for all proteins in the list.
Workflow Instructions
Repair PDB files using process.py and save the fixed PDB files in ./pdb_fix/${input_list_name}.
Calculate hemisphere exposure and contact density using hse2.py, and store the features in ./feat/biopython_networkx.feat/${input_list_name}.
Calculate secondary structure features using dssp_biopython.py, and store the features in ./feat/dssp_biopython.feat/${input_list_name}.
Run FEATURE_WORKFLOW.sh to call Feature_flow.sh three times (once for each shell) to calculate FEATURE-class features, and store the features in ./feat/FEATURE.feat/$input_list_name under the corresponding model directory.
Execute the full workflow directly using workflow.sh.
