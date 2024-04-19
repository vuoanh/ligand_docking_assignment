#!/bin/bash
prefix=$1
native=`readlink -f $2`
#seed=1
mkdir -p ./analysis/${prefix}/
ROSETTA=~/apps/rosetta.source.release-334/

$ROSETTA/main/source/bin/rosetta_scripts.linuxgccrelease \
     -parser:script_vars native=${native} \
     -l ${prefix}_dock.lst \
     -score:analytic_etable_evaluation true \
     -out:file:scorefile ${prefix}_scores.out \
     -out:file:silent ${prefix}.silent \
     -out:file:silent_struct_type binary \
     -extra_res_fa ../ligand_prep/${prefix}/UNK.fa.params \
     -ex1 \
     -ex2 \
     -no_optH false \
     -flip_HNQ true \
     -ignore_ligand_chi true \
     -parser:protocol ./scripts/score_rmsd.xml \
     -restore_pre_talaris_2013_behavior true \
     -out:path:all ./analysis/${prefix}/ \
     #-constraints:cst_fa_file c1.cst 
