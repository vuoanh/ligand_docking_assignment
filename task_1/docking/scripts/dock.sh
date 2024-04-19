#!/bin/bash
prefix=$1
#seed=1
mkdir -p ./output_files_${prefix}/
ROSETTA=~/apps/rosetta.source.release-334/

seq 1 8| xargs -n1 -I% -P4 $ROSETTA/main/source/bin/rosetta_scripts.linuxgccrelease \
     -parser:script_vars x=24.274 y=2.696 z=-63.992 maxddg=10\
     -out:prefix ${prefix}_% \
     -nstruct 625 \
     -s templates/%_${prefix}.pdb \
     -score:analytic_etable_evaluation true \
     -out:file:scorefile ${prefix}_%_scores.out \
     -out:pdb_gz \
     -extra_res_fa ../ligand_prep/${prefix}/UNK.fa.params \
     -ex1 \
     -ex2 \
     -no_optH false \
     -flip_HNQ true \
     -ignore_ligand_chi true \
     -constraints:cst_fa_weight 50 \
     -parser:protocol ./scripts/dock.xml \
     -restore_pre_talaris_2013_behavior true \
     -out:path:all ./output_files_${prefix}/
     #-constraints:cst_fa_file c1.cst 
