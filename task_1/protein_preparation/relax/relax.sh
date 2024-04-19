#!/bin/bash
START=`readlink -e $1`
XML=`readlink -e $2`
PREFIX=$3
JOB_NUM=$4

cpu_num=`grep "^cpu\\scores" /proc/cpuinfo | uniq |  awk '{print $4}'`

echo $cpu_num
ROSETTA=~/apps/rosetta.source.release-334/

seq 1 $JOB_NUM | xargs -n1 -I% -P $((cpu_num-1)) $ROSETTA/main/source/bin/rosetta_scripts.linuxgccrelease \
    -parser:protocol $XML\
    -out:pdb_gz \
    -out:prefix ${PREFIX}_%_ \
    -nstruct 1 \
    -out:file:scorefile ${PREFIX}_%_scores.out \
    -s $START \
    -native $START \
    -relax:constrain_relax_to_start_coords \
    -relax:constrain_relax_to_native_coords 
