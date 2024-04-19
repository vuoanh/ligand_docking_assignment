#!/bin/bash
###############################################################################################
# Oanh Vu
# This script clusters the rosettaCM models that have total scores smaller than a cutoff value
################################################################################################


if [ "$#" -eq 3 ]
then
    r_e_cutoff=$1
    cluster_num=$2
    fasta=$3
    output_dir=./
    cd ${output_dir}
    tail -n+3  ./*.out| grep S|awk -v n=${r_e_cutoff} '$2 < n {print $2,$NF".pdb.gz"}'| sort -g > top_models.sc
    awk '{print $2}' top_models.sc > top_models.list
    echo "before removing corrupted models, there are `wc -l top_models.list` models"
    cat top_models.list |xargs -n1 -I@ -P 10 gunzip @
    sed -i 's/.gz//g' top_models.*
    echo removing incomplete pdb file before the analysis
    for i in `cat top_models.list` ; do echo ${i} `cat ${i} | sed '/^[[:space:]]*$/d'|tail -n1` | grep 'END';done | awk '{print $1}' > top_models_temp.list
    mv top_models_temp.list top_models.list
    echo "after removing corrupted models, there are `wc -l top_models.list` models"
    cat top_models.list|xargs -n1 -I@ grep ' @' top_models.sc > top_models_temp.sc
    mv top_models_temp.sc top_models.sc
    echo superimpose all models to the best scored msdodel
    python ~vuot2/bin/superimpose_rmsd.py -c A -o rmsd.out -s `tail -1 ${fasta}` -u `tail -1 ${fasta}` -l top_models.list -p 10
    echo clustering the model with clusco_cpu
    clusco_cpu -l top_models.list -s rmsd -o distance.matrix 0 ${cluster_num}
    paste rmsd.out <(tail -n+5  top_models.list.clustering1|awk '{print $1, $2}' ) top_models.sc | awk '{print $1, $2, $5, $4}'| column -t| sort -gk3 > rmsd.score.cluster
    # output graph with cluster labels
    ~vuot2/bin/plotting/multiple_scatter_plot.py -i rmsd.score.cluster -o rmsd.score.cluster.png -x 1 -X RMSD -y 2 -Y 'Total Score' -g 3
    echo "The final graph is rmsd.score.cluster.png"
    cat top_models.list |xargs -n1 -I@ -P 10 gzip @
else
    echo "This script clusters the rosettaCM models that have total scores smaller than a cutoff value, then did a nice visualization of the clustering results"
    echo "Wrong command. Usage post_rosettaCM_clustering.sh <Rosetta score cutoff> <number of cluster> <fasta_file>"
fi
