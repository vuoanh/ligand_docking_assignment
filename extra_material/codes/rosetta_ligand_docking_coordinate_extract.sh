#!/bin/bash
# Check the number of the parameters
if [ "$#" -ne 2 ]
then
    echo "Wrong number of arguements. Usage: rosetta_ligand_docking_coordinate_extract.sh <list of Rosetta Docking models> <name of output>"
    exit 1
else
    file=$1
    if [ -s "$file" ]
    then 
	OUT=$2
	if [ -s "$OUT" ]
	then
	    rm $OUT
	fi
	for i in `cat $file`; do echo `grep HETATM ${i} |grep -v "           H"| awk '{print $7, $8, $9}'| tr '\n' ' '` >> $OUT ; done
	exit 0
    else
	echo "Error: the input list file does not exist"
	exit 1
    fi
fi

