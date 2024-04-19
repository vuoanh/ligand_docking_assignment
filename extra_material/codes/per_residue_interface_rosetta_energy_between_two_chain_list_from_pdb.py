#!/home/ubuntu/anaconda3/envs/pyrosetta/bin/python

import pyrosetta
import numpy as np
from pyrosetta.teaching import *
from argparse import ArgumentParser
import os

# create the command line parser
parser = ArgumentParser(
    """This script calculate per-residue between chains 
    interface scores (or rosetta predicted binding free energy) broken down accross input PDB files. 
    The PDB file name, the order numbers of two lists of chains that the user need to investigate 
    the interface rosetta score, and the name of the file containing the b-factor values need to be given.
    The interface between 2 same chains will not be calculated"""
)

parser.add_argument(
    "-p",
    "--pdb_list",
    dest="pdb_file_list",
    required=True,
    help="""The list of PDB filenames. this file won\"t be changed by this script."""
)
parser.add_argument(
    "-b",
    "--b_factors",
    dest="b_factors_file",
    required=True,
    help="""An output file containing the per-residue interface score of the input protein complex ensemble. 
    The file format is expected as <chain_id> <sequence_id> <interaction e mean> <interaction e std>"""
)
parser.add_argument(
    "-f",
    "--fist_chain_list",
    dest="first_chain_list",
    nargs='+',
    type = int,
    required=True,
    help="Chain order (a natural number) of the first chain list (first chain start from 1). e.g. -f 1 4"
)
parser.add_argument(
    "-s",
    "--second_chains",
    nargs='+', type=int,
    dest="second_chain_list",
    required=True,
    help="Chain orders (natural number) of the second chain list . e.g. -s 2 3"
)
parser.add_argument(
    "-a",
    "--print_all_pdb",
    type=bool,
    dest="print_all_pdb",
    default="False",
    help="set to True if the user want to print all values for indivisual input pdb files"
)
# parse the command line
args = parser.parse_args()

# start pyrosetta
pyrosetta.init()

## Helper functions
# Calculate the per-residue interface score and stores the values as an array
def calculate_ddg(pose, f_chains, s_chains, ddg_arr):
    # the default scoring function is ref2015
    sfxn = get_score_function(True)
    # compute and show the components of the score
    sfxn.show(pose)
    # compute the pair-wise interaction energies and then divide it evenly to 2 residues
    for f_order in f_chains:
        for s_order in s_chains:
            if f_order != s_order:
                f_start = pose.chain_begin( f_order )
                f_last = pose.chain_end( f_order )
                s_start = pose.chain_begin( s_order )
                s_last = pose.chain_end( s_order )
                for res1 in range( f_start, f_last + 1 ):
                    for res2 in range( s_start, s_last + 1 ):
                        edge = pose.energies().energy_graph().find_energy_edge( res1, res2 )
                        if edge:
                            ddg = edge.dot( sfxn.weights() )
                            ddg_arr[ res1 - 1 ] += ddg / 2
                            ddg_arr[ res2 - 1 ] += ddg / 2
    return ddg_arr

# calculate and report mean of std of per-residue interface score
def print_ddg_bfactor(pose, ddg_arr, bfactor_filename):
    f = open( bfactor_filename, 'w' )
    f.write("chain_id res_id mean std\n")
    for res in range( 1, len( ddg_arr ) + 1 ):
        f.write( "%s %s %.2f %.2f\n" % ( pose.pdb_info().chain(res), 
                                        pose.pdb_info().number(res), 
                                        np.mean( ddg_arr [ res - 1 ] ), 
                                        np.std( ddg_arr [ res - 1 ] ) 
                                       ) )
    f.close()

# report mean of std of per-residue interface score of every complexes in the input protein ensemble
def print_full_ddg_b_factor(pose, pdb_list, ddg_arr, bfactor_filename):
    f = open( "full_version_" + bfactor_filename, 'w' )
    f.write("chain_id res_id ")
    for pdb_name in pdb_list:
        f.write( "%s " % pdb_name)
    f.write("\n")
    for res in range( 1, len( ddg_arr ) + 1 ):
        f.write( "%s %s " % ( pose.pdb_info().chain(res), pose.pdb_info().number( res ) ) )
        for ddg_value in ddg_arr[ res - 1 ]:
            f.write("%.2f " % ddg_value)
        f.write("\n")
    f.close()

# remove pdb file names of the non-existent pdb files from the input pdb list
def clean_pdb_list(input_pdb_list, pdb_list):
    for pdb in input_pdb_list:
        if os.path.exists(pdb):
            pdb_list.append( pdb )
        else:
            print(" The input file %s does not exist. Remove this filename from the input pdb list" % pdb)
    return pdb_list
    
## Computes the interface scores and then writes out the output files
bfactor_filename = args.b_factors_file
input_pdb_list = np.loadtxt( args.pdb_file_list, dtype='S200' )
pdb_list = []

clean_pdb_list(input_pdb_list, pdb_list)
struct_num = len( pdb_list )

num_res = len( pyrosetta.pose_from_pdb( pdb_list[0] ) )
ddgs = np.zeros( ( struct_num, num_res ) )
struct_id = 0
for pdb in pdb_list:
    pdb_pose = pyrosetta.pose_from_pdb( pdb )
    ## Generate the scoring functions for the strutucture
    calculate_ddg(pdb_pose, args.first_chain_list, args.second_chain_list, ddgs[struct_id] )
    struct_id += 1
trans_ddgs = np.transpose(ddgs)
assert len(trans_ddgs) == num_res
print_ddg_bfactor(pdb_pose, trans_ddgs, bfactor_filename)
if args.print_all_pdb:
    print_full_ddg_b_factor(pdb_pose, pdb_list,trans_ddgs, bfactor_filename)
