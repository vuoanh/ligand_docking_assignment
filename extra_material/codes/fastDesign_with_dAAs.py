#!/home/oanhvu/anaconda3/bin/python python
# import packages
import logging
logging.basicConfig(level=logging.INFO)
import pyrosetta
import pyrosetta.toolbox
import os
import argparse
import sys
parser = argparse.ArgumentParser(description='Author: Oanh Vu. '
                                 'This script runs PyRosetta to design an input protein using all D and L canonical amino acids')
parser.add_argument("-n","--n_struct", default=1, type=int, help="number of output models generated. Default is 1")
parser.add_argument("-i","--input", required=True, help="name of the amber initial coordinate file")
parser.add_argument("-p","--prefix", required=True, help="prefix for names of the output files")
parser.add_argument("-r","--resfile", required=True, help="name of resfile for design")
parser.add_argument("-d","--out_dir", default=os.path.abspath(os.curdir), help="name of the output directory, default is the current directory")
parser.add_argument("-g","--interface_score", default=False, type=bool, help="output interface score, default: False")
parser.add_argument("-j","--interface_jump", required='--interface_score' in sys.argv, type=int, default=1, 
                    help="jump number for calculating the interface score, only required if the interface_score option is set to True. Default: 1")
parser.add_argument("-f","--output_format", default="silent", choices=['silent', 'pdb'], help="output format: pdb or silent")
parser.add_argument("-b","--bb_fix", nargs='+',type=int, help="indices of residues that need to fix bb atoms")
args = parser.parse_args()
n_struct = args.n_struct
prefix = args.prefix
input_pdb = os.path.abspath(args.input)
resfile=os.path.abspath(args.resfile)
os.makedirs(args.out_dir, exist_ok=True)
output_dir = os.path.abspath(args.out_dir)

#### HELPER FUNCTION
def write_file(filename, string):
    f = open(filename, "w")
    f.write(string)
    f.close()
    
#### RUNNING DESIGN
pyrosetta.init("-ignore_unrecognized_res 1 -ex1 -ex2aro -detect_disulf 0")
start_pose = pyrosetta.pose_from_pdb(input_pdb)
# read the list of modified residues from resfile
f = open(resfile, "r")
lines = f.readlines()
design_residues = []
for line in lines:
    line = line.rstrip()
    if line:
        res = line.split(' ')[0]
        if res.isdigit():
            design_residues.append(int(res))
f.close()
print("design_residues are: ", design_residues)
pose = start_pose.clone()
scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")

# Build packer palettes with L and D amino acids
palette = pyrosetta.rosetta.core.pack.palette.CustomBaseTypePackerPalette()
palette.parse_additional_residue_types("DALA,DASP,DGLU,DPHE,DHIS,DILE,DLYS,DLEU,DMET,DASN,DPRO,DGLN,DARG,DSER,DTHR,DVAL,DTRP,DTYR")
# task operation on DAA design
posPhi = pyrosetta.rosetta.core.select.residue_selector.PhiSelector()
posPhi.set_select_positive_phi(True)
negPhi = pyrosetta.rosetta.core.select.residue_selector.PhiSelector()
negPhi.set_select_positive_phi(False)
# <ReadResfile name="l_res" filename="./l_res.txt" selector="negPhi_pep"/>
d_res_design = """
PIKAA X[DALA]X[DCYS]X[DASP]X[DGLU]X[DPHE]X[DHIS]X[DILE]X[DLYS]X[DLEU]X[DMET]X[DASN]X[DPRO]X[DGLN]X[DARG]X[DSER]X[DTHR]X[DVAL]X[DTRP]X[DTYR]
start
"""
write_file("%s/d_res.txt" % output_dir, d_res_design)
l_res_design = """
PIKAA ACDEFHIKLMNPQRSTVWY
start
"""
write_file("%s/l_res.txt" % output_dir, l_res_design)
daa_design_restrictions = pyrosetta.rosetta.core.pack.task.operation.ReadResfile("%s/d_res.txt" % output_dir)
daa_design_restrictions.set_residue_selector(posPhi)
laa_design_restrictions = pyrosetta.rosetta.core.pack.task.operation.ReadResfile("%s/l_res.txt" % output_dir)
laa_design_restrictions.set_residue_selector(negPhi)

# The task factory accepts all the task operations
tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
tf.set_packer_palette(palette) ## Only needed for noncanonical design
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())

# Include the resfile
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.ReadResfile(resfile))
tf.push_back(daa_design_restrictions)
tf.push_back(laa_design_restrictions)

# Convert the task factory into a PackerTask to take a look at it
packer_task = tf.create_task_and_apply_taskoperations(pose)

# prevent movement of the residues whose bb atoms needs to be fixed
cyclization_fix = pyrosetta.rosetta.core.select.movemap.move_map_action(0)
terminal_res = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
for res in args.bb_fix :
    terminal_res.append_index(res)
    
# View the PackerTask
print("details of the PackerTask:")
print(packer_task)

# Set up a MoveMapFactory
mmf = pyrosetta.rosetta.core.select.movemap.MoveMapFactory()
mmf.all_bb(setting=True)
mmf.all_bondangles(setting=True)
mmf.all_bondlengths(setting=True)
mmf.all_chi(setting=True)
mmf.all_jumps(setting=True)

# fix terminal residue bb atoms
mmf.add_bb_action(cyclization_fix, terminal_res)
mmf.set_cartesian(setting=True)
display_pose = pyrosetta.rosetta.protocols.fold_from_loops.movers.DisplayPoseLabelsMover()
display_pose.tasks(tf)
display_pose.movemap_factory(mmf)
display_pose.apply(pose)
fr = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn_in=scorefxn, standard_repeats=1)
fr.cartesian(True)
fr.set_task_factory(tf)
fr.set_movemap_factory(mmf)
fr.min_type("lbfgs_armijo_nonmonotone") # For non-Cartesian scorefunctions, use "dfpmin_armijo_nonmonotone"

# Write output silent file and score file
out_put_score_file = "%s/%s_per_res.sc" % (output_dir, prefix)
sc_f = open(out_put_score_file, "w", buffering=1)
if args.interface_score:
    sc_f.write("model,rmsd,delta total score,interface score")
else:
    sc_f.write("model,rmsd,delta total score")
for res in design_residues:
    sc_f.write(",res %d,name,score" % res) 
sc_f.write("\n")

# write the checkpoint file
checkpointfile = os.path.abspath("%s/%s.ckp" % (output_dir, prefix))
checkpoint = 0
if os.path.isfile(checkpointfile) and os.path.getsize(checkpointfile) > 0 :
    with open(checkpointfile) as f:
        first_line = f.readline().strip('\n')
        assert first_line.isdigit() and int(first_line) < n_struct, "the checkpoint number must be a integer number less than %d" % n_struct
        checkpoint = int(first_line)
        
### BEGIN SOLUTION
for i in range(checkpoint, n_struct):  
    # apply fast design instruction on the pose
    fr.apply(pose)
    # calculate CA RMSD of the resulting pose to the starting pose
    rmsd = pyrosetta.rosetta.core.scoring.CA_rmsd(start_pose, pose)
    # calculate interace score
    delta_total_score = scorefxn(pose) - scorefxn(start_pose)
    if args.interface_score:
        interface_analyzer = pyrosetta.rosetta.protocols.analysis.InterfaceAnalyzerMover()
        interface_analyzer.set_interface_jump(args.interface_jump)
        interface_analyzer.apply(pose)
        sc_f.write("%s_%.4d,%.2f,%.2f,%.2f" % (prefix, i, rmsd, delta_total_score, interface_analyzer.get_interface_dG()))
    else:
        sc_f.write("%s_%.4d,%.2f,%.2f" % (prefix, i, rmsd, delta_total_score))
    pose.pdb_info().name("%s_%.4d" % (prefix, i))
    # output to silent structure file
    if args.output_format == 'silent':
        pyrosetta.io.poses_to_silent(pose, "%s/%s.silent" % (output_dir, prefix))
    else :
        pose.dump_pdb("%s/%s_%.4d.pdb" % (output_dir,prefix, i))
        print("here")
    # calculate per-residue energy improvement for each of redesigned residue
    for j in design_residues:
        pose_total_score = pyrosetta.rosetta.protocols.relax.get_per_residue_scores(pose, pyrosetta.rosetta.core.scoring.ScoreType.total_score)[j]
        start_pose_total_score = pyrosetta.rosetta.protocols.relax.get_per_residue_scores(start_pose, pyrosetta.rosetta.core.scoring.ScoreType.total_score)[j]
        sc_f.write(",%s,%s,%.2f" % ( start_pose.residue(j).name(), pose.residue(j).name(), pose_total_score - start_pose_total_score))
    sc_f.write("\n")
    write_file(checkpointfile, "%d" % (i+1))
sc_f.close()
