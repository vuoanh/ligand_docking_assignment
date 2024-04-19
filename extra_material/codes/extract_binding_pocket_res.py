import numpy as np
import argparse

parser = argparse.ArgumentParser(description='filter the docked model that contact to at least a number of residues')
parser.add_argument('-i', '--input_pdb', required=True,help='input pdb file of the target-ligand complex to extract the'
                    + 'binding pocket residue.', dest='input_pdb')
parser.add_argument('-rc', '--target_chain', dest='target_chain', default='A', help='Chain name of target. Default: A')
parser.add_argument('-lc', '--ligand_chain', dest='ligand_chain', default='X', help='Chain name of ligand. Default: X')
parser.add_argument('-d', '--distance', default=6, type=float, dest='distance',
                    help='max distance (in angstrong) to be considered a contact. Default: 6')
parser.add_argument('-o', '--output', dest='output', default="output.txt",
                    help='output file with the list of IDs of the binding pocket residues. Default: output.txt')
args = parser.parse_args()

from Bio.PDB.PDBParser import PDBParser
parser = PDBParser()
three2one = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
    'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
    'THR':'T','TRP':'W','TYR':'Y','VAL':'V'
}

# Return the distance between two atoms
def atm_distance(atom1, atom2):
    sum = 1000
    if (atom1.element != 'H') & (atom2.element != 'H'):
        sum = 0
        for i in range(0,3):
            sum = sum + (atom1.get_coord()[i] -
                         atom2.get_coord()[i])**2
        sum = sum**0.5
    return sum

# return the min atom pairwise distance between two residues
def res_distance(res1, res2):
    min_dis = 1000
    for atom1 in res1:
        for atom2 in res2:
            min_dis = min(min_dis, atm_distance(atom1, atom2))
    return min_dis

# returns the number of contacted important residues of a model (start id = 1)
def contacted_res_num(model):
    pdb = parser.get_structure("pdb", model)[0]
    ligand = pdb[args.ligand_chain].get_list()[0]
    target = pdb[args.target_chain].get_list()
    f = open(args.output, 'w')
    print( "output the contacting residues into %s" % args.output)
    for resi in range(len(target)): 
        residue = target[resi]
        if res_distance(ligand, residue) <= args.distance:
            f.write("%d %s %s\n" % (resi + 1, args.target_chain, three2one[residue.get_resname()] ))
    f.close()
    print( "done")

def main():
    contacted_res_num(args.input_pdb)

if __name__ == "__main__":
    main()
