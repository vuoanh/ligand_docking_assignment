# Structure preparation
Since the compound #13 and #22 is closest to the compound #18, I decided to use the structure of the compound #18 with ERK3 as the starting structure for docking.
Download the structure PDB-ID 6YLL to the current directory and then use pymol to superimpose four crystal structures (chain A - D) to chain A, and save those structure into 4 PDB files (6yll_A.pdb, 6yll_B.pdb, 6yll_C.pdb, and 6yll_D.pdb)
Extract only the atom coordinate lines from the pdb files
```
python2.7 clean_pdb.py 6yll_A.pdb A
grep "OWQ A 401" 6yky_A.pdb > 6yky_lig_A.pdb 
sed -i 's/OWQ A 401/OWQ X 401/g' 6yky_lig_A.pdb
cat 6yky_A_A.pdb 6yky_lig_A.pdb > 6yky_1.pdb
```
repeat the same process with structure B, C, D to create 6yky_2.pdb, 6yky_3.pdb, and 6yky_4.pdb

## Adding missing loops
Unising Rosetta remodling application, I added the missing loops into 4 crystal structures
## Optimiaze the complex structures with Rosetta
Remove residues at the beginning and the end of each structure to make sure that all 4 structures have the same sequence ERK3(17-319)
Run the repacking and minimization with Rosetta Relax application
```
for i in `seq 1 4`; do bash relax.sh 6yky_${i}.pdb relax_cart_to_starting_coord_no_selection.xml 6yky_${i} 10;done
```
For starting structure, generate 10 relaxed models, and select top 2 models as the template for docking.

## Ligand preparation
### Compound #18
Extract the coordinate of compound 18 from the crystal structure
```
grep OWQ 6yky_1_0001.pdb > 6yky_lig_1.pdb
obabel -i pdb 6yky_lig_1.pdb -o sdf -O lig18_1.sdf
```

Generate a conformer library with 300 conformations for compound #18 and align them at the central triazolo[4,5-d]pyrimidin-5-amine ring
```
~/apps/bcl/build/linux64_release/bin/bcl-apps-static.exe molecule:ConformerGenerator -conformation_comparer SymmetryRMSD 0.0 -max_iterations 8000 -top_models 300 -cluster -ensemble_filenames lig18_1.sdf -conformers_single_file lig18_conformer.sdf -explicit_aromaticity -add_h -neutralize
~/apps/bcl/build/linux64_release/bin/bcl-apps-static.exe molecule:AlignToScaffold lig18_1.sdf lig18_conformer.sdf lig18_conformer_aligned.sdf -align_ensemble_atoms 1 17 18 -explicit_aromaticity -align_scaffold_atoms 1 17 18
```
Generate Rosetta params file for the ligand for RosettaLigand dockinh
```
~/apps/rosetta.source.release-334/main/source/scripts/python/public/molfile_to_params.py -p UNK -n UNK --clobber --extra_torsion_output --centroid --conformers-in-one-file lig18_conformer_aligned.sdf

~/apps/rosetta.source.release-334/main/source/scripts/python/public/molfile_to_params.py -p UNK -n UNK --clobber --extra_torsion_output --centroid lig18_1.sdf

echo "PDB_ROTAMERS UNK.fa_conformers.pdb" >> UNK.fa.params
```
### Compound #22
Use ChemDraw to make the sdf file for compound #22
Use obabel to generate 3D structure file SDF
```
obabel -i sdf c13.sdf -o sdf -O c13-3d.sdf -h --gen3d
```
Generate a conformer library with 300 conformations for compound #18 and align them at the central triazolo[4,5-d]pyrimidin-5-amine ring
```
~/apps/bcl/build/linux64_release/bin/bcl-apps-static.exe molecule:ConformerGenerator -conformation_comparer SymmetryRMSD 0.0 -max_iterations 8000 -top_models 300 -cluster -ensemble_filenames c22-3d.sdf -conformers_single_file c22_conformer.sdf -explicit_aromaticity -add_h -neutralize

~/apps/bcl/build/linux64_release/bin/bcl-apps-static.exe molecule:AlignToScaffold c22-3d.sdf c22_conformer.sdf c22_conformer_aligned.sdf -align_ensemble_atoms 2 6 8 -explicit_aromaticity -align_scaffold_atoms 2 6 8
```

Generate Rosetta params file for the ligand for RosettaLigand docking
```
id=c22; ~/apps/rosetta.source.release-334/main/source/scripts/python/public/molfile_to_params.py -p UNK -n UNK --clobber --extra_torsion_output --centroid --conformers-in-one-file ${id}_conformer_aligned.sdf

~/apps/rosetta.source.release-334/main/source/scripts/python/public/molfile_to_params.py -p UNK -n UNK --clobber --extra_torsion_output --centroid ${id}-3d.sdf

echo "PDB_ROTAMERS UNK.fa_conformers.pdb" >> UNK.fa.params
```
### Compound #13

Use ChemDraw to make the sdf file for compound #22
Use obabel to generate 3D structure file SDF
```
obabel -i sdf c13.sdf -o sdf -O c13-3d.sdf -h --gen3d
```
Generate a conformer library with 300 conformations for compound #18 and align them at the central triazolopyrimidin-5-amine ring
```
~/apps/bcl/build/linux64_release/bin/bcl-apps-static.exe molecule:ConformerGenerator -conformation_comparer SymmetryRMSD 0.0 -max_iterations 8000 -top_models 300 -cluster -ensemble_filenames c13-3d.sdf -conformers_single_file c13_conformer.sdf -explicit_aromaticity -add_h -neutralize

~/apps/bcl/build/linux64_release/bin/bcl-apps-static.exe molecule:AlignToScaffold c13-3d.sdf c13_conformer.sdf c13_conformer_aligned.sdf -align_ensemble_atoms 1 3 8 -explicit_aromaticity -align_scaffold_atoms 1 3 8
```
Generate Rosetta params file for the ligand for RosettaLigand docking
```
id=c22; ~/apps/rosetta.source.release-334/main/source/scripts/python/public/molfile_to_params.py -p UNK -n UNK --clobber --extra_torsion_output --centroid --conformers-in-one-file ${id}_conformer_aligned.sdf

~/apps/rosetta.source.release-334/main/source/scripts/python/public/molfile_to_params.py -p UNK -n UNK --clobber --extra_torsion_output --centroid ${id}-3d.sdf

echo "PDB_ROTAMERS UNK.fa_conformers.pdb" >> UNK.fa.params
```

