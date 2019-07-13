
rna="2g5k"

/home/itssahil/PROJECTS/Mol_Feturizer/Test_systems_manuscript/RNAPoser-master/bin/rna_poser ../${rna}/complex.pdb -trj ../${rna}/complexes.dcd -mol2 ../${rna}/lig_${rna}.mol2 -mode L -rdock ../${rna}/Scores.txt

#awk '{if ( $3 <= 2.5) printf("%d %s\n", 1, $3); else printf("%d %s\n", 0, $3) }' /home/itssahil/PROJECTS/Mol_Feturizer/Test_systems_manuscript/decoys_for_sahil/${rna}/rmsd.txt > /home/itssahil/PROJECTS/Mol_Feturizer/Test_systems_manuscript/decoys_for_sahil/${rna}/tags.txt

#paste dump.txt /home/itssahil/PROJECTS/Mol_Feturizer/Test_systems_manuscript/decoys_for_sahil/${rna}/tags.txt > prediction_R_${rna}.txt

