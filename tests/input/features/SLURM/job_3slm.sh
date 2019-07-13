#!/bin/bash
#SBATCH --job-name=grp3slm
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 03:00:00
#SBATCH -p frank

module load gcc/4.9.4
module load babel


# get ligand info
#cat ../3slm/complex.pdb | sed 's/O1P/OP1/g' | sed 's/O2P/OP2/g' | sed 's/C.2/C2 /g' | sed 's/N.1/N1 /g' | sed 's/N.2/N2 /g' | sed 's/N.3/N3 /g' | sed 's/N.4/N4 /g' | sed 's/O.2/O2 /g' | sed 's/P.3/P  /g' | sed 's/OP3/OP3 /g' > ..//3slm/complex_fix.pdb
#mv ..//3slm/complex_fix.pdb ..//3slm/complex.pdb
#babel -ipdb ..//3slm/lig_3slm.pdb -omol2 ..//3slm/lig_3slm.mol2
#grep 'UNK' ..//3slm/lig_3slm.mol2| awk -v rna=3slm '{print rna $0}' > ml_data/info_3slm.txt
#cp ..//3slm/rmsd.txt ml_data/rmsd_3slm.txt
# featurize trajectory
selatms=":ADE.C1' :ADE.C2 :ADE.C2' :ADE.C3' :ADE.C4 :ADE.C4' :ADE.C5 :ADE.C5' :ADE.C6 :ADE.C8 :ADE.N1 :ADE.N3 :ADE.N6 :ADE.N7 :ADE.N9 :ADE.O2' :ADE.O3' :ADE.O4' :ADE.O5' :ADE.OP1 :ADE.OP2 :ADE.P :CYT.C1' :CYT.C2 :CYT.C2' :CYT.C3' :CYT.C4 :CYT.C4' :CYT.C5 :CYT.C5' :CYT.C6 :CYT.N1 :CYT.N3 :CYT.N4 :CYT.O2 :CYT.O2' :CYT.O3' :CYT.O4' :CYT.O5' :CYT.OP1 :CYT.OP2 :CYT.P :GUA.C1' :GUA.C2 :GUA.C2' :GUA.C3' :GUA.C4 :GUA.C4' :GUA.C5 :GUA.C5' :GUA.C6 :GUA.C8 :GUA.N1 :GUA.N2 :GUA.N3 :GUA.N7 :GUA.N9 :GUA.O2' :GUA.O3' :GUA.O4' :GUA.O5' :GUA.O6 :GUA.OP1 :GUA.OP2 :GUA.P :URA.C1' :URA.C2 :URA.C2' :URA.C3' :URA.C4 :URA.C4' :URA.C5 :URA.C5' :URA.C6 :URA.N1 :URA.N3 :URA.O2 :URA.O2' :URA.O3' :URA.O4 :URA.O4' :URA.O5' :URA.OP1 :URA.OP2 :URA.P"

cat ..//3slm/complex.pdb | awk '/UNK/{printf "%s \n", $3}' | grep -v "H" | awk '{printf ":%s.%s ", "UNK", $1}' > lig_pattern_3slm.txt
sed 's/ *$//' lig_pattern_3slm.txt > row_pattern_3slm.txt
rowatms="`cat row_pattern_3slm.txt`"

#awk '{if ( $3 <= 2.5) printf("%d %s\n", 1, $3); else printf("%d %s\n", 0, $3) }' /home/itssahil/PROJECTS/Mol_Feturizer/Test_systems_manuscript/decoys_for_sahil/3slm/rmsd.txt > /home/itssahil/PROJECTS/Mol_Feturizer/Test_systems_manuscript/decoys_for_sahil/3slm/tags.txt

awk '{if ( $3 <= 1) printf("%d %s\n", 1, $3); else printf("%d %s\n", 0, $3) }' /home/itssahil/PROJECTS/Mol_Feturizer/Test_systems_manuscript/decoys_for_sahil/3slm/rmsd.txt > /home/itssahil/PROJECTS/Mol_Feturizer/Test_systems_manuscript/decoys_for_sahil/3slm/tags_1.txt

#sh /home/itssahil/PROJECTS/Mol_Feturizer/Test_systems_manuscript/RNAPoser-master/make_predictor.sh R
#/home/itssahil/PROJECTS/Mol_Feturizer/Test_systems_manuscript/RNAPoser-master/bin/rna_poser ..//3slm/complex.pdb -trj ..//3slm/complexes.dcd -mol2 ..//3slm/lig_3slm.mol2 -mode R -rdock ..//3slm/Scores.txt
#sed -i -e "1d" prediction.txt
#paste prediction.txt /home/itssahil/PROJECTS/Mol_Feturizer/Test_systems_manuscript/decoys_for_sahil/3slm/tags.txt > prediction_R_3slm.txt

#sh /home/itssahil/PROJECTS/Mol_Feturizer/Test_systems_manuscript/RNAPoser-master/make_predictor.sh L
#/home/itssahil/PROJECTS/Mol_Feturizer/Test_systems_manuscript/RNAPoser-master/bin/rna_poser ..//3slm/complex.pdb -trj ..//3slm/complexes.dcd -mol2 ..//3slm/lig_3slm.mol2 -mode L -rdock ..//3slm/Scores.txt
#mv prediction.txt prediction_3slm.txt
#sed -i -e "1d" prediction_3slm.txt
#paste prediction_3slm.txt /home/itssahil/PROJECTS/Mol_Feturizer/Test_systems_manuscript/decoys_for_sahil/3slm/tags.txt > prediction_L_3slm.txt
#rm prediction_3slm.txt

#sh /home/itssahil/PROJECTS/Mol_Feturizer/Test_systems_manuscript/RNAPoser-master/make_predictor.sh RL
#/home/itssahil/PROJECTS/Mol_Feturizer/Test_systems_manuscript/RNAPoser-master/bin/rna_poser ..//3slm/complex.pdb -trj ..//3slm/complexes.dcd -mol2 ..//3slm/lig_3slm.mol2 -mode RL -rdock ..//3slm/Scores.txt
#mv prediction.txt prediction_3slm.txt
#sed -i -e "1d" prediction_3slm.txt
#paste prediction_3slm.txt /home/itssahil/PROJECTS/Mol_Feturizer/Test_systems_manuscript/decoys_for_sahil/3slm/tags.txt > prediction_RL_3slm.txt
#rm prediction_3slm.txt

/home/itssahil/PROJECTS/Mol_Feturizer/New_Systems/AtomicFeaturizer/bin/featurize -molecular 1 -normalization 1 -etaStartPow 1 -numEta 3 -cutoff 20.0 -outfile 3slm_  -mol2 ..//3slm/lig_3slm.mol2 -rowatm "`echo ${rowatms}`" -selatm "`echo ${selatms}`" ..//3slm/complex.pdb -trj ..//3slm/complexes.dcd
#/home/itssahil/PROJECTS/Mol_Feturizer/New_Systems/AtomicFeaturizer/bin/featurize -molecular 1 -normalization 1 -etaStartPow 1 -numEta 3 -cutoff 20.0 -outfile 3slm_  -mol2 ..//3slm/lig_3slm.mol2 -rowatm "" -selatm "" ..//3slm/complex.pdb -trj ..//3slm/complexes.dcd

#awk '{if ( $3 <= 2.5) printf("%d %s\n", 1, $3); else printf("%d %s\n", 0, $3) }' /home/itssahil/PROJECTS/Mol_Feturizer/New_Systems/decoys_for_sahil/featurize/ml_data/rmsd_3slm.txt > /home/itssahil/PROJECTS/Mol_Feturizer/New_Systems/decoys_for_sahil/featurize/ml_data/tags_3slm.txt
paste 3slm_traj1.txt /home/itssahil/PROJECTS/Mol_Feturizer/New_Systems/decoys_for_sahil/3slm/Scores.txt /home/itssahil/PROJECTS/Mol_Feturizer/New_Systems/decoys_for_sahil/3slm/tags_1.txt > features_3slm_1.txt
# clean up 

