#!/bin/bash
#SBATCH --job-name=rnaposer
#SBATCH -N 1
#SBATCH -n 1

# setup environment
module load gcc
module load anaconda

# intitialize 
base=/home/afrankz/local_software/repo/RNAPoser/tests/input/rna-protein/ten_3drpc/
rnaposer=/home/afrankz/local_software/repo/RNAPoser
cd ${rnaposer}

# loop over test cases
cutoffs="1 1.5 2 2.5"
etas="2 24 248"
pdbs="1B7F 1DFU 1JBS 1P6V 1WPU 1WSU 2ASB 2BH2 3BX2"
for pdb in $pdbs
do
    for cutoff in $cutoffs
    do
        for eta in $etas
        do
            ./src/rna_poser_validation.sh ${base}/ ${pdb} ${pdb}_complex.pdb ${pdb}_protein_from_pdb.mol2 ${pdb}_complex.dcd 100 2.5 248 ${base}/${pdb}_classifications_${cutoff}_${eta}.txt
        done
    done
done
