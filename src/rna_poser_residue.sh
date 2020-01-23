#!/bin/bash
#SBATCH --job-name=rnaposer
#SBATCH -N 1
#SBATCH -n 1

# setup environment
module load gcc
module load anaconda

# intitialize 
pdb=3BO2
base=/home/afrankz/local_software/repo/3dRPC-2.0/data3dRPC/rnaposer/${pdb}/bound
rnaposer=/home/afrankz/local_software/repo/RNAPosers

# goto working directory
cd ${rnaposer}

# get array of residues
resids=(`grep UNK ${base}/rnaposer_complex.pdb | awk '{print $6}' | uniq`)

# check if SLURM_ARRAY_TASK_ID greater than length of residue array
if [[ ${SLURM_ARRAY_TASK_ID} > ${#resids[@]} ]]
then
    exit
else
    resid=${resids[${SLURM_ARRAY_TASK_ID}]}
    echo $resid
    ./src/rna_poser.sh ${base}/ ${pdb} rnaposer_complex.pdb rnaposer_protein_from_pdb.mol2 rnaposer_complex.dcd 100 2.5 248 ${base}/features_${resid} ${base}/classifications_test_${resid} `echo ":UNK/${resid}."`
fi


# sbatch --time=00:20:00 --array=0-95%95 -N1 -p frank /home/afrankz/local_software/repo/RNAPosers/src/rna_poser_residue.sh
