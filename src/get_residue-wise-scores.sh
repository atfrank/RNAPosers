#!/bin/bash
if [[ $# -ne 1 ]]
then
    pdb=3BO2
    echo "usage: $0 <pdbid>"
    echo " e.g.: $0 ${pdb}"
else
    # initialize variable
    module load r
    pdb=$1
    base=/home/afrankz/local_software/repo/3dRPC-2.0/data3dRPC/rnaposer/${pdb}/bound
    rnaposer=/home/afrankz/local_software/repo/RNAPosers
    
    # get array of residues
    resids=`grep UNK ${base}/rnaposer_complex.pdb | awk '{print $6}' | uniq`

    # check if SLURM_ARRAY_TASK_ID greater than length of residue array
    rm -rf ${base}/classifications_all.txt
    for resid in $resids
    do
        if [[  -f ${base}/classifications_test_${resid}.txt ]]
        then
            cat -n ${base}/classifications_test_${resid}.txt | awk -v pdb=${pdb} -v resid=${resid} '{print pdb, resid, $0}' >> ${base}/classifications_all.txt
        fi
    done
    
    # get composite scores
    ${rnaposer}/src/./get_composite_class_scores.R --data_transform="quantile" -o ${base}/${pdb}_scores.txt ${base}/classifications_all.txt
    
fi
