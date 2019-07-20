#!/bin/bash
if [[ $# -ne 11 ]]
then
    pdb=2b57
    echo "usage: $0 [base-directory] 
                          [pdbid] 
                          [ref-pdb-name] 
                          [ref-mol2-name] 
                          [ref-dcd-name] 
                          [last-frame] 
                          [rmsd-threshold: 1., 1.5, 2., or 2.5] 
                          [eta: 2, 24, or 248] 
                          [feature-file-prefix] 
                          [class-scores-prefix]
                          [selection]"
    
    echo "example: ${0} tests/input/${pdb}/ ${pdb} complex.pdb lig_${pdb}.mol2 complexes.dcd 2.5 248 output/${pdb}_features output/${pdb}_class_scores"
else
    source_dir=${1}
    pdbid=${2}
    pdb="${1}/${3}"
    mol2="${1}/${4}"
    dcd="${1}/${5}"
    stop=${6}
    rmsd=${7}
    eta=${8}
    feature_file=${9}
    class_score_file=${10}
    selection=${11}
    
    # check if rmsd valid
    if ! [[ ":1:1.5:2:2.5:" = *:$rmsd:* ]]
    then
        echo "ERROR: rmsd value not valid (need to be 1, 1.5, 2, 2.5). Exit."
        exit
    fi
    # check if eta valid
    if [ "$eta" == "2" ]
    then
        nEta=1
    elif [ "$eta" == "24" ]
    then
        nEta=2
    elif [ "$eta" == "248" ]
    then
        nEta=3
    else
        echo "ERROR: eta value not valid (need to be 2, 24 or 248). Exit."
        exit
    fi
    
    # featurize trajectory
    selatms=":ADE.C1' :ADE.C2 :ADE.C2' :ADE.C3' :ADE.C4 :ADE.C4' :ADE.C5 :ADE.C5' :ADE.C6 :ADE.C8 :ADE.N1 :ADE.N3 :ADE.N6 :ADE.N7 :ADE.N9 :ADE.O2' :ADE.O3' :ADE.O4' :ADE.O5' :ADE.OP1 :ADE.OP2 :ADE.P :CYT.C1' :CYT.C2 :CYT.C2' :CYT.C3' :CYT.C4 :CYT.C4' :CYT.C5 :CYT.C5' :CYT.C6 :CYT.N1 :CYT.N3 :CYT.N4 :CYT.O2 :CYT.O2' :CYT.O3' :CYT.O4' :CYT.O5' :CYT.OP1 :CYT.OP2 :CYT.P :GUA.C1' :GUA.C2 :GUA.C2' :GUA.C3' :GUA.C4 :GUA.C4' :GUA.C5 :GUA.C5' :GUA.C6 :GUA.C8 :GUA.N1 :GUA.N2 :GUA.N3 :GUA.N7 :GUA.N9 :GUA.O2' :GUA.O3' :GUA.O4' :GUA.O5' :GUA.O6 :GUA.OP1 :GUA.OP2 :GUA.P :URA.C1' :URA.C2 :URA.C2' :URA.C3' :URA.C4 :URA.C4' :URA.C5 :URA.C5' :URA.C6 :URA.N1 :URA.N3 :URA.O2 :URA.O2' :URA.O3' :URA.O4 :URA.O4' :URA.O5' :URA.OP1 :URA.OP2 :URA.P"
    rowatms=${selection}
    echo ${selection}
    ./bin/featurize \
    -etaStartPow 1 \
    -numEta $nEta \
    -cutoff 20 \
    -outfile ${feature_file} \
    -scalar 1 \
    -molecular 1 \
    -mol2 $mol2 \
    -normalization 1 \
    -rowatm "`echo ${rowatms}`" \
    -selatm "`echo ${selatms}`" \
    -stop $stop \
    -trj $dcd $pdb
    
    # get score
    python src/rna_poser.py --classifier classifier/eta${eta}/RF_${rmsd}.predictor.pkl  --features ${feature_file}_traj1.txt --output ${class_score_file}.txt

fi