#!/bin/bash

if [[ $# -ne 5 ]]
then
    echo "usage: $0 [ref-pdb-name]
                    [ref-mol2-name]
                    [ref-dcd-name]
                    [feature-file]
                    [stop]"

    echo "example: $0 test/complex.pdb test/lig.mol2 test/complexes.dcd 2.5 248 test/feature test/score.txt 10"
else
    pdb=$1
    mol2=$2
    dcd="$3"
    feature_file=$4
    stop_frame=$5

    # featurize trajectory
    echo "[RNAPosers Debugging] Running Featurization... This may take several minutes."
    selatms=":ADE.C1' :ADE.C2 :ADE.C2' :ADE.C3' :ADE.C4 :ADE.C4' :ADE.C5 :ADE.C5' :ADE.C6 :ADE.C8 :ADE.N1 :ADE.N3 :ADE.N6 :ADE.N7 :ADE.N9 :ADE.O2' :ADE.O3' :ADE.O4' :ADE.O5' :ADE.OP1 :ADE.OP2 :ADE.P :CYT.C1' :CYT.C2 :CYT.C2' :CYT.C3' :CYT.C4 :CYT.C4' :CYT.C5 :CYT.C5' :CYT.C6 :CYT.N1 :CYT.N3 :CYT.N4 :CYT.O2 :CYT.O2' :CYT.O3' :CYT.O4' :CYT.O5' :CYT.OP1 :CYT.OP2 :CYT.P :GUA.C1' :GUA.C2 :GUA.C2' :GUA.C3' :GUA.C4 :GUA.C4' :GUA.C5 :GUA.C5' :GUA.C6 :GUA.C8 :GUA.N1 :GUA.N2 :GUA.N3 :GUA.N7 :GUA.N9 :GUA.O2' :GUA.O3' :GUA.O4' :GUA.O5' :GUA.O6 :GUA.OP1 :GUA.OP2 :GUA.P :URA.C1' :URA.C2 :URA.C2' :URA.C3' :URA.C4 :URA.C4' :URA.C5 :URA.C5' :URA.C6 :URA.N1 :URA.N3 :URA.O2 :URA.O2' :URA.O3' :URA.O4 :URA.O4' :URA.O5' :URA.OP1 :URA.OP2 :URA.P"
    rowatms=":UNK."
    if [ "$stop_frame" -gt 0 ]
    then
        $RNAPOSERS_PATH/bin/featurize \
        -etaStartPow 1 \
        -etaBase 2 \
        -numEta 1 \
        -cutoff 20 \
        -outfile ${feature_file} \
        -scalar 1 \
        -molecular 1 \
        -mol2 $mol2 \
        -normalization 0 \
        -rowatm "`echo ${rowatms}`" \
        -selatm "`echo ${selatms}`" \
        -stop ${stop_frame} \
        -trj $dcd $pdb  > /dev/null
    else
        $RNAPOSERS_PATH/bin/featurize \
        -etaStartPow 1 \
        -etaBase 2 \
        -numEta 1 \
        -cutoff 20 \
        -outfile ${feature_file} \
        -scalar 1 \
        -molecular 1 \
        -mol2 $mol2 \
        -normalization 0 \
        -rowatm "`echo ${rowatms}`" \
        -selatm "`echo ${selatms}`" \
        -trj $dcd $pdb  > /dev/null
    fi
    
    mv ${feature_file}_traj1.txt ${feature_file}.txt
fi
