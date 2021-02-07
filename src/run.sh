#!/bin/bash

if [ "$1" == "-h" ]; then
  echo "usage: $0 [receptor mol2: .mol2 file of receptor structure. Default: tests/input/receptor_ligand/receptor.mol2]
                  [ligand poses sd: .sd file containing all ligand posees. Default: tests/input/receptor_ligand/poses.sd]
                  [output file: where to save scores. Default: tests/score.txt]
                  [rmsd: 1, 1.5, 2, 2.5. Default: 2]
                  [eta: 2, 24, or 248 (2A, 2A and 4A, 2A 4A and 8A). Default: 248]
                  [stop frame: only score the first several poses. Default: -1 (using all frames)]"

  echo "example: $0 tests/input/1AM0/receptor.mol2 tests/input/1AM0/poses.sd tests/output/1AM0.txt 10 2.5 248"
  exit 0
else
    receptor=${1:-tests/input/1AM0/receptor.mol2}
    ligand=${2:-tests/input/1AM0/poses.sd}
    score=${3:-tests/output/1AM0.txt}
    stop_frame=${4:-"-1"}
    rmsd=${5:-2}
    eta=${6:-248}
    python ${RNAPOSERS_PATH}/py/rnaposer.py $receptor $ligand $score $rmsd $eta $stop_frame
fi
