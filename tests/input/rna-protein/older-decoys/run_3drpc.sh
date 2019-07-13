#!/bin/bash
source ~/.bashrc 

# go home
cd /home/afrankz/local_software/repo/RNAPoser/tests/input/rna-protein

# generate 3dRPC parameter file
cat << EOF >| parameter_3drpc
RPDock.receptor = protein_original_3dRPC.pdb
RPDock.receptor.chain = A
RPDock.ligand = rna_original_3dRPC.pdb
RPDock.ligand.chain = B
RPDock.outfile = docking.out
RPDock.grid_step = 0.5
RPDock.out_pdb = 100
EOF


# loop over RNAs
rnas="1B7F 1DFU 1JBS 1P6V 1WPU 1WSU 2ASB 2BH2 2QUX 3BX2"
rnas="1B7F"
for rna in $rnas
do
    # set up and copy need files
    rm -rf ${rna}
    mkdir -p ${rna}
    cd ${rna}
    cp ../${rna}_protein_original_3dRPC.pdb protein_original_3dRPC.pdb
    cp ../${rna}_rna_original_3dRPC.pdb rna_original_3dRPC.pdb
    cp ../parameter_3drpc .

    # execute 3dRPC
    ${HOME_3dRPC}/3dRPC -mode 9 -system 9 -par parameter_3drpc
    cd ../
done

