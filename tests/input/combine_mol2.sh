#This combine all mol2 for Ligand Atom Types Frequency


rnas="2o3w 2xnw 3fu2 3mum 3q50 3sd3 3slm 4erj 4fe5 4jf2 4lx5 4nya 4xwf 4yb0 5c7w 2b57 1f1t 2ydh 3npn 4aob 4kqy 4l81 5kpy"

for rna in ${rnas}
do

cat ${rna}/lig_${rna}.mol2 >> l.mol2

done
