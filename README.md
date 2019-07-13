# RNAPosers
Machine-Learning Pose Classifiers for RNA-Ligand Complexes

## Install RNPosers
```shell
$ cd /path/to/RNAPoser/
$ make clean
$ make
```

## Install Scikit Learn
```shell
$ conda install -c anaconda scikit-learn 
```

## Usage manual
```shell
$ ./src/rna_poser.sh $0 [base-directory]  [pdbid] [ref-pdb-name] [ref-mol2-name] [ref-dcd-name] [stop] [rmsd-threshold: 1., 1.5, 2., or 2.5] [eta: 2, 24, or 248] [feature-file prefix][class-score file prefix]"
Arguments: [base-directory: the path to folder that contain the receptor-ligand pdb file, ligand mol2 file, receptor-ligand poses dcd file]
           [pdbid: identifier, the output file will be saved in working_dir/${id}/]
           [ref-pdb-name: name of reference RNA-ligand pdb file]
           [ref-mol2-name: name of reference ligand mol2 file]
           [ref-dcd-name: name of receptor-ligand pose dcd file]
           [rmsd-threshold: definition of nativeness. Either: 1., 1.5, 2., or 2.5]
           [eta: widths of gaussian used for generating pose fingerprint. Either 2={2 Å}[used in the manuscript], 24={2 and 4 Å}, or 248 = {2, 4, and 8 Å}]
           [feature-file prefix: prefix file to which pose features should outputted to]
           [class-score file prefix: prefix file to which classification scores should outputted to]
```
## Example
```shell
$ pdb=2b57
$ ./src/rna_poser.sh tests/input/${pdb}/ ${pdb} complex.pdb lig_${pdb}.mol2 complexes.dcd 2.5 248 output/${pdb}_features output/${pdb}_class_scores

# columns: prediction probability(0) probability(1)
cat tests/output/class_scores_${pdb}.txt
1.000000 0.014000 0.986000
1.000000 0.006000 0.994000
1.000000 0.007000 0.993000
1.000000 0.005000 0.995000
1.000000 0.006000 0.994000
1.000000 0.006000 0.994000
1.000000 0.005000 0.995000
1.000000 0.002000 0.998000
1.000000 0.002000 0.998000
1.000000 0.003000 0.997000
1.000000 0.003000 0.997000
1.000000 0.003000 0.997000
1.000000 0.002000 0.998000
1.000000 0.003000 0.997000
  ...

## License
```
  Copyright University of Michigan.
  Author: Jingru Xie and Aaron T. Frank

```
