# RNAPoser
Machine-Learning Pose Classifier for RNA-Ligand Complexes

## Install
```shell
$ cd /path/to/RNAPoser/
$ make clean
$ make
```

## Usage manual
```shell
$ ./src/rna_poser.sh [input directory] [id] [receptor] [poses] [rmsd] [eta]
Options: [id: identifier, the output file will be saved in working_dir/${id}/]
         [input directory: the path to folder that contain the input pdb and sd file]
         [receptor: receptor coor file in pdb format, should be in the input directory]
         [poses: poses coordinates in sd format containing multiple frames, should be in the input directory]
         [rmsd: predictors trained with different rmsd (possible values: 1, 1.5, 2, 2.5)]
         [eta: eta values to use for featurization (possible values: 2, 24, 248)]

```
## Example
```shell
$ pdb=2b57
$ ./src/rna_poser.sh tests/input/${pdb}/ ${pdb} receptor.pdb poses.sd

# columns: prediction probability(0) probability(1)
file: working_dir/2b57/prediction.txt
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
```
To validate paper results:
```shell
$ pdb=2b57
$ ./src/rna_poser_validation.sh tests/input/${pdb}/ ${pdb} complex.pdb lig_${pdb}.mol2 complexes.dcd

```

## License
```
  Copyright University of Michigan.
  Author: Jingru Xie and Aaron T. Frank

```
