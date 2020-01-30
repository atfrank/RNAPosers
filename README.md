# RNAPosers
Machine-Learning Pose Classifiers for RNA Containing Complexes

## Quick Start
### Install
```
git clone git@github.com:atfrank/RNAPosers.git
cd RNAPoser/
make clean
make
echo "export RNAPOSERS_PATH=$(pwd)" >> ~/.bashrc
source ~/.bashrc
```

### Create python environment
```
conda install -c anaconda scikit-learn
TBC
```

### PyMOL plugin setup notes
See PyMOL website: https://pymolwiki.org/index.php/Plugin_Manager




## Usage Notes
- RNAPosers expect that the ordering of ligand atoms in the reference receptor-ligand pdb matches the order in the reference ligand mol2 file
- RNAPosers reads atomtypes from the reference ligand mol2 file
- Atomtypes should to be one of the following SYBYL atom types:

Description | Type
--- | ---
Carbon sp3 | C.3
Carbon sp2 | C.2
Carbon sp | C.1
Carbon aromatic | C.ar
Carbocation (guanadinium) | C.cat
Nitrogen sp3 | N.3
Nitrogen sp2 | N.2
Nitrogen sp | N.1
Nitrogen aromatic | N.ar
Nitrogen amide | N.am
Nitrogen trigonal planar | N.pl3
Nitrogen sp3 positively charged | N.4
Oxygen sp3 | O.3
Oxygen sp2 | O.2
Oxygen in carboxylates and phosphates | O.co2
Sulphur sp3 | S.3
Sulphur sp2 | S.2
Sulphoxide sulphur | S.o
Sulphone sulphur | S.o2
Phosphorus sp3 | P.3


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
$ cd /path/to/RNAPoser/
$ pdb=2b57
$ mkdir tests/output/
$ ./src/rna_poser.sh tests/input/${pdb}/ ${pdb} complex.pdb lig_${pdb}.mol2 complexes.dcd 2.5 2 tests/output/${pdb}_features tests/output/${pdb}_class_scores

cat tests/output/class_scores_${pdb}.txt
# columns: prediction probability(0) probability(1)

1.000000 0.004000 0.996000
1.000000 0.001000 0.999000
1.000000 0.004000 0.996000
1.000000 0.002000 0.998000
1.000000 0.001000 0.999000
1.000000 0.001000 0.999000
1.000000 0.000000 1.000000
1.000000 0.001000 0.999000
1.000000 0.001000 0.999000
1.000000 0.000000 1.000000
1.000000 0.000000 1.000000
...      ...      ...     
...      ...      ...     
...      ...      ...     
0.000000 0.981000 0.019000
0.000000 0.973000 0.027000
0.000000 0.971000 0.029000
0.000000 0.974000 0.026000
0.000000 0.968000 0.032000
0.000000 0.980000 0.020000
0.000000 0.963000 0.037000
0.000000 0.967000 0.033000
0.000000 0.975000 0.025000
0.000000 0.983000 0.017000
0.000000 0.984000 0.016000
```
## License
```
  Copyright University of Michigan.
  Author: Jingru Xie and Aaron T. Frank

```
