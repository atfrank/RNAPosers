# RNAPosers
RNAPosers: Machine-Learning Pose Classifiers for RNA Containing Complexes.

This repo contains source code for RNAPosers' pose fingerprint and prediction (classification) modules. Pose fingerprint module is written in C++, and will generate executable `bin/featurize` once compiled. Prediction (classification) is done using Python script `src/rna_poser.py` and pre-trained classifiers with various parameter setting in `classifier/` . All the classifiers were trained with `Python 3.5.3` and `sklearn v0.19.2`, but should be compatible with both Python2 and Python 3 and any later version sklearn. The combined fingerprinting and classification process can be done by running `src/run.sh`.

**RNAPosers-plugin**: We also offer a PyMOL plugin version of RNAPosers (`rnaposerplugin.zip`) with graphical user interface as a supplement to the source code version. See [Using RNAPosers PyMOL plugin](#Using-RNAPosers-PyMOL-plugin) for installation instructions.

Manuscript (In submission): Chhabra, Sahil, Jingru Xie, and Aaron T. Frank. "RNAPosers: Machine Learning Classifiers For RNA-Ligand Poses." bioRxiv (2019): 702449.

## Prerequisite
* [Anaconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)

## Quick Start
### Install
```
git clone --depth=1 git@github.com:atfrank/RNAPosers.git
cd RNAPoser/
make clean
make
echo "export RNAPOSERS_PATH=$(pwd)" >> ~/.bashrc
source ~/.bashrc
```

### Setup environment
```
conda env create -f env/rnaposers.yml
conda activate rnaposers
```
<!---
# conda create --name rnaposers
# conda activate rnaposers
# conda install -c schrodinger pymol
# conda install -c schrodinger pymol-psico
# conda install -c openbabel openbabel
# conda install pandas
# conda install -c anaconda scikit-learn
-->

### Using RNAPosers

```
./src/run.sh -h
```

#### Arguments
- **receptor mol2**: .mol2 file of receptor structure. Default: tests/input/1AM0/receptor.mol2
- **ligand poses sd**: .sd file containing all ligand posees. Default: tests/input/1AM0/poses.sd
- **output file**: where to save scores. Default: tests/output/1AM0.txt
- **rmsd**: 1, 1.5, 2, 2.5. Default: 2
- **eta**: 2, 24, or 248 (2A, 2A and 4A, 2A 4A and 8A). Default: 248
- **stop frame**: only score the first several poses. Default: -1 (using all frames)

### Example
```
cd $RNAPOSERS_PATH
# Test 1: score the first 10 poses of 1AM0
./src/run.sh tests/input/1AM0/receptor.mol2 tests/input/1AM0/poses.sd tests/output/1AM0.txt 10 2.5 248
# Test 2: score the full set of poses of 2B57
./src/run.sh tests/input/2B57/receptor.mol2 tests/input/2B57/poses.sd tests/output/2B57.txt
```
#### Output
```
cat tests/output/1AM0.txt
# columns: prediction probability(0) probability(1)
1.000000 0.027000 0.973000
1.000000 0.023000 0.977000
1.000000 0.009000 0.991000
1.000000 0.062000 0.938000
1.000000 0.008000 0.992000
1.000000 0.050000 0.950000
1.000000 0.031000 0.969000
1.000000 0.010000 0.990000
1.000000 0.284000 0.716000
0.000000 0.687000 0.313000

cat tests/output/2B57.txt
1.000000 0.025000 0.975000
1.000000 0.015000 0.985000
1.000000 0.025000 0.975000
1.000000 0.005000 0.995000
1.000000 0.004000 0.996000
1.000000 0.004000 0.996000
1.000000 0.004000 0.996000
1.000000 0.005000 0.995000
1.000000 0.005000 0.995000
1.000000 0.004000 0.996000
...
```

### PyMOL plugin
The RNAPosers PyMOL plugin is compressed it as `rnaposerplugin.zip`. See PyMOL website: https://pymolwiki.org/index.php/Plugin_Manager for an instruction on installing PyMOL plugin from local file. Note that you still have to install this repo and setup environment to use the plugin.


### Additional Notes
- Ligand Atom types should to be one of the following SYBYL atom types:

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


## License
```
  Copyright University of Michigan.
  Author: Jingru Xie and Aaron T. Frank

```
