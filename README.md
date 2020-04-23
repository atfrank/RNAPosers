# RNAPosers
RNAPosers: Machine-Learning Pose Classifiers for RNA Containing Complexes.

This repo contains source code for RNAPosers' pose fingerprint and prediction (classification) modules. Pose fingerprint module is written in C++, and will generate executable `bin/featurize` once compiled. Prediction (classification) is done using Python script `src/rna_poser.py` and pre-trained classifiers with various parameter setting in `classifier/` . All the classifiers were trained with `Python 3.5.3` and `sklearn v0.19.2`, but should be compatible with both Python2 and Python 3 and any later version sklearn. The combined fingerprinting and classification process can be done by running `rna_posers.sh`.

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
- **receptor mol2**: .mol2 file of receptor structure. Default: tests/input/receptor_ligand/receptor.mol2
- **ligand poses sd**: .sd file containing all ligand posees. Default: tests/input/receptor_ligand/poses.sd
- **output file**: where to save scores. Default: tests/score.txt
- **rmsd**: 1, 1.5, 2, 2.5. Default: 2
- **eta**: 2, 24, or 248 (2A, 2A and 4A, 2A 4A and 8A). Default: 248
- **stop frame**: only score the first several poses. Default: -1 (using all frames)

### Example
```
cd $RNAPOSERS_PATH
./src/run.sh tests/input/receptor_ligand/receptor.mol2 tests/input/receptor_ligand/poses.sd tests/score.txt 10 2.5 248
```
#### Output
```
cat tests/score.txt
# columns: prediction probability(0) probability(1)
0.000000 0.697000 0.303000
0.000000 0.701000 0.299000
0.000000 0.678000 0.322000
0.000000 0.678000 0.322000
0.000000 0.701000 0.299000
0.000000 0.678000 0.322000
0.000000 0.694000 0.306000
0.000000 0.684000 0.316000
0.000000 0.682000 0.318000
0.000000 0.745000 0.255000
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
