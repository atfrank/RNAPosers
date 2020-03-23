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
git clone --depth=1 https://github.com/atfrank/RNAPosers.git
cd RNAPosers/
make clean
make
echo "export RNAPOSERS_PATH=$(pwd)" >> ~/.bashrc
source ~/.bashrc
```

### Setup environment
```
conda create --name rnaposers
conda activate rnaposers
conda install -c anaconda scikit-learn
```

### Using RNAPosers

```
./src/rna_poser.sh -h
usage: ./rna_poser.sh [base-directory]
                      [pdbid]
                      [ref-pdb-name]
                      [ref-mol2-name]
                      [ref-dcd-name]
                      [rmsd-threshold: 1., 1.5, 2., or 2.5]
                      [eta: 2, 24, or 248]
                      [feature-file-prefix]
                      [class-scores-prefix]
example: ./rna_poser.sh tests/input/2b57/ 2b57 complex.pdb lig_2b57.mol2 complexes.dcd 2.5 248 output/2b57_features output/2b57_class_scores
```

#### Arguments
- **base-directory**: the path to folder that contain the receptor-ligand pdb file, ligand mol2 file, receptor-ligand poses dcd file
- **pdbid**: identifier of current structure
- **ref-pdb-name**: name of reference RNA-ligand pdb file]
- **ref-mol2-name**: name of reference ligand mol2 file]
- **ref-dcd-name**: name of receptor-ligand pose dcd file]
- **rmsd-threshold**: definition of nativeness. Either: `1`, `1.5`, `2`, or `2.5`
- **eta**: a Guassian width parameter for pose fingerprint. `eta=2` means {2 Å} (used in the manuscript), `eta=24` means {2Å and 4 Å}, and `eta=248` means {2, 4, and 8 Å}. The higher the eta values, the more complex the fingerprint and longer the time will take for the computation of pose fingerprint.
- **feature-file prefix**: prefix to pose fingerprint file
- **class-score file prefix**: prefix to which classification scores output file

### Example
```
cd $RNAPOSERS_PATH
pdb=2b57
mkdir tests/output/
./src/rna_poser.sh tests/input/${pdb}/ ${pdb} complex.pdb lig_${pdb}.mol2 complexes.dcd 2.5 2 tests/output/${pdb}_features tests/output/${pdb}_class_scores
```
#### Output
```
cat tests/output/class_scores_${pdb}.txt
# columns: prediction probability(0) probability(1)
1.000000 0.004000 0.996000
1.000000 0.001000 0.999000
1.000000 0.004000 0.996000
1.000000 0.002000 0.998000
1.000000 0.001000 0.999000
...      ...      ...     
...      ...      ...     
...      ...      ...     
0.000000 0.963000 0.037000
0.000000 0.967000 0.033000
0.000000 0.975000 0.025000
0.000000 0.983000 0.017000
0.000000 0.984000 0.016000
```

### PyMOL plugin
The RNAPosers PyMOL plugin is compressed it as `rnaposerplugin.zip`. See PyMOL website: https://pymolwiki.org/index.php/Plugin_Manager for an instruction on installing PyMOL plugin from local file. Note that you still have to install this repo and setup environment to use the plugin.


### Additional Notes
- RNAPosers expect that the ordering of ligand atoms in the reference receptor-ligand pdb matches the order in the reference ligand mol2 file
- RNAPosers reads atom types from the reference ligand mol2 file
- Atom types should to be one of the following SYBYL atom types:

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
