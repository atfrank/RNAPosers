# ssMD Unbinding Profiles
RNAPosers-ssMD: A method that use a ML regression to predict the ligand unbinding profile.

This folder contains a script that predicts ssMD unbinding profile directly from structure using features similar to those used in this work, Chhabra, Sahil, Jingru Xie, and Aaron T. Frank. "RNAPosers: Machine Learning Classifiers For RNA-Ligand Poses." J. Phys. Chem. B 2020, 124, 22, 4436–4445; Publication Date:May 19, 2020; https://doi.org/10.1021/acs.jpcb.0c02322. 

Accompanies the manuscript: Yichen Liu and Aaron T. Frank. "Instigating Unbinding Events to Identify Near-Native RNA-Ligand Poses" in preparation. 

## Prerequisite
* RNAPosers
* R
```
conda install -c r r=3.6.0 r-pls r-optparse -m -n my-r
```

### Using RNAPosers-ssMD
```
conda activate my-r
Rscript rna_poser_ssMD.R -h
```
#### Arguments
- **pretrained PLS model**: e.g.: final_model.RData
- **pose feature file**: e.g.: example/features.txt

#### Options
- **output**: e.g.: examples/predicted_profile.txt

### Example 1: From Structure
```
export RNAPOSERS_PATH=/Users/aaronfranklab/Documents/GitHub/RNAPosers
export PYTHONPATH=$RNAPOSERS_PATH/py:$PYTHONPATH
```

#### (1) generate features
```
cd $RNAPOSERS_PATH/ssMD
conda activate rnaposers
python $RNAPOSERS_PATH/ssMD/rna_poser_ssMD.py example/receptor.mol2 example/poses.sd example/test_features 10
conda deactivate
```

#### (2) predict profile from features
```
conda activate my-r
Rscript $RNAPOSERS_PATH/ssMD/rna_poser_ssMD.R final_model.RData example/test_features.txt -o example/predicted_profile_1.txt
conda deactivate
```

### Example 2: From Feature Files
```
conda activate my-r
Rscript R/rna_poser_ssMD.R final_model.RData example/features.txt -o examples/predicted_profile_2.txt
conda deactivate
```
### Output
```
cat examples/predicted_profile_1.txt
# columns: pose_number A tau B
1 0.8310509 10.26660 0.001373989
2 0.8315107 10.24533 0.001372110
3 0.8195701 10.83536 0.001430587
4 0.8210322 10.76239 0.001423038
5 0.8308009 10.27905 0.001375514
6 0.8255692 10.54354 0.001400799
...
```
## COMMERCIAL USE LICENSE:

If you are interested in commercial licensing of these applications (clinical, operational, etc.) please contact the University of Michigan Office of Technology Transfer for a quote and licensing options.

Drew Bennett - https://techtransfer.umich.edu/team/drew-bennett/

or

techtransfer@umich.edu
