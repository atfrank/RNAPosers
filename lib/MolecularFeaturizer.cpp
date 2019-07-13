//
//  MolecularFeaturizer.cpp
//  
//
//  Created by J.X on 4/9/18.

#include "Coor.hpp"
#include "Atom.hpp"
#include "Analyze.hpp"
#include "Molecule.hpp"
#include "Select.hpp"
#include "AtomicFeaturizer.hpp"
#include "MolecularFeaturizer.hpp"

// getAtmtypes
Molecule* getAtmType(string mol2file, vector<string>& atmType){
  Molecule* mol2 = Molecule::readMol2(mol2file);
  std::vector<Atom*> atmVec = mol2->getAtmVec();
  for (int k = 0; k < atmVec.size(); k ++){
    atmType.push_back(atmVec.at(k)->getAtmType());
  }
  return mol2;
}

// Transfer atomic features into molecular features
bool MolecularFeaturizer(
    Molecule* mol,
    Molecule* mol2,
    std::vector<std::vector<double> >& features,
    std::vector<std::vector<double> >& mf_traj,
    string fname,
    const std::vector<string>& atmTypes,
    const std::vector<string>& SYBYL,
    std::vector<double>& mfeats, bool verbose)
{
  cout << "Generating molecular features..." << endl;
  if (features[0].size() != atmTypes.size()){
    cout << "Mol2 molecule size does not match with selected molecule! " << endl;
    return false;
  }
  

  Molecule* molc = mol->clone();
  cout << molc->getAtmVecSize() << endl;

  std::vector<Atom*> atmVec = molc->getAtmVec();
  std::vector<Atom*> atmVec2 = mol2->getAtmVec();
  
  
  for (int k = 0; k < atmVec.size(); k ++){
    if ((atmVec.at(k)->getResName() != atmVec2.at(k)->getResName()) || (atmVec.at(k)->getAtmName() != atmVec2.at(k)->getAtmName())){
      cout << "Mol2 atoms do not match with selected molecule!" << endl;
	  if ((atmVec.at(k)->getResName() != atmVec2.at(k)->getResName())){
		cout << "Residue name mismatch at atom " << k << ": " << atmVec.at(k)->getResName() << " and " << atmVec2.at(k)->getResName() << endl;
	  }
	  if ((atmVec.at(k)->getAtmName() != atmVec2.at(k)->getAtmName())){
		cout << "Atom name mismatch at atom " << k << ": " << atmVec.at(k)->getAtmName() << " and " << atmVec2.at(k)->getAtmName() << endl;
	  }
      return false;
    }
  }
  
  unsigned int i, j, k;
  
  
  // transpose features
  vector<vector<double> > transfeatures(features[0].size(),
                                  vector<double>(features.size()));
  for (i = 0; i < features.size(); ++i){
    for (j = 0; j < features[0].size(); ++j){
      transfeatures[j][i] = features[i][j];
    }
  }
  
  std::vector<std::vector<double> > mol_features;
  mfeats.clear();
  
  // define Natom = Number of atoms
  int Natom = atmTypes.size();
  
  for (i=0; i< SYBYL.size(); i++){
    // initialize mol_features to all zeros
    mol_features.push_back(std::vector<double> (transfeatures[0].size(), 0.0));
  }

  for (i=0; i< SYBYL.size(); i++){
    for (j=0; j< Natom; j++){
      if (atmTypes[j] == SYBYL[i]){
        // add features to mol_features if atmType match
        for (k=0; k < mol_features[0].size(); k++){
          mol_features[i][k] += transfeatures[j][k];
        }
      }
    }
  }

  if (verbose){
    // open file to write
    ofstream outFile;
    outFile.open(fname.c_str());
    for (i = 0; i < mol_features.size(); i ++){
      for (j = 0; j < mol_features[0].size(); j ++){
        outFile << mol_features[i][j] << " ";
      }
    }
  }
  
  for (i = 0; i < mol_features.size(); i ++){
    for (j = 0; j < mol_features[0].size(); j ++){
      mfeats.push_back(mol_features[i][j]);
    }
  }
  mf_traj.push_back(mfeats);
  cout << "Molecular features done." << endl;
  return true;

}

// write traj molecular features to file
bool writeMFtraj(std::vector<std::vector<double> >& mf_traj, string fname, vector<int>& frames){
  ofstream outFile;
  outFile.open(fname.c_str());
  if (mf_traj.size() != frames.size()){
    cout << "ERROR: Number of frames in traj and feature length do not match!" << endl;
    return 0;
  }
  for (unsigned int i=0; i < mf_traj.size(); i++){
    outFile << frames[i] << " ";
    for (unsigned int j=0; j < mf_traj[0].size(); j++){
      outFile << mf_traj[i][j] << " ";
    }
    outFile << std::endl;
  }
  return true;
}

// normalize mf_traj by column: divide each column by its median
bool normalize(std::vector<std::vector<double> >& mf_traj){

  unsigned int i, j;
  double median;
  int numcols, veclen;
  numcols = mf_traj[0].size();
  veclen = mf_traj.size();
  
  // transpose features
  vector<vector<double> > trans_mf(mf_traj[0].size(),
                                        vector<double>(mf_traj.size()));
  for (i = 0; i < mf_traj.size(); ++i){
    for (j = 0; j < mf_traj[0].size(); ++j){
      trans_mf[j][i] = mf_traj[i][j];
    }
  }

  for (i = 0; i < numcols; ++i){
    sort(trans_mf[i].begin(), trans_mf[i].end());
    // find median
    if (veclen % 2 == 0){
      median = (trans_mf[i][veclen / 2 - 1] + trans_mf[i][veclen / 2]) / 2;
    }
    else{
      median = trans_mf[i][veclen / 2];
    }
    if (median != 0){
      for (j = 0; j < veclen; ++j){
        // divide the whole column in the original vector matrix by column median
        mf_traj[j][i] = mf_traj[j][i] / median;
      }
    }
  }
  return 1;
}

