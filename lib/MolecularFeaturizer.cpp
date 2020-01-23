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

// bool check mol2 file provided
bool checkMol2(std::vector<std::string> mol2s){
  if (!mol2s.size()){
    cout << endl <<  "Error: Please provide a mol2 file. [-mol2 MOL2file] " << endl << endl;
    return false;
  }
  return true;
}

// process mol2 file
Molecule* processMol2(std::string mol2file,
                      vector<string> row_atmname,
                      vector<string>& atmType){

  Molecule *mol2=NULL;
  mol2 = Molecule::readMol2(mol2file);
  mol2->selAll();

  if (!row_atmname.empty()){
    stringstream ss;
    string selKeyAtm;
    for (unsigned int j = 0; j < row_atmname.size() - 1; j ++){
      ss << row_atmname[j] << "_";
    }
    ss << row_atmname[row_atmname.size() - 1];
    selKeyAtm = ss.str();
    mol2->select(selKeyAtm);
  }
  
  // getAtmtypes from mol2
  atmType.clear();
  for (unsigned int i=0; i< mol2->getAtmVecSize(); i++){
    if (mol2->getAtom(i)->getSel() == true){
      atmType.push_back(mol2->getAtom(i)->getAtmType());
    }
  }
  return mol2;
}


// Transfer atomic features into molecular features
bool MolecularFeaturizer(
    Molecule* molall,
    Molecule* mol2all,
    std::vector<std::vector<double> >& features,
    std::vector<std::vector<double> >& mf_traj,
    string fname,
    const std::vector<string>& atmTypes,
    const std::vector<string>& SYBYL,
    bool verbose)
{
  Molecule* mol = molall->clone();
  Molecule* mol2 = mol2all->clone();
  cout << "Generating molecular features..." << endl;
  for (int k = 0; k < mol2->getNAtom(); k ++){
    if ((mol->getAtom(k)->getResName() != mol2->getAtom(k)->getResName()) || (mol->getAtom(k)->getAtmName() != mol2->getAtom(k)->getAtmName())){
      cout << "Mol2 atoms do not match with selected molecule!" << endl;
	  if ((mol->getAtom(k)->getResName() != mol2->getAtom(k)->getResName())){
		cout << "Residue name mismatch at atom " << k << ": " << mol->getAtom(k)->getResName() << " and " << mol2->getAtom(k)->getResName() << endl;
	  }
	  if ((mol->getAtom(k)->getAtmName() != mol2->getAtom(k)->getAtmName())){
		cout << "Atom name mismatch at atom " << k << ": " << mol->getAtom(k)->getAtmName() << " and " << mol2->getAtom(k)->getAtmName() << endl;
	  }
      return false;
    }
  }
  
  if (mol2->getNAtom() != mol->getNAtom() ){
    cout << "mol2 size " << mol2->getNAtom() << " does not match mol size " << mol->getNAtom() << endl;
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
  std::vector<double> mfeats;
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
  
  for (i = 0; i < mol_features.size(); i ++){
    for (j = 0; j < mol_features[0].size(); j ++){
      mfeats.push_back(mol_features[i][j]);
    }
  }
  mf_traj.push_back(mfeats);
  
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

