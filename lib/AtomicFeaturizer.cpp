//
//  AtomicFeaturizer.cpp
//  
//
//  Created by J.X on 9/29/17.
//
//  Pending: tests on revised -rowatm and -selatm functionals

#include "Coor.hpp"
#include "Atom.hpp"
#include "Analyze.hpp"
#include "Molecule.hpp"
#include "Select.hpp"
#include "AtomicFeaturizer.hpp"


// value ctor
AtomicFeaturizer::AtomicFeaturizer(Molecule *mol){
  this->mol = mol;
  
  mol->selAll();
  
  // set atminx for easier 2-D lookup tables
  mol->assignAtmInx();
}

// another value ctor: build from a pdb or mol2 file
AtomicFeaturizer::AtomicFeaturizer(string pdb){
  if (pdb.find(".pdb") != std::string::npos){
    mol = Molecule::readPDB(pdb);
  }
  else{
    cout << "PDB extension not found; try with Mol2." << endl;
    mol = Molecule::readMol2(pdb);
  }
  /* Print out current molecule info */
  cout << "============================================================" << endl;
  string fn = pdb;
  cout << "Input: " << fn << endl;
  cout << "Number of Atoms: " << mol->getAtmVecSize() << endl;
  cout << "Number of Residues: " << mol->getResVecSize() << endl;
  cout << "Number of Chains: " << mol->getChnVecSize() << endl;
  
  mol->selAll();
  
  // set atminx for easier 2-D lookup tables
  mol->assignAtmInx();
}



// featurizer for single eta of single atom
bool AtomicFeaturizer::featurizeEtAtm(const vector<vector<Atom*> >& neighbors,
                                      const vector<vector<double> >& pdin,
                                      const vector<vector<double> >& pdx,
                                      const vector<vector<double> >& pdy,
                                      const vector<vector<double> >& pdz,
                                      double cutoff, double eta, Atom *ai,
                                      string ngh_select, vector<double>& Atmfeatures){
  // ngh_select: NEIGHBOR nucleus to be included in features
  
  
  // initialize parameters
  double dx;
  double dy;
  double dz;
  double multiplier;
  double dist;
  
  Atom* aj;
  int atmInx;
  int ajAtmInx;
  string chain = "";
  string resname = "";
  string atmname = "";
  double xfeature = 0.0;
  double yfeature = 0.0;
  double zfeature = 0.0;
  bool isamatch = false;
  bool selempty = false;
  
  // basic information of current atom
  atmInx = ai->getAtmInx();
  // transform ngh_select into chain, resname & atmname
  size_t pos;
  if (ngh_select.empty()){
    selempty = true;
  }
  else{
    if ((pos=ngh_select.find(":")) != std::string::npos){
      if (pos > 0){
        chain = ngh_select.substr(0, pos);
      }
      ngh_select.erase(0, pos + 1);
    }
    if ((pos=ngh_select.find(".")) != std::string::npos){
      if (pos > 0){
        resname = ngh_select.substr(0, pos);
      }
      ngh_select.erase(0, pos + 1);
    }
  }
  atmname = ngh_select;
  /* calculate features*/
  for (unsigned int j = 0; j < neighbors[atmInx].size(); j++){
    aj = neighbors[atmInx][j];
    isamatch = selempty;
    if (!selempty){
      if (chain.empty()||aj->getChainId() == chain){
        if (resname.empty()||aj->getResName() == resname){
          if (atmname.empty()||aj->getAtmName() == atmname){
            isamatch = true;
          }
        }
      }
    }
    if (isamatch){
      ajAtmInx = aj->getAtmInx();
      dist = pdin[atmInx][ajAtmInx];
      dx = pdx[atmInx][ajAtmInx];
      dy = pdy[atmInx][ajAtmInx];
      dz = pdz[atmInx][ajAtmInx];
      multiplier = 1/ dist * exp((- dist/eta) * (dist/eta)) * (cos(PI * dist/ cutoff) + 1)/2;
      xfeature += dx * multiplier;
      yfeature += dy * multiplier;
      zfeature += dz * multiplier;
    }
  }
  
  Atmfeatures.push_back(xfeature);
  Atmfeatures.push_back(yfeature);
  Atmfeatures.push_back(zfeature);
  return 0;
}



// featurizer for single eta & single selected_neighbor_atomname (if selected)
bool AtomicFeaturizer::featurizeEta(const vector<vector<Atom*> >& neighbors,
                                    const vector<vector<double> >& pdin,
                                    const vector<vector<double> >& pdx,
                                    const vector<vector<double> >& pdy,
                                    const vector<vector<double> >& pdz,
                                    double cutoff, double eta, string ngh_select,
                                    Molecule* selMol,
                                    vector<double>& features){
  // ngh_select: NEIGHBOR atmname to be considered to calculate features;
  //             leave empty when not specified, all atoms will be taken into consideration
  
  // initialize feature vectors
  int natom = selMol->getNAtom();
  vector<double> xfeatures (natom);
  vector<double> yfeatures (natom);
  vector<double> zfeatures (natom);
  vector<double> featureOneAtm;
  
  // initialize parameters
  Atom* ai;
  string atmname;
  
  for (unsigned int i=0; i < natom; i++){
    featureOneAtm.clear();// cleanup for current atom
    ai = selMol->getAtom(i);
    atmname = ai->getAtmName();
    featurizeEtAtm(neighbors, pdin, pdx, pdy, pdz, cutoff, eta, ai, ngh_select, featureOneAtm);
    xfeatures[i] = featureOneAtm[0];
    yfeatures[i] = featureOneAtm[1];
    zfeatures[i] = featureOneAtm[2];
  }
  
  // push back into feature vector in the order of
  // 1x 1y 1z 2x 2y 2z 3x 3y 3z ...
  for (unsigned j = 0; j < natom; j ++){
    features.push_back(xfeatures[j]);
    features.push_back(yfeatures[j]);
    features.push_back(zfeatures[j]);
  }
  return 0;
}




bool AtomicFeaturizer::getNeighAtmList(vector<vector<Atom*> >& neighbors, double cutoff,
                                       const vector<vector<double> >& pdin){
  cout << "Generating neighboring atom list..." << endl;
  Atom* ai;
  Atom* aj;
  int atmInx;
  string atmname;
  int nghatmInx;
  int dist;
  vector<Atom*> neigh;
  
  for (unsigned int i=0; i < mol->getNAtom(); i ++){
    neigh.clear();
    ai = mol->getAtom(i);
    // basic information of current atom
    atmname = ai->getAtmName();
    atmInx = ai->getAtmInx();
    for (unsigned int j=0; j < mol->getNAtom(); j ++){
      aj = mol->getAtom(j);
      nghatmInx = aj->getAtmInx();
      dist = pdin[atmInx][nghatmInx];
      if ((nghatmInx != atmInx) && (dist < cutoff)){
        neigh.push_back(aj);
      }
    }
    neighbors.push_back(neigh);
  }
  cout << "Neighboring atom list done" << endl;
  return 0;
}



// featurizer - vectorial version
bool AtomicFeaturizer::featurize(double cutoff, vector<double> etaList,
                                 vector<vector<double> >& features, string fname,
                                 vector<string>& sel_atmname,
                                 vector<string>& row_atmname){
  
  // declaration of temporary parameters
  stringstream ss;
  string selKeyAtm;
  vector<double> featuresOneEta;
  vector<vector<Atom*> > neighbors;
  Molecule* selMol;
  int nAtomSel;
  
  // open file to write
  ofstream outFile;
  outFile.open(fname.c_str());
  
  // column label
  outFile << "atmid  resname  atmname  ";
  
  /* get pairwise distance matrix*/
  vector<vector<double> > pdin;
  vector<vector<double> > pdx;
  vector<vector<double> > pdy;
  vector<vector<double> > pdz;
  
  // assign atom indices and pairwise distance matrix
  mol->assignAtmInx();
  Analyze::pairwiseDistance(mol, pdin);
  cout << "Pairwise distance matrix done " << endl;
  
  // pairwise distance component vectors
  Analyze::pairwisedv('x', mol, pdx);
  Analyze::pairwisedv('y', mol, pdy);
  Analyze::pairwisedv('z', mol, pdz);
  
  // select row atoms to be featurized
  if (!row_atmname.empty()){
    for (unsigned int j = 0; j < row_atmname.size() - 1; j ++){
      ss << row_atmname[j] << "_";
    }
    ss << row_atmname[row_atmname.size() - 1];
    selKeyAtm = ss.str();
    cout << "Select atom list: " << selKeyAtm << endl;
    mol->select(selKeyAtm);
    selMol = mol->copy();// make a copy of selected row atoms
  }
  else{
    mol->selAll();
    selMol = mol;
  }
  nAtomSel = selMol->getNAtom();
  cout << "Number of atoms selected: " << nAtomSel  << endl;
  
  getNeighAtmList(neighbors, cutoff, pdin);
  if (sel_atmname.empty()){
    // if no atom type specified: equally sum over all neighboring atoms (within cutoff)
    for (unsigned int j=0; j< etaList.size(); j++){
      outFile << "eta=" << etaList[j] << " ";//column labels
      featuresOneEta.clear();//clean up at each loop
      cout << "eta " << j << " = " << etaList[j] << endl;
      featurizeEta(neighbors, pdin, pdx, pdy, pdz, cutoff, etaList[j], "",
                   selMol, featuresOneEta);
      features.push_back(featuresOneEta);
    }
  }
  else{
    // if sel_atmname not empty: only neighboring atoms belong to sel_atmname will be
    // considered for featurization
    // one feature for each type of neighboring atoms
    for (unsigned int j=0; j< etaList.size(); j++){
      cout << "eta " << j << " = " << etaList[j] << endl;
      for (unsigned int i=0; i< sel_atmname.size(); i++){
        outFile << etaList[j] << "," << sel_atmname[i] << " ";//column labels
        vector<double> featuresOneD;
        featurizeEta(neighbors, pdin, pdx, pdy, pdz, cutoff, etaList[j],
                     sel_atmname[i], selMol,featuresOneD);
        features.push_back(featuresOneD);
      }
    }
  }
  outFile << "comp" << endl;
  
  // selected row atoms
  vector<Atom*> atomList = selMol->getAtmVec();
  char vec3d [] = {'x', 'y', 'z'};
  int atmnum;
  string resname;
  string atmname;
  
  for (unsigned int j = 0; j < features[0].size(); j ++){
    atmnum = atomList[j/3]->getAtmNum();
    resname = atomList[j/3]->getResName();
    atmname = atomList[j/3]->getAtmName();
    
    outFile << atmnum << "  " << resname << "  " << atmname << "  ";
    for (unsigned int i = 0; i < features.size(); i ++){
      outFile << features[i][j] << " ";
    }
    outFile << vec3d[ j % 3];
    outFile << endl;
  }
  return 0;
}



/****************************************
 
		Scalar version of featurizer
 
 ****************************************/


// featurizer for single eta of single atom
double AtomicFeaturizer::featurizeEtAtm(const vector<vector<Atom*> >& neighbors,
                                      const vector<vector<double> >& pdin,
                                      double cutoff, double eta, Atom *ai,
                                      string ngh_select){
  // ngh_select: NEIGHBOR nucleus to be included in features
  
  
  // initialize parameters

  double dist;
  
  Atom* aj;
  int atmInx;
  int ajAtmInx;
  string chain = "";
  string resname = "";
  string atmname = "";
  double feature = 0.0;
  bool isamatch = false;
  bool selempty = false;
  
  // basic information of current atom
  atmInx = ai->getAtmInx();
  // transform ngh_select into chain, resname & atmname
  size_t pos;
  if (ngh_select.empty()){
    selempty = true;
  }
  else{
    if ((pos=ngh_select.find(":")) != std::string::npos){
      if (pos > 0){
        chain = ngh_select.substr(0, pos);
      }
      ngh_select.erase(0, pos + 1);
    }
    if ((pos=ngh_select.find(".")) != std::string::npos){
      if (pos > 0){
        resname = ngh_select.substr(0, pos);
      }
      ngh_select.erase(0, pos + 1);
    }
  }
  atmname = ngh_select;
  /* calculate features*/
  for (unsigned int j = 0; j < neighbors[atmInx].size(); j++){
    aj = neighbors[atmInx][j];
    isamatch = selempty;
    if (!selempty){
      if (chain.empty()||aj->getChainId() == chain){
        if (resname.empty()||aj->getResName() == resname){
          if (atmname.empty()||aj->getAtmName() == atmname){
            isamatch = true;
          }
        }
      }
    }
    if (isamatch){
      ajAtmInx = aj->getAtmInx();
      dist = pdin[atmInx][ajAtmInx];
      feature += exp((- dist/eta) * (dist/eta)) * (cos(PI * dist/ cutoff) + 1)/2;
    }
  }
  return feature;
}



// featurizer for single eta & single selected_neighbor_atomname (if selected)
bool AtomicFeaturizer::featurizeEta(const vector<vector<Atom*> >& neighbors,
                                    const vector<vector<double> >& pdin,
                                    double cutoff, double eta, string ngh_select,
                                    Molecule* selMol,
                                    vector<double>& features){
  // ngh_select: NEIGHBOR atmname to be considered to calculate features;
  //             leave empty when not specified, all atoms will be taken into consideration
  
  // initialize feature vectors
  int natom = selMol->getNAtom();
  double featureOneAtm;
  
  // initialize parameters
  Atom* ai;
  string atmname;
  
  features.clear();
  
  // get feature of each atom and push back into the features vector
  for (unsigned int i=0; i < natom; i++){
    ai = selMol->getAtom(i);
    atmname = ai->getAtmName();
    featureOneAtm = featurizeEtAtm(neighbors, pdin, cutoff, eta, ai, ngh_select);
    features.push_back(featureOneAtm);
  }
  
  return 0;
}


// featurizer - scalar
bool AtomicFeaturizer::featurizeScalar(double cutoff, vector<double> etaList,
                                 vector<vector<double> >& features, string fname,
                                 vector<string>& sel_atmname,
                                 vector<string>& row_atmname, bool verbose){
  
  
  // declaration of temporary parameters
  stringstream ss;
  string selKeyAtm;
  vector<double> featuresOneEta;
  vector<vector<Atom*> > neighbors;
  Molecule* selMol;
  int nAtomSel;
  

  // open file to write
  ofstream outFile;
  
  if (verbose){
    outFile.open(fname.c_str());
    // column label
    outFile << "atmid  resname  atmname  ";
  }
  
  /* get pairwise distance matrix*/
  vector<vector<double> > pdin;
  
  // assign atom indices and pairwise distance matrix
  mol->assignAtmInx();
  Analyze::pairwiseDistance(mol, pdin);
  cout << "Pairwise distance matrix done " << endl;
  
  // select row atoms to be featurized
  if (!row_atmname.empty()){
    for (unsigned int j = 0; j < row_atmname.size() - 1; j ++){
      ss << row_atmname[j] << "_";
    }
    ss << row_atmname[row_atmname.size() - 1];
    selKeyAtm = ss.str();
    cout << "Select atom list: " << selKeyAtm << endl;
    mol->select(selKeyAtm);
    selMol = mol->copy();// make a copy of selected row atoms
  }
  else{
    mol->selAll();
    selMol = mol;
  }
  nAtomSel = selMol->getNAtom();
  cout << "Number of atoms selected: " << nAtomSel  << endl;
  
  getNeighAtmList(neighbors, cutoff, pdin);
  if (sel_atmname.empty()){
    // if no atom type specified: equally sum over all neighboring atoms (within cutoff)
    for (unsigned int j=0; j< etaList.size(); j++){
      if (verbose){
        outFile << "eta=" << etaList[j] << " ";//column labels
      }
      featuresOneEta.clear();//clean up at each loop
      cout << "eta " << j << " = " << etaList[j] << endl;
      featurizeEta(neighbors, pdin, cutoff, etaList[j], "", selMol, featuresOneEta);
      features.push_back(featuresOneEta);
    }
  }
  else{
    // if sel_atmname not empty: only neighboring atoms belong to sel_atmname will be
    // considered for featurization
    // one feature for each type of neighboring atoms
    for (unsigned int j=0; j< etaList.size(); j++){
      cout << "eta " << j << " = " << etaList[j] << endl;
      for (unsigned int i=0; i< sel_atmname.size(); i++){
        if (verbose){
          outFile << etaList[j] << "," << sel_atmname[i] << " ";//column labels
        }
        featuresOneEta.clear();
        featurizeEta(neighbors, pdin, cutoff, etaList[j], sel_atmname[i], selMol,featuresOneEta);
        features.push_back(featuresOneEta);
      }
    }
  }
  
  /* write features to file */
  if (verbose){
    outFile << endl;
    // selected row atoms
    vector<Atom*> atomList = selMol->getAtmVec();
    int atmnum;
    string resname;
    string atmname;
  
    for (unsigned int j = 0; j < features[0].size(); j ++){
      atmnum = atomList[j]->getAtmNum();
      resname = atomList[j]->getResName();
      atmname = atomList[j]->getAtmName();
      
      outFile << atmnum << " " << resname << " " << atmname << " ";
      for (unsigned int i = 0; i < features.size(); i ++){
        outFile << features[i][j] << " ";
      }
      outFile << endl;
    }
  }
  return 0;
}
