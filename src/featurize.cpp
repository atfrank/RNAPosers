/*

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.

  Author: Jingru Xie and Aaron T. Frank

*/

#include "Molecule.hpp"
#include "Residue.hpp"
#include "Atom.hpp"
#include "Misc.hpp"
#include "Analyze.hpp"
#include "Trajectory.hpp"
#include "AtomicFeaturizer.hpp"
#include "MolecularFeaturizer.hpp"


#include <iostream>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <math.h>
#include <algorithm>
#include <time.h> // keep track of processing time of program

using namespace std;

void usage(){
  std::cerr << "====================================================" << std::endl;
  std::cerr << "ATOMIC Featurizer v1.00" << std::endl;
  std::cerr << "(c) 2017 Jingru Xie, Aaron T. Frank and University of Michigan." << std::endl;
  std::cerr << "====================================================" << std::endl;
  std::cerr << "Usage:   featurize [-options] <PDBfile>" << std::endl;
  std::cerr << "Options: [-cutoff CUToff]" << std::endl;
  std::cerr << "         [-numEta NumEta]" << std::endl;
  std::cerr << "         [-etaBase etaStartPow]" << std::endl;
  std::cerr << "         [-etaStartPow etaStartPow]" << std::endl;
  std::cerr << "         [-selatm sel_atmname]" << std::endl;
  std::cerr << "         [-rowatm row_atmname (e.g. ':UNK.', ':UNK/2.')]" << std::endl;
  std::cerr << "         [-scalar scalar option: enter 1 for scalar featurizer]" << std::endl;
  std::cerr << "         [-molecular molecular option: enter 1 for molecular featurizer]" << std::endl;
  std::cerr << "         [-normalization normalization option: enter 1 to enable feature normalization. Default: 0. This option only works for trajectory molecular featurization]" << std::endl;
  std::cerr << "         [-outfile path and name (without extension) of output feature file. If multiple trajs or pdbs are provided, the frame id will automatically be appended to filename]" << std::endl;
  std::cerr << "         [-mol2 MOL2file]" << std::endl;
  std::cerr << "         [-trj TRAJfile]" << std::endl;
  std::cerr << "         [-skip frames] [-start frame] [-stop frame]" << std::endl;
  std::cerr << "         [-identification ID]" << std::endl;
  std::cerr << std::endl;
  exit(0);
}

int main (int argc, char **argv){
  int i;
  unsigned int f;
  clock_t t_start, t_end;

  t_start = clock();

  std::stringstream resid;
  std::vector<std::string> pdbs;
  std::vector<std::string> mol2s;
  std::string currArg;
  std::string outf;
  std::string outfile;

  std::vector<std::string> trajs;
  int start = 0;
  int stop=std::numeric_limits<int>::max();
  int skip = 0;
  bool startFlag=false;
  unsigned long long int itrj;
  std::ifstream trjin;
  Trajectory *ftrjin;
  unsigned long long int nframe = 0;

  bool scalar = false;
  bool molecular = false;
  bool normalization = false;
  double cutoff;
  double eta;
  int numEta;
  int etaStartPow;
  double etaBase;
  string select_atm;
  vector<string> sel_atmname;
  vector<string> row_atmname;

  cutoff = 5.;
  numEta = 8;
  etaStartPow = -1;
  etaBase = 2.0;

  Molecule *neighbormol;
  neighbormol=NULL;

  Atom *ai, *aj;
  ai=NULL;
  aj=NULL;
  outfile="features";

  pdbs.clear();

  for (i = 1; i < argc; i++){
    currArg = argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0)
    {
      usage();
    }
    else if (currArg.compare("-cutoff") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> cutoff;
    }
    else if (currArg.compare("-numEta") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> numEta;
    }
    else if (currArg.compare("-etaBase") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> etaBase;
    }
    else if (currArg.compare("-etaStartPow") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> etaStartPow;
    }
    else if (currArg.compare("-selatm") == 0){
      currArg=argv[++i];
      select_atm = currArg;
      size_t pos = 0;
      string atm;
      string delim = " ";
      while ((pos = select_atm.find(delim)) != string::npos){
        atm = select_atm.substr(0, pos);
        sel_atmname.push_back(atm);
        select_atm.erase(0, pos + delim.length());
      }
      sel_atmname.push_back(select_atm);
    }
    else if (currArg.compare("-rowatm") == 0){
      currArg=argv[++i];
      select_atm = currArg;
      size_t pos = 0;
      string atm;
      string delim = " ";
      while ((pos = select_atm.find(delim)) != string::npos){
        atm = select_atm.substr(0, pos);
        row_atmname.push_back(atm);
        select_atm.erase(0, pos + delim.length());
      }
      row_atmname.push_back(select_atm);
    }
    else if (currArg.compare("-scalar") == 0)
    {
      currArg=argv[++i];
      stringstream(currArg) >> scalar;
    }
    else if (currArg.compare("-molecular") == 0)
    {
      currArg=argv[++i];
      stringstream(currArg) >> molecular;
    }
    else if (currArg.compare("-normalization") == 0)
      {
      currArg=argv[++i];
      stringstream(currArg) >> normalization;
      }
    else if (currArg.compare("-outfile") == 0)
    {
      currArg=argv[++i];
      outfile = currArg;
    }
    else if (currArg.compare("-mol2") == 0)
    {
      currArg=argv[++i];
      mol2s.push_back(currArg);
    }
    else if (currArg.compare("-trj") == 0 || currArg.compare("-traj") == 0)
    {
      currArg=argv[++i];
      trajs.push_back(currArg);
    }
    else if (currArg.compare("-skip") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> skip;
    }
    else if (currArg.compare("-start") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> start;
      start--;
      startFlag=true;
    }
    else if (currArg.compare("-stop") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> stop;
    }
    else if (currArg.compare(0,1,"-") == 0)
    {
      std::cerr << "Warning: Skipping unknown option \"" << currArg << "\"" << std::endl;
    }
    else{
      pdbs.push_back(currArg);
    }
  }
  if (pdbs.size() == 0)
  {
    std::cerr << std::endl << "Error: Please provide an input file" << std::endl << std::endl;
    usage();
  }


  //initialize
  Molecule *mol=NULL;
  Molecule *mol2=NULL;
  std::vector<Atom*> atmVec;
  AtomicFeaturizer *atomicfeature;
  atomicfeature=NULL;
  std::vector<string> atmType;

  std::vector<double> etalist;
  std::vector<std::vector<double> > features;
  std::vector<std::vector<double> > mol_features;
  std::vector<int> frames; // frameids being featurized in trajectory analysis

  // initialize SYBYL sequence
  std::string str[] = {"C.1","C.2","C.3","C.ar","C.cat","H","N.1","N.2","N.3","N.4","N.am","N.ar","N.pl3","O.2","O.3","O.co2","P.3","S.2","S.3","S.o","S.o2"};
  const std::vector<string> SYBYL(str, str+(sizeof(str)/sizeof(string)));

  //create list of eta values
  for (int i = 0; i < numEta; i ++){
    eta = pow(etaBase, i + etaStartPow);
    etalist.push_back(eta);
  }

  if (trajs.size() > 0){
    if (pdbs.size() > 1){
      std::cerr << std::endl << "Warning: Only the first PDB structure is used for trajectory analysis" << std::endl << std::endl;
    }
    /* Trajectory analysis */
    mol=Molecule::readPDB(pdbs.at(0));
    mol->selAll();

    // initialize mol2 atmtypes for molecular featurization
    if (molecular){
      if (!checkMol2(mol2s)){
        return 0;
      }
      mol2 = processMol2(mol2s.at(0), row_atmname, atmType);
    }

    /* Process trajectories */
    for (itrj=0; itrj< trajs.size(); itrj++){
      trjin.open(trajs.at(itrj).c_str(), std::ios::binary);
      frames.clear();
      mol_features.clear();
      if (trjin.is_open()){
        ftrjin=new Trajectory;
        ftrjin->setMolecule(mol);
        if (ftrjin->findFormat(trjin) == true){
          ftrjin->readHeader(trjin);
          if (skip > 0 && startFlag == false){
            start=skip;
          }
          /* Print out current molecule info */
          cout << trajs.at(itrj) << endl;
          cout << "Number of Atoms: " << mol->getAtmVecSize() << endl;
          cout << "Number of Residues: " << mol->getResVecSize() << endl;
          cout << "Number of Chains: " << mol->getChnVecSize() << endl;
          cout << "Number of frames in traj: " << ftrjin->getNFrame() << endl;
          /* Loop through desired frames */
          for (i=start; i< ftrjin->getNFrame() && i< stop; i=i+1+skip){
            if (ftrjin->readFrame(trjin, i) == false){ // read frame and setmol here
              std::cerr << "Warning: EOF found before the next frame could be read" << std::endl;
              break;
              }
            nframe++;
            frames.push_back(i + 1);
            cout << endl << "Frame " << i+1 << ":" << endl;
            outf = outfile + "_traj" + to_string(itrj + 1) + "_pdb" + to_string(i+1) + ".txt";
            atomicfeature = new AtomicFeaturizer(mol);
            if (molecular){
              atomicfeature->featurizeScalar(cutoff, etalist, features, outf, sel_atmname, row_atmname, false);
              MolecularFeaturizer(mol, mol2, features, mol_features, outf, atmType, SYBYL, false);
            }
            else if (scalar){
              atomicfeature->featurizeScalar(cutoff, etalist, features, outf, sel_atmname, row_atmname);
            }
            else{
              atomicfeature->featurize(cutoff, etalist, features,
                                       outf, sel_atmname, row_atmname);
            }
            features.clear();
          }
        }
        else{
          std::cerr << "Warning: Skipping unknown trajectory format \"";
          std::cerr << trajs.at(itrj) << "\"" << std::endl;
        }
        if (ftrjin != NULL){
          delete ftrjin;
        }
        if (atomicfeature != NULL){
          delete atomicfeature;
        }
        if (molecular){
          if (normalization){
            normalize(mol_features);
          }
          outf = outfile + "_traj" + to_string(itrj + 1) + ".txt";
          writeMFtraj(mol_features, outf, frames);
        }
      }
      trjin.close();
    }
  }

  else{
    if (molecular){
       /* molecular featurization */
      if (!checkMol2(mol2s)){
        return 0;
      }
      mol2 = processMol2(mol2s.at(0), row_atmname, atmType);
    }
    for (f=0; f< pdbs.size(); f++){
      mol=Molecule::readPDB(pdbs.at(f));
      /* Print out file and atom info */
      atomicfeature = new AtomicFeaturizer(mol);
      cout << "Input: " << pdbs.at(f) << endl;
      cout << "Number of Atoms: " << mol->getAtmVecSize() << endl;
      cout << "Number of Residues: " << mol->getResVecSize() << endl;
      cout << "Number of Chains: " << mol->getChnVecSize() << endl;
      outf = outfile + ".txt";
      cout << "Output: " << outf << endl;

      if (molecular){
        if (!checkMol2(mol2s)){
          return 0;
        }
        if (mol2s.size() > 1){
          mol2 = processMol2(mol2s.at(0), row_atmname, atmType);
        }
        atomicfeature->featurizeScalar(cutoff, etalist, features, outf, sel_atmname, row_atmname, false);

        MolecularFeaturizer(mol, mol2, features, mol_features, outf, atmType, SYBYL);
      }
      else if (scalar){
        atomicfeature->featurizeScalar(cutoff, etalist, features, outf, sel_atmname, row_atmname);
      }
      else{
        atomicfeature->featurize(cutoff, etalist, features,
                                   outf, sel_atmname, row_atmname);
      }
      features.clear();
    }
    if (atomicfeature!=NULL){
      delete atomicfeature;
    }
  }

  t_end = clock();
  float diff ((float)t_end - (float)t_start);
  cout << "Processing time: " << diff/CLOCKS_PER_SEC << " seconds." << endl << endl << endl;

  return 0;
}
