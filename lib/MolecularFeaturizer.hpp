//
//  MolecularFeaturizer.hpp
//
//
//  Created by Jingru Xie on 4/9/18.
//
//

#include <iostream>
#include <fstream>

#include <stdio.h>
#include <sstream>
#include <math.h>
#include <algorithm>

#define PI 3.14159265

using namespace std;

// getAtmtypes from mol2 file
Molecule* getAtmType(
  string mol2file, std::vector<string>& atmType
  );

// reshape atomic features into molecular features
bool MolecularFeaturizer(
  Molecule* mol,
  Molecule* mol2,
  std::vector<std::vector<double> >& features,
  std::vector<std::vector<double> >& mf_traj,
  string fname,
  const std::vector<string>& atmTypes,
  const std::vector<string>& SYBYL,
  std::vector<double>& mfeats,
  bool verbose = true
  );

bool writeMFtraj(vector<vector<double> >& mf_traj,
                 string outFile,
                 vector<int>& frames);

// normalize mf_traj by column: divide each column by its median
bool normalize(std::vector<std::vector<double> >& mf_traj);
