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

// check if mol2 file provided
bool checkMol2(std::vector<std::string> mol2s);

// process mol2 file
Molecule* processMol2(std::string mol2file,
                      vector<string> row_atmname,
                      vector<string>& atmType);

// reshape atomic features into molecular features
bool MolecularFeaturizer(
  Molecule* mol,
  Molecule* mol2,
  std::vector<std::vector<double> >& features,
  std::vector<std::vector<double> >& mf_traj,
  string fname,
  const std::vector<string>& atmTypes,
  const std::vector<string>& SYBYL,
  bool verbose = true
  );

bool writeMFtraj(vector<vector<double> >& mf_traj,
                 string outFile,
                 vector<int>& frames);

// normalize mf_traj by column: divide each column by its median
bool normalize(std::vector<std::vector<double> >& mf_traj);
