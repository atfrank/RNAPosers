//
//  AtomicFeaturizer.hpp
//
//
//  Created by Jingru Xie on 9/29/17.
//
//

#include <iostream>
#include <fstream>

#include <stdio.h>
#include <sstream>
#include <math.h>
#include <algorithm>

#ifndef AtomicFeaturizer_hpp
#define AtomicFeaturizer_hpp

#define PI 3.14159265

using namespace std;

class Molecule;

class AtomicFeaturizer{
public:
  AtomicFeaturizer();//default constructor
  AtomicFeaturizer(Molecule *mol=NULL);//value ctor
  AtomicFeaturizer(string pdb);//another value ctor
  bool getNeighAtmList(vector<vector<Atom*> >& neighbors, double cutoff,
                       const vector<vector<double> >& pdin);
  bool featurizeEtAtm(const vector<vector<Atom*> >& neighbors,
                      const vector<vector<double> >& pdin,
                      const vector<vector<double> >& pdx,
                      const vector<vector<double> >& pdy,
                      const vector<vector<double> >& pdz,
                      double cutoff, double eta, Atom *ai,
                      string ngh_select, vector<double>& Atmfeatures);
  bool featurizeEta(const vector<vector<Atom*> >& neighbors,
                    const vector<vector<double> >& pdin,
                    const vector<vector<double> >& pdx,
                    const vector<vector<double> >& pdy,
                    const vector<vector<double> >& pdz,
                    double cutoff, double eta, string ngh_select,
                    Molecule* selMol,
                    vector<double>& features);
  bool featurize(double cutoff, vector<double> etaList, vector<vector<double> >& features, string fname, vector<string>& sel_atmname, vector<string>& row_atmname);
  
  // scalar overload of above functions
  double featurizeEtAtm(const vector<vector<Atom*> >& neighbors,
                        const vector<vector<double> >& pdin,
                        double cutoff, double eta, Atom *ai,
                        string ngh_select);
  bool featurizeEta(const vector<vector<Atom*> >& neighbors,
                    const vector<vector<double> >& pdin,
                    double cutoff, double eta, string ngh_select,
                    Molecule* selMol,
                    vector<double>& features);
  bool featurizeScalar(double cutoff, vector<double> etaList,
          						 vector<vector<double> >& features, string fname,
                       vector<string>& sel_atmname, vector<string>& row_atmname,
                       bool verbose=true);

private:
  Molecule* mol;
};


#endif /* AtomicFeaturizer_hpp */
