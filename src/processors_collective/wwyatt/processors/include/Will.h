

/*
 * Created by William Wyatt
 * On Aug 30th 2017
 * Make some vector utilities.
 */
#ifndef WILLIAMS_FUN_TIME
#define WILLIAMS_FUN_TIME 1

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <cmath>

#include "lcio.h"
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>

using namespace lcio;
using namespace std;
static vector<MCParticle *> AVAL;
class Will{
 public:
  //Specific function used in prediction algorithm.
  //Finds the highest energy particle
  static map<int,double> maxEnergy(LCCollection*, 
				   initializer_list<int> ids, 
				   vector<MCParticle*>& final_state);

  static double* getVector(MCParticle*);

  //Returns the sum of the two; assumes a 4 vector
  static double* addVector(double*, double*, const int SIZE=4);

  //Returns transverse momentum magnitude
  static double getTMag(double*);

  static double getMag(double*);

 private:

};

#endif
