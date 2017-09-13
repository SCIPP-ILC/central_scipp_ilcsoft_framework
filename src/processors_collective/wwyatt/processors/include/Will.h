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
namespace Will{
  struct fourvec{
    union{ double X; double x=0.0; };
    union{ double Y; double y=0.0; };
    union{ double Z; double z=0.0; };
    union{ double E; double e=0.0; };
    union{ double T; double t=0.0; };
    fourvec operator+(const fourvec& a) const{
      fourvec ret;
      ret.x=a.x+x;
      ret.y=a.y+y;
      ret.z=a.z+z;
      ret.t=a.t+t;
      return ret;
    }
    fourvec operator+=(const fourvec& a){
      *this=a+*this;
      return *this;
    }
  };
  struct prediction{
    fourvec hadronic;
    fourvec electronic;
    fourvec electron;
    fourvec positron;
    double mag=0.0;
    bool scattered=false;
  };
  //Specific function used in prediction algorithm.
  //Finds the highest energy particle
   map<int,double> maxEnergy(LCCollection*, 
			     initializer_list<int> ids, 
			     vector<MCParticle*>& final_state);

   double* getVector(MCParticle*);
   fourvec getFourVector(MCParticle*);

  //Returns the sum of the two; assumes a 4 vector
   double* addVector(double*, double*, const int SIZE=4);

  //Returns transverse momentum magnitude
   double getTMag(const double*);

  //Returns momentum from a momentum vector
   double getMag(const double*);
   double getMag(const fourvec);

   double* legacy(fourvec);


  /*Returns a map with a few four vectors in it.
   * - hadronic vector
   * - electronic vector
   * - electron vector
   * - positron vector
   * This should be used to calculate a prediction vector.
   */
   prediction getPrediction(LCCollection*);
}

#endif
