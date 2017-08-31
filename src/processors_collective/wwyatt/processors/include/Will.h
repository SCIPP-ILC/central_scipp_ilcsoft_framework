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
  static map<int,double> maxEnergy(LCCollection* col, initializer_list<int> ids, vector<MCParticle*>& fs){
    const int SIZE=ids.size();
    //double* max=new double[SIZE];
    map<int, double>* max=new map<int,double>;
    for(int i=0; i < col->getNumberOfElements(); ++i){
      MCParticle* particle=dynamic_cast<MCParticle*>(col->getElementAt(i));
      if(particle->getGeneratorStatus()!=1)continue;
      fs.push_back(particle);
      int pid=particle->getPDG();
      for(auto id:ids){
	if(pid==id&&particle->getEnergy()>(*max)[id])
	  (*max)[id]=particle->getEnergy();
      }
    }
    return *max;
  }

  static double* getVector(MCParticle* particle){
    double* output=new double[4];
    const double* mom=particle->getMomentum();
    output[0]=mom[0];
    output[1]=mom[1];
    output[2]=mom[2];
    output[3]=particle->getEnergy();
    return output;
  }
  
  //Returns the sum of the two; assumes a 4 vector
  static double* addVector(double* a, double* b, const int SIZE=4){
    double* output=new double[SIZE];
    for(int i=0; i<SIZE; ++i) output[i]=a[i]+b[i];
    return output;
  }

  //Returns transverse momentum magnitude
  static double getTMag(double* input){
    return sqrt(pow(input[0], 2) + pow(input[1], 2));
  }
 private:

};

#endif
