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
#include "scipp_ilc_utilities.h"

using namespace lcio;
using namespace std;
namespace Will{
  struct fourvec{
    union{ double X; double x=0.0; };
    union{ double Y; double y=0.0; };
    union{ double Z; double z=0.0; };
    union{ double E; double e=0.0; };
    union{ double T; double t=0.0; };
    fourvec operator+(const fourvec& a) const;
    double operator*(const fourvec& a) const;
    fourvec operator+=(const fourvec& a);
    fourvec();
    fourvec(const double,const double);
    fourvec(const double,const double,const double);
    fourvec(const double,const double,const double,const double);
    fourvec(const double,const double,const double,const double,const double);
    fourvec(const double*,const unsigned short SIZE);
  };
  
  
  struct prediction{
    fourvec electron;
    fourvec positron;
    prediction(double x,double y);
  };
  
  
  struct measure{
    fourvec hadronic;
    fourvec electronic;
    fourvec electron;
    fourvec positron;
    double mag=0.0;
    bool scattered=false;
    bool p_scatter=false;
    bool e_scatter=false;
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
   double getTMag(const fourvec);
   double getTMag(const double*);
   
   //Returns momentum from a momentum vector
   double getMag(const double*);
   double getMag(const fourvec);

   //Returns angle off of the z-axis, theta
   double getTheta(const double*);
   double getTheta(const fourvec);

   //Returns dot product of two vectors
   double getDot(const double*, const double*);
   double getDot(const fourvec, const fourvec);

   //Retuns anglebetween vectors or doubles in rads
   double getTheta(const double*, const double*);
   double getTheta(const fourvec, const fourvec);

   //Returns hit status
   // 1 - hit Beamcal
   // 2 - outside Beamcal radius
   // 3 - outgoing beampipe hole
   // 4 - incoming beampipe hole
   int get_hitStatus(const fourvec);
   static int STATS[5]={0};

   double* legacy(fourvec);


  /*Returns a map with a few four vectors in it.
   * - hadronic vector
   * - electronic vector
   * - electron vector
   * - positron vector
   * This should be used to calculate a prediction vector.
   */
   measure getMeasure(LCCollection*);
   
   void print(string );
   void print(string , string );


}

#endif
