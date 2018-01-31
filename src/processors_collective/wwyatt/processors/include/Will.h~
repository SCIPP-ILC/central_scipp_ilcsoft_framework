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
#include <sstream>

#include <TFile.h>
#include <TH2D.h>

#include "lcio.h"
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include "scipp_ilc_utilities.h"

using namespace lcio;
using namespace std;
namespace Will{

  /* ==== fourvec ====
   * I made fourvec to make all the code more readable.
   * To make a fourvec just pass it components
   * fourvec example(momentum[0], momentum[1], momentum[2], particle.getEnergy()); 
   * To get momentum:
   * example.x //return x momentum.
   * example.T //returns transverse momentum
   * Also it is case insensitive so:
   * example.y = example.Y
   */
  struct fourvec{
    union{ double X; double x=0.0; };
    union{ double Y; double y=0.0; };
    union{ double Z; double z=0.0; };
    union{ double E; double e=0.0; };
    union{ double T; double t=0.0; };
    fourvec operator+(const fourvec& a) const;
    double operator*(const fourvec& a) const;
    fourvec operator*(const double) const;
    fourvec operator+=(const fourvec& a);
    fourvec operator*=(const double a);
    fourvec operator*=(const fourvec& a);
    fourvec operator/=(const double a);
    fourvec operator/(const double a)const;
    fourvec();
    fourvec(const double,const double);
    fourvec(const double,const double,const double);
    fourvec(const double,const double,const double,const double);
    fourvec(const double,const double,const double,const double,const double);
    fourvec(const double*,const unsigned short SIZE);
  };
  
  //Used for analyzing efficency
  //This object makes analyzing energy cuts easy.
  struct Bundle{
    fourvec predicted;
    fourvec actual;
    double system_energy;
  };

  //This is used to store all the vectors needed in two-photon analysis.
  struct measure{
    fourvec hadronic_nopseudo;
    fourvec hadronic;
    fourvec electronic;
    fourvec electron;
    fourvec positron;
    fourvec pseudo;
    
    double mag=0.0;
    bool scattered=false;
    bool p_scatter=false;
    bool e_scatter=false;
  };

  //This will contain and calculate a prediction vector from a filled measure stucture.
  struct prediction{
    fourvec electron;
    fourvec positron;
    double alpha=0;
    double beta=0;
    prediction(double x,double y);
    prediction( measure );
  };
  
  //Used to store data about a HM Grid
  //This is an analysis tool.
  struct hmgrid{
    int hh=0;
    int hm=0;
    int mh=0;
    int mm=0;
  };

  //Keep track of numbers for debugging.
  struct META{
    int SCATTERS=0;
    int NOSCATTERS=0;
    int STATS[5]={0};
    int MSC=0;
    fourvec tmp;
    fourvec pred_e, pred_p, real_e, real_p;
    double fv[4];
    fourvec* _tmp;
    int hh=0,hm=0,mh=0,mm=0,total=0;
    int err_direction=0;
    static const int BEAMCAL = 3265; //Distance to beamcal
  };
  static META meta;

  META getMETA();

  //Specific function used in prediction algorithm.
  //Finds the highest energy particle
   map<int,MCParticle*> maxParticle(LCCollection*, initializer_list<int>);
   map<int,double> maxEnergy(LCCollection*, 
			     initializer_list<int> ids, 
			     vector<MCParticle*>& final_state);
   
   //Returns string in x:y:z format
   //string str(fourvec in);

   //Gets momentum vector as non constant
   double* getVector(MCParticle*);

   //casting function for MCParticle to fourvec
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
   int get_hitStatus(const fourvec, const bool=true);
   int get_hitStatus(MCParticle*);
   
   //Like the ilc version but it supports MCParticle and fourvectors. 
   //Also it returns a new foucvec that has been transformed.
   fourvec transform_to_lab(MCParticle*);
   fourvec transform_to_lab(const fourvec);

   //Used to debug old code.
   double* legacy(fourvec);


  /*Returns a map with a few four vectors in it.
   * - hadronic vector
   * - electronic vector
   * - electron vector
   * - positron vector
   * This should be used to calculate a prediction vector.
   */
   measure getMeasure(LCCollection*);
   
   //Returns a position fourvec, of the particle on the face of the beamcal.
   fourvec getBeamcalPosition(fourvec, signed short = 0);

   //Calculates a HM Grid and stores it in a hmgrid object.
   hmgrid getHMGrid(vector<fourvec> predicted, vector<fourvec> actual);
   void printHMGrid(vector<fourvec> predicted, vector<fourvec> actual);
   void printHMGrid(hmgrid);

   //Calculates a HM grid with the option of an energy cut.
   hmgrid getHMGrid(vector<Bundle> input, double energy_cut=0.0);   
   void printHMGrid(vector<Bundle> input, double energy_cut=0.0);
   
   //Helper function when calculating the HM Grid,
   //this is the code that checks to see if the particle hit the beamcal.
   void recordHMValue(hmgrid &output, fourvec predicted, fourvec actual);

   //Gets a graph of the angle distribution.
   TH1F* getDistribution(string name,vector<Bundle> input, double energy_cut=0.0, double upper_bound=0.006);

   //Python version of cout. I don't use these.
   void print(string );
   void print(string , string );
}
#endif
