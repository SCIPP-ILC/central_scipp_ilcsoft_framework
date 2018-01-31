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
#include <TwoPhoton.h>

#include <TFile.h>
#include <TH2D.h>

#include "lcio.h"
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include "scipp_ilc_utilities.h"

using namespace lcio;
using namespace std;
namespace Will{

  

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
   
   //Returns string in x:y:z format
   //string str(fourvec in);


   //Used to debug old code.
   double* legacy(fourvec);


   //Gets a graph of the angle distribution.
   TH1F* getDistribution(string name,vector<TwoPhoton::Result> input, double energy_cut=0.0, double upper_bound=0.006);

   //Python version of cout. I don't use these.
   void print(string );
   void print(string , string );
}
#endif
