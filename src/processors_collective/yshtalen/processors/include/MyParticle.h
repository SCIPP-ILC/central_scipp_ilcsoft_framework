//created by Jane Shtalenkova
//January 10, 2016
//extends LCIO MCParticle class
//

#ifndef EVENT_MYPARTICLE_H
#define EVENT_MYPARTICLE_H 1

//#include "/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-09/lcio/v02-07/src/cpp/include/lcio.h"
#include <string>
#include <vector>
#include <iostream>
#include <cmath>

#include "scipp_ilc_utilities.h"
#include "scipp_ilc_globals.h"
#include "polar_coords.h"

#include "lcio.h"
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>

using namespace lcio;
using namespace std;

class MyParticle;

typedef std::vector<MyParticle*> MyParticleVec;
class MyParticle{
    
    MCParticle* source;
    double* momentum; 
    double en; 
    bool had, elec, detect, able = false;

public: 

    MyParticle();

    ~MyParticle();
    
    void setSource(MCParticle* in_part); 

    double* mom();
    
    double energy();

    int id();

    int stat(); 
    
    bool isHadronic();

    bool isElectronic();

    bool isDetectable();
    
    bool isDetected();

    void setHadronic(bool set);
    
    void setElectronic(bool set);

    void setDetectable(bool set);
    
    void setDetected(bool set);
    
    double* getLorentzMom();
    
}; //class    

#endif 
