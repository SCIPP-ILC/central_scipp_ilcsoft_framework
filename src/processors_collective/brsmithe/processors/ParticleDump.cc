#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0
/* 
 * Ok, so I like C++11. Unfortunately,
 * Marlin is built with ansi C, so the processor
 * constructor freaks out about the string that is
 * passed to it as an argument. The above two lines
 * fix that issue, allowing our code to be compatible
 * with ansi C class declarations.
 * Big thanks to Daniel Bittman for helping me fix this.
 */

/*
 * author Jane Shtalenkovae
 * August 5, 2016
 */

#include "ParticleDump.h"
#include "scipp_ilc_utilities.h"
#include "scipp_ilc_globals.h"
#include "polar_coords.h"
#include <iostream>
#include <cmath>

#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>

#include <TFile.h>
#include <TH2D.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"



using namespace lcio;
using namespace marlin;
using namespace std;

ParticleDump ParticleDump;

static TFile* _rootfile;
//static TH2F* _hitmap;
static TH1F* _mass;
//static TH1F* _scalar;
static TH1F* _vector;
//static TH2F* _hitmap
static TH1F* _endpoints;
static TH2F* _hitmap;

static TH1F* _xSum;
static TH1F* _ySum;

ParticleDump::ParticleDump() : Processor("ParticleDump") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
}



void ParticleDump::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("ParticleDump.root","RECREATE");
    _vector = new TH1F("vector", "Deflected Particle Momentum Magnitude, sqrt(pX^2+pY^2)", 2000.0, 0.0, 20.0);
    _mass = new TH1F("mass", "Deflected Particle sqrt(Q^2) = sqrt(E^2 - <del_p>^2)", 2000.0, 0.0, 3.0);
    _endpoints = new TH1F("endpoints", "Endpoint Distribution", 4000.0, -2000.0, 2000.0);
    _hitmap = new TH2F("hitmap", "Hit Distribution", 200.0, -10.0, 10.0, 200.0, -10.0, 10.0);

    // usually a good idea to
    //printParameters() ;

    _nRun = 0 ;
    _nEvt = 0;

}



void ParticleDump::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 


void ParticleDump::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...


    LCCollection* col = evt->getCollection( _colName ) ;

    
    double scatter_vec[] = {0, 0, 0};
    double mag = 0;
    double energy = 0;
    double theta;
    int id, stat;
    
    const double* mom;
    
    const bool onlyleptons = true;

    bool leftDetector;
    
    double Zero = 0.000000001;

    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
        
        //first, find last electron and positron in the event
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
    
           id = hit->getPDG(); 
           stat = hit->getGeneratorStatus();
   
	   


           if(stat==1){
	     if( id==11 || id ==-11 ){
	       const double * mom = hit->getMomentum();
	       leftDetector = hit->hasLeftDetector();
	       
	       if(leftDetector){}

	     }
	     
	     
	     
           }//end final state
	   
	   
        }//end for loop
        
    }//end collection
    _nEvt ++ ;
}//end process



void ParticleDump::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void ParticleDump::end(){
    _rootfile->Write();
}

