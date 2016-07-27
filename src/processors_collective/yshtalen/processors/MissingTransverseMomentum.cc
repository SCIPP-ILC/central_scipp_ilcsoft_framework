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
 * author Christopher Milke
 * April 5, 2016
 */

#include "MissingTransverseMomentum.h"
#include "scipp_ilc_utilities.h"
#include <iostream>
#include <cmath>

#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/MCParticle.h>

#include <TFile.h>
#include <TH2D.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"



using namespace lcio;
using namespace marlin;
using namespace std;


MissingTransverseMomentum MissingTransverseMomentum;

static TFile* _rootfile;
static TH2F* _hitmap;

MissingTransverseMomentum::MissingTransverseMomentum() : Processor("MissingTransverseMomentum") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
}



void MissingTransverseMomentum::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("hitmapeBpW.root","RECREATE");
    _hitmap = new TH2F("hitmap","Hit Distribution",300.0,-150.0,150.0,300.0,-150.0,150.0);
    // usually a good idea to
    //printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;

    _no_def_count = 0;
    _e_def_count = 0;
    _p_def_count = 0;
}



void MissingTransverseMomentum::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void MissingTransverseMomentum::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    LCCollection* col = evt->getCollection( _colName ) ;
    cout << endl;
    cout << endl;
    cout << endl;
    cout << "event = " << _nEvt << endl;
    

    MCParticle* high_e;
    MCParticle* high_p;


    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
        
        
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
    
           int id = hit->getPDG(); 
           int stat = hit->getGeneratorStatus();
           
           if(stat==1){
                if(id==11){
                    high_e = hit;
                }
                if(id==-11){
                    high_p = hit;
                }
           }//end final state
        }//end for loop
        
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
    
           int id = hit->getPDG(); 
           int stat = hit->getGeneratorStatus();
           
           if(stat==1){

                //find high energy electron
                if(id==11){
                    if(hit->getEnergy()>high_e->getEnergy()){
                        high_e = hit;
                    }               
                }    
                //find high energy positron
                if(id==-11){
                    if(hit->getEnergy()>high_p->getEnergy()){
                        high_p = hit;
                    }               
                }    
           }//end final state
        }//end for loop

    _nEvt ++ ;
}



void MissingTransverseMomentum::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void MissingTransverseMomentum::end(){ 

    _rootfile->Write();
    cout << "Ratio of events with no deflection: " << _no_def_count/_nEvt << endl;
    cout << "Electron deflected: " << _e_def_count/_nEvt << endl;
    cout << "Positron deflected: " << _p_def_count/_nEvt << endl;
}
