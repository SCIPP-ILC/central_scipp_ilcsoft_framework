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

#include "Bruce.h"
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

Bruce Bruce;

Bruce::Bruce() : Processor("Bruce") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
}

void Bruce::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _nRun = 0 ;
    _nEvt = 0 ;
}



void Bruce::processRunHeader( LCRunHeader* run) { 

} 

void Bruce::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...
    LCCollection* col = evt->getCollection( _colName ) ;
    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
           id = hit->getPDG();
           stat = hit->getGeneratorStatus();
           mom = hit->getMomentum();
           if(stat==1){

           }//end final state
        }//end for loop
    }
}

void Bruce::check( LCEvent * evt ) { 

}

void Bruce::end(){
    _rootfile->Write();
}


