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

#include "EnergyFinder.h"
#include "scipp_ilc_utilities.h"
#include <iostream>

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


EnergyFinder EnergyFinder;

static TFile* _rootfile;
static TH2F* _hitmap;


EnergyFinder::EnergyFinder() : Processor("EnergyFinder") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}



void EnergyFinder::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile(_root_file_name.c_str(),"RECREATE");
    _hitmap = new TH2F("hitmap","Hit Distribution",300.0,-150.0,150.0,300.0,-150.0,150.0);

    // usually a good idea to
    //printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;

}



void EnergyFinder::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void EnergyFinder::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    LCCollection* col = evt->getCollection( _colName ) ;

    // this will only be entered if the collection is available
    if( col != NULL ){
        cout << "if statement entered"<< endl; 
        int nElements = col->getNumberOfElements()  ;
        
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           SimCalorimeterHit* hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt(hitIndex) );

           const float* pos = hit->getPosition();
           _hitmap->Fill(pos[0],pos[1]);
        } 
    }
    cout << "Finished Event " << _nEvt << endl; 
    _nEvt ++ ;
}



void EnergyFinder::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void EnergyFinder::end(){ 
    _rootfile->Write();
}
