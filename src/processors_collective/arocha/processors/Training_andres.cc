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

#include "Training_andres.h"
#include "scipp_ilc_utilities.h"
#include <iostream>
#include <math.h>

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


Training_andres Training_andres;

static TFile* _rootfile;
static TH2F* _Energy_Teta;
static TH1F* _Energy;


Training_andres::Training_andres() : Processor("Training_andres") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}



void Training_andres::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile(_root_file_name.c_str(),"RECREATE");
    _Energy_Teta = new TH2F("EvsT","Energy vs. Theta",300.0,0.0,0.1,300.0,0.0,250.0);
    _Energy = new TH1F("Eng","Energy Distribution",300.0,0.0,250.0);

    // usually a good idea to
    //printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;

}



void Training_andres::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void Training_andres::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    LCCollection* col = evt->getCollection( _colName ) ;

    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;

        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );


           double mom[3];
	   mom[0] = hit->getMomentum()[0];
	   mom[1] = hit->getMomentum()[1];
	   mom[2] = hit->getMomentum()[2];
	   
	   double momT = sqrt(pow(mom[0],2) + pow(mom[1],2));
	   double teta = atan(momT / mom[2]);


	   double energy = hit->getEnergy();
	   _Energy->Fill(energy);
           _Energy_Teta->Fill(teta,energy);
        } 
    }

    _nEvt ++ ;
}



void Training_andres::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void Training_andres::end(){ 
    _rootfile->Write();
}
