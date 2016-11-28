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

#include "first.h"
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


first first;

static TFile* _rootfile;
static TH2F* _hitmap;
static TH1F* _vector;

first::first() : Processor("first") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}



void first::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;
    cout << "Initialized" << endl;
    _rootfile = new TFile("eBpB_vector.root","RECREATE");
    // usually a good idea to
    //printParameters() ;
    _nEvt = 0 ;

    //Setting up Plots\\
    _vector = new TH1F("test_1", "This is a test.", 200.0, 0.0, 2000.0);
}



void first::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void first::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...
    cout << "Event started." << endl;
    LCCollection* col = evt->getCollection( _colName );

    int stat, id =0;
    double tot_mom[]={0, 0};
    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;

        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
        
            id = hit->getPDG();
            stat = hit->getGeneratorStatus();
            if(stat==1){
                cout << "Particle " << hitIndex << " with ID: " << id << endl;
		
		//Setting up Test Plots\\
		double thisEnergy = hit->getEnergy();
		double thisMomentum = *hit->getMomentum();
		_vector->Fill(id,hit->getEnergy());

            }//end final state   
        }//end for
         
    }

    _nEvt ++ ;
}



void first::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void first::end(){ 
    _rootfile->Write();
}
