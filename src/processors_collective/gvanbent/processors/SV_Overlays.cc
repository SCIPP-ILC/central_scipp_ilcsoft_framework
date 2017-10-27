#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0
/* 
 * C. Milke: 
 * Ok, so I like C++11. Unfortunately,
 * Marlin is built with ansi C, so the processor
 * constructor freaks out about the string that is
 * passed to it as an argument. The above two lines
 * fix that issue, allowing our code to be compatible
 * with ansi C class declarations.
 * Big thanks to Daniel Bittman for helping me fix this.
 */

/*
 * author Yevgeniya Shtalenkova
 * June 6, 2017
 */

#include "SV_Overlays.h"
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

#define _USE_MATH_DEFINES


using namespace lcio;
using namespace marlin;
using namespace std;
using namespace scipp_ilc;

SV_Overlays SV_Overlays;

static TFile* _rootfile;
static TH1D* _S;
static TH1D* _V;

SV_Overlays::SV_Overlays() : Processor("SV_Overlays") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}



void SV_Overlays::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("SV_test.root","RECREATE");
    _S = new TH1D("S", "Total Scalar Magnitude", 200, 0.0, 20.0);
    _V = new TH1D("V", "Total Vector Magnitude", 200, 0.0, 20.0);
    
    // usually a good idea to
    //printParameters() ;
    _nEvt = 0 ;

    _tot = 0;
}



void SV_Overlays::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void SV_Overlays::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...
    LCCollection* col = evt->getCollection( _colName ) ;

    int stat, id = 0;
    //total scalar magnitude of hadronic system
    double S = 0;
    //total vector magnitude of system
    double tot_mom[] = {0, 0};
    double V;

    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
        cout << endl;
        cout << endl;
        cout << "***************************EVENT: " << _nEvt << "****************************" << endl;

        vector<MCParticle*> final_system;
            //create final state subsystem
            //determine beam particle energies for identification
            for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){

                MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
                stat = hit->getGeneratorStatus();
                if(stat==1){
                    double mom[2];
                    mom[0] = hit->getMomentum()[0];
                    mom[1] = hit->getMomentum()[1];

                    double mag = sqrt(pow(mom[0], 2)+pow(mom[1], 2));
                    S+=mag;

                    tot_mom[0]+=mom[0];
                    tot_mom[1]+=mom[1]; 
                }
            }
            V = sqrt(pow(tot_mom[0], 2)+pow(tot_mom[1], 2)); 
            _S->Fill(S);
            _V->Fill(V);
    }//end collection

    _nEvt ++ ;
}



void SV_Overlays::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void SV_Overlays::end(){ 
    _rootfile->Write();
}
