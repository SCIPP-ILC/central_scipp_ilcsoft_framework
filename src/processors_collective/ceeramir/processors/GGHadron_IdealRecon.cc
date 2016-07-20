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

#include "GGHadron_IdealRecon.h"
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


GGHadron_IdealRecon GGHadron_IdealRecon;

static TFile* _rootfile;
static TH2F* _hitmap;
static TH1F* _pVec;


GGHadron_IdealRecon::GGHadron_IdealRecon() : Processor("GGHadron_IdealRecon") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
}



void GGHadron_IdealRecon::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("hitmap.root","RECREATE");
    _hitmap = new TH2F("hitmap","Hit Distribution",300.0,-150.0,150.0,300.0,-150.0,150.0);
    _pVec = new TH1F("pVec", "Scatter Transverse Momentum Magnitude", 1000.0, 0.0, 260.0); 

    // usually a good idea to
    //printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;

}



void GGHadron_IdealRecon::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void GGHadron_IdealRecon::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    LCCollection* col = evt->getCollection( _colName ) ;

    const double* comp;
    //const double* pTot;
    MCParticle* high;

    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
        
        MCParticle* first = dynamic_cast<MCParticle*>( col->getElementAt(0) );
        comp = first->getMomentum();
        high = first;
         
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
    
           const double* mom = hit->getMomentum();
           cout << mom << endl;   
           if (mom > comp){
                high = hit;
           }
        }
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );

           const double* pos = hit->getEndpoint();
           _hitmap->Fill(pos[0],pos[1]);
        
            /*if(hit!=high){
                pTot+= hit->getMomentum();    
                
            }*/
        }
        
    cout << "highest momentum: " << high->getMomentum() << endl;
    }
    //_pVec->Fill(pTot);
    _nEvt ++ ;
}



void GGHadron_IdealRecon::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void GGHadron_IdealRecon::end(){ 
    _rootfile->Write();
}
