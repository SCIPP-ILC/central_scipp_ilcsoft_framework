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

#include "CesarTest.h"
#include "scipp_ilc_utilities.h"
#include "scipp_ilc_globals.h"
#include "polar_coords.h"
#include <iostream>
#include <cmath>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/MCParticle.h>

#include <TFile.h>
#include <TH2D.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"



using namespace lcio;
using namespace marlin;
using namespace std;


CesarTest CesarTest;

static TFile* _rootfile;
static TH2F* _hitmap_beamcal;
static TH2F* _hitmap_outgoing;
static TH2F* _hitmap_incoming;
static TH2F* _hitmap_lumical;


CesarTest::CesarTest() : Processor("CesarTest") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
}



void CesarTest::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("hittest.root","RECREATE");
    _hitmap_beamcal = new TH2F("hitmap_beamcal","Hit Distribution",300.0,-150.0,150.0,300.0,-150.0,150.0);
    _hitmap_outgoing = new TH2F("hitmap_outgoing","Hit Distribution",300.0,-150.0,150.0,300.0,-150.0,150.0);
    _hitmap_incoming = new TH2F("hitmap_incoming","Hit Distribution",300.0,-150.0,150.0,300.0,-150.0,150.0);
    _hitmap_lumical = new TH2F("hitmap_lumical","Hit Distribution",300.0,-150.0,150.0,300.0,-150.0,150.0);

    // usually a good idea to
    //printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;

}



void CesarTest::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void CesarTest::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    int id, stat, gen;
    const double* mom;
    double pos[] = {0, 0, 0};
    double in_energy, out_energy, out_x;
    
    LCCollection* col = evt->getCollection( _colName ) ;

    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;

        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
            
           id = hit->getPDG();
           gen = hit->getGeneratorStatus(); 
           
           mom = hit->getMomentum();
           
           if(gen==1){
                pos[2] = scipp_ilc::_BeamCal_zmin;
                pos[0] = mom[0]*pos[2]/mom[2];      
                pos[1] = mom[1]*pos[2]/mom[2];  

                in_energy = hit->getEnergy();

                scipp_ilc::transform_to_lab(mom[0], in_energy, out_x, out_energy);   

                pos[0] = out_x*pos[2]/mom[2];

                scipp_ilc::z_to_beam_out(pos[0], pos[1], pos[2]);

                stat = scipp_ilc::get_hitStatus(pos[0], pos[1]);
           }
           
           switch(stat){
      
         case 1:
            _hitmap_beamcal->Fill(pos[0], pos[1]);
            break;
         case 2:
            _hitmap_lumical->Fill(pos[0], pos[1]);
            break;
         case 3:
            _hitmap_outgoing->Fill(pos[0], pos[1]);
            break;
         case 4: 
            _hitmap_incoming->Fill(pos[0], pos[1]);
            break;
           }
             
        } 
    }

    _nEvt ++ ;
}



void CesarTest::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void CesarTest::end(){ 
    _rootfile->Write();
}
