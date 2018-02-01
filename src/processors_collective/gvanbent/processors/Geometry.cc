
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

#include "Geometry.h"
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


Geometry Geometry;

static TFile* _rootfile;
static TH2F*  _hitmap;  
static TH1F* _energy;

Geometry::Geometry() : Processor("Geometry") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
}



void Geometry::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;
    cout << "I have run init" << endl ;
    _rootfile = new TFile("hitmap.root","RECREATE"); //output file
    _energy = new TH1F("energy", "Energy Distribution", 300.0, 0, 275);
    _hitmap = new TH2F("hitmap","Hit Distribution",300.0,-250.0,250.0,300.0,-250.0,250.0);  //creates scatterplot with name "hitmap", #of bins, xmin, xmax, # of y bins, ymin, ymax

    // usually a good idea to
    //printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;

}



void Geometry::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void Geometry::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...
    
    cout << "Event " << _nEvt++ << endl;
    LCCollection* col = evt->getCollection( _colName ) ;

    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;

        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
	  
	   MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) ); //col->getElementAt() same as  (*col).getElementAt()
           const double* pos = hit->getMomentum();
	   cout << "particl: " << hit->getPDG() << endl;
           _hitmap->Fill(pos[0],pos[1]);
	   _energy->Fill(hit->getEnergy());
        } 
    }


}



void Geometry::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void Geometry::end(){ 
    _rootfile->Write();
}
