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

#include "example.h"
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


example example;

static TFile* _rootfile;
static TH2F* _plot;
static TH1F* _histo;

static int _nEvt=0;


example::example() : Processor("example") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
    
    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}


void example::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;
    cout << "Initialized " << endl;

    _rootfile = new TFile("example.root","RECREATE");

    _plot = new TH2F("hh", "Hit-Hit HeatMap", 300.0, -150.0, 150.0, 300.0, -150.0, 150.0);
    _histo = new TH1F("energy", "Energy Distribution", 300, 0, 550);

    _nEvt = 0 ;
}



void example::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 
void example::processEvent( LCEvent * evt ) { 
    LCCollection* col = evt->getCollection( _colName );
    _nEvt++;
    int end_particles=0;
    cout << " +++ New Event +++" << endl;
    for(int i=0; i < col->getNumberOfElements(); ++i){
      MCParticle* particle=dynamic_cast<MCParticle*>(col->getElementAt(i));
      int pid=particle->getPDG();
      int state=particle->getGeneratorStatus();
      if(state==1) ++end_particles;
    }
}


void example::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void example::end(){ 
  cout << endl << "Total # of events: " << _nEvt << endl;
  
  _rootfile->Write();
}
