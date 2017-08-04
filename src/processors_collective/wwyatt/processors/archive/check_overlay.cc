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

#include "check_overlay.h"
#include "scipp_ilc_utilities.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/MCParticle.h>

#include <TFile.h>
#include <TH2D.h>

#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <algorithm>
// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

using namespace lcio;
using namespace marlin;
using namespace std;

check_overlay check_overlay;

static TFile* _rootfile;

static TH1F* _energy;
static TH1F* _energy1;
static TH1F* _energy2;
static TH1F* _energy3;
static TH1F* _energy4;
static TH1F* _energy5;
static TH1F* _energy6;
static TH1F* _energy7;
static TH1F* _energy8;
static double totalEnergy=0;
static double[9] ergy = {0};
check_overlay::check_overlay() : Processor("check_overlay") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
    
    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}

void check_overlay::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;
    _rootfile = new TFile("check_overlay.root","RECREATE");
    _energy = new TH1F("energy", "Energy Distribution", 300,0,6000);
    _energy1 = new TH1F("energy1", "Energy Distribution 200-600", 300,200,600);
    _energy2 = new TH1F("energy2", "Energy Distribution 650-1050", 300,650,1050);
    _energy3 = new TH1F("energy3", "Energy Distribution 1050-1600", 300,1050,1600);
    _energy4 = new TH1F("energy4", "Energy Distribution 1600-2100", 300,1600,2100);
    _energy5 = new TH1F("energy5", "Energy Distribution 2100-2600", 300,2100,2600);
    _energy6 = new TH1F("energy6", "Energy Distribution 2600-3200", 300,2600,3200);
    _energy7 = new TH1F("energy7", "Energy Distribution 3200-3600", 300,3200,3600);
    _energy8 = new TH1F("energy8", "Energy Distribution 3800-4200", 300,3800,4200);
    _nEvt = 0 ;
}

void check_overlay::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 


//Bruce - Add energy of fianl state particles to verify if the particles are the same.
void check_overlay::processEvent( LCEvent * evt ) { 
    LCCollection* col = evt->getCollection( _colName );
    int stat, id =0;
    double tEnergy = 0;
    bool v=false;
    if( col != NULL ){
      //    printAllEvents(col);
        int nElements = col->getNumberOfElements();
	for(int i = 0; i < col->getNumberOfElements(); ++i){
	  MCParticle* hit = dynamic_cast<MCParticle*>(col->getElementAt(i));
	  if(hit->getGeneratorStatus()==1)tEnergy+=hit->getEnergy();
	}
	if(v && tEnergy<=10.0){
	  parser::printAllEvents(col);
	}
	_energy->Fill(tEnergy);
	_energy1->Fill(tEnergy);
	_energy2->Fill(tEnergy);
	_energy3->Fill(tEnergy);
	_energy4->Fill(tEnergy);
	_energy5->Fill(tEnergy);
	_energy6->Fill(tEnergy);
	_energy7->Fill(tEnergy);
	_energy8->Fill(tEnergy);
    }    
}

void check_overlay::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void check_overlay::end(){ 
  cout << "Total energy of all events: " << totalEnergy << endl;
  _rootfile->Write();
}
