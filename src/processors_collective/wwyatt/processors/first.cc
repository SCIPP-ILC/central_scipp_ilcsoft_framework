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

//Header with Bhabah Processing Code
#include "Bhabha_util.h"


using namespace lcio;
using namespace marlin;
using namespace std;


first first;

static TFile* _rootfile;
static TH2F* _hitmap;
static TH1F* _energy;
static TH1F* _phi;
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
    _rootfile = new TFile("will_test.root","RECREATE");
    _energy = new TH1F("energy", "Energy", 520.0,  0.0, 260.0);
    _hitmap = new TH2F("pos", "Position Distribution", 300.0, -150.0, 150.0, 300.0, -150.0, 150.0);
    _phi = new TH1F("phi", "Angle Difference", 300, .1, 1);

    //printParameters() ;
    _nEvt = 0 ;

    //Setting up Plots\\
    cout << "Vector allocated." << endl;
}



void first::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 

class BhabahBundle:BB::Bundle{
public:
  //Where I do my math to output to graph.
  void initialize(){
    double p_phi = getPhi(positron);
    double e_phi = getPhi(electron);
    double del_phi = p_phi-e_phi;

    //Put in graph
    _phi->Fill(del_phi);
  }

private:
}

void first::processEvent( LCEvent * evt ) { 
    LCCollection* col = evt->getCollection( _colName );
    
    BhabahBundle* bundle = new BhabahBundle; //Created a bundle to process the data for the bahbah event.

    int stat, id =0;
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );        
            id = hit->getPDG();
            stat = hit->getGeneratorStatus();
            if(stat==1){
	      //hit is an end particle:
	      bundle->addParticle(hit); //I add it to a Bahbah class to do the rest of the work.

            }//end final state   
        }//end for
    }
    delete[] bundle;
    _nEvt ++ ;
}



void first::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void first::end(){ 
    _rootfile->Write();
}
