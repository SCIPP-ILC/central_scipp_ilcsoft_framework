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

#include "PopulationChecker.h"
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


PopulationChecker PopulationChecker;

static TFile* _rootfile;


PopulationChecker::PopulationChecker() : Processor("PopulationChecker") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}



void PopulationChecker::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile(_root_file_name.c_str(),"RECREATE");

    // usually a good idea to
    //printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;

}



void PopulationChecker::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void PopulationChecker::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    bool pion = false;
    cout << "*********************************** " << _nEvt << " ***************************************" << endl;
    LCCollection* col = evt->getCollection( _colName ) ;

    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;

        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );


           int id = hit->getPDG();
           int stat = hit->getGeneratorStatus();

           if(stat==1){
                cout << id << endl;    
                if(abs(id)==211){pion = true;}
           }
        } 
    }
    if(pion = false){
        pion_event++;
        pions += to_string(_nEvt);
        pions += " ";
    }

    _nEvt ++ ;
}



void PopulationChecker::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void PopulationChecker::end(){ 
    cout << pions << endl;
    _rootfile->Write();
}
