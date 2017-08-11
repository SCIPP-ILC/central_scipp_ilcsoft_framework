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

#include "SVM_Test.h"
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


SVM_Test SVM_Test;

static TFile* _rootfile;
static TH1D* _S;
static TH1D* _V;

SVM_Test::SVM_Test() : Processor("SVM_Test") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}



void SVM_Test::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile(_root_file_name.c_str(),"RECREATE");
    _S = new TH1D("S","S Hit Distribution",200,0.0,20.0);
    _V = new TH1D("V","V Hit Distribution",200,0.0,20.0);


    // usually a good idea to
    //printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;

}



void SVM_Test::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void SVM_Test::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    LCCollection* col = evt->getCollection( _colName ) ;

    double S = 0.0;
    double V = 0.0;
    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements() ;
	double v_px = 0.0;
	double v_py = 0.0;
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
	    MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
	  if(1==(hit->getGeneratorStatus())){
	    const double* momentum = hit->getMomentum();
	    double s_px = momentum[0];
	    double s_py = momentum[1];
	    S += sqrt((s_px*s_px)+(s_py*s_py));
	    v_px+=s_px;
	    v_py+=s_py;
	  }
        } 
	V=sqrt((v_px*v_px)+(v_py*v_py));


    }
	_S->Fill(double_t(S));
	_V->Fill(double_t(V));
    _nEvt ++ ;
}



void SVM_Test::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void SVM_Test::end(){ 
    _rootfile->Write();
}
