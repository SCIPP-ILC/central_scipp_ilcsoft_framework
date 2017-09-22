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

#include "Prediction.h"
#include "scipp_ilc_utilities.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/MCParticle.h>

#include <TFile.h>
#include <TH2D.h>
#include <MyParticle.h>
#include <Will.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"



using namespace lcio;
using namespace marlin;
using namespace std;
using namespace Will;

Prediction Prediction;


static TFile* _rootfile;
static TH2F* _prediction;
static TH1F* _vector;
static TH1F* _p_theta;
static TH1F* _e_theta;

static vector<double> spread_e;
static vector<double> spread_p;



Prediction::Prediction() : Processor("Prediction") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}



void Prediction::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("ppredict.root","RECREATE");
    // usually a good idea to
    //printParameters() ;
    _prediction = new TH2F("predict", "Predicted Angle of Scatter, Correct vs Incorrect Kinematics", 1000, 0.0, 0.01, 1000, 0.0, 0.01);
    _p_theta = new TH1F("p_theta", "Theta between positron and hadronic system", 360, 0, 3.5);
    _e_theta = new TH1F("e_theta", "Theta between positron and hadronic system", 360, 0, 3.5);
    _vector = new TH1F("vector", "Vector", 200, 0.0, 0.05);
    _nEvt = 0 ;
}



void Prediction::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 

void Prediction::processEvent( LCEvent * evt ) { 
    LCCollection* col = evt->getCollection( _colName ) ;
    _nEvt++;
    if( col == NULL )return;
    vector<MCParticle*> hadronic_system;
    vector<MCParticle*> final_system;
    int stat, id =0;

    Will::measure data = Will::getMeasure(col);
    if(data.scattered == true){
      fourvec
	electron=data.electron,
	positron=data.positron,
	hadronic=data.hadronic,
	electronic=data.electronic;
      double mag=data.mag;
      double alpha=(500-hadronic.E -hadronic.z);
      double  beta=(500-hadronic.E +hadronic.z);

      prediction p(-hadronic.x,-hadronic.y);
      p.electron.z = -(pow(electron.T, 2)-pow(alpha, 2))/(2*alpha);
      p.positron.z = (pow(positron.T, 2)-pow(beta, 2))/(2*beta);

      double electronTheta=getTheta(p.electron);
      double positronTheta=getTheta(p.positron);

      if(mag>1.0){
	_prediction->Fill(electronTheta,positronTheta);
      }
      
      double p_mag = getMag(electronic);
      double e_theta=getTheta(electronic,p.electron);
      double p_theta=getTheta(electronic,p.positron);

      _p_theta->Fill(e_theta);
      _e_theta->Fill(p_theta);
    }
}

void Prediction::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void Prediction::end(){ 
    cout << interest << endl;
    _rootfile->Write();
}
