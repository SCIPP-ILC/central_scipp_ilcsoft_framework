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
#include <iomanip> 

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
static TH1F* zmom;
static TH1F* tmom;


static vector<double> spread_e;
static vector<double> spread_p;

static int p_scatter=0;
static int e_scatter=0;
static int ph_th = 0;
static int ph_tm = 0;
static int pm_th = 0;
static int pm_tm = 0;
static int total = 0;

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
  //usually a good idea to
  //printParameters() ;
  _prediction = new TH2F("predict", "Predicted Angle of Scatter, Correct vs Incorrect Kinematics", 1000, 0.0, 0.01, 1000, 0.0, 0.01);
  _p_theta = new TH1F("p_theta", "Theta between positron and hadronic system", 360, 0, 3.5);
  _e_theta = new TH1F("e_theta", "Theta between positron and hadronic system", 360, 0, 3.5);
  _vector = new TH1F("vector", "Vector", 200, 0.0, 0.05);
  zmom=new TH1F("zmom", "Z-Momentum Distribution", 300, 0, 500);
  tmom=new TH1F("tmom", "T-Momentum Distribution", 300, 0, 300);

  _nEvt = 0 ;
}



void Prediction::processRunHeader( LCRunHeader* run) { 
  //    _nRun++ ;
} 

void Prediction::processEvent( LCEvent * evt ) { 

  // DISCLAIMER: THERE IS MORE DOCUMENTATION IN THE HEADER FILE (Will.h).
  LCCollection* col = evt->getCollection( _colName ) ;
  if( col == NULL )return;
  //Run janes code for comparison.
  Will::measure data = Will::getMeasure(col);  
  /* "getMeasure(LCCollection)" Takes the collection, calculates the hadronic system.
   * A measure is a struct with several fourvectors in it:
   * - electron - Highest energy final state electron
   * - positron - Highest energy final state positron
   * - hadronic - Calculated from all non final state high energy electron/positron
   * - electronic
   *   |=> Either the electron or the positron above WITH transverse momentum (e.g. scatter)
   *   |=> Transverse momentum means this is the particle that scattered.
   */
  if(!data.scattered || data.mag<=1) return;//If there was no scatter, then there is nothing to see.
  prediction p(data); //Store prediction vector as variable 'p';
  /* "prediction" will calculate and return two prediction vectors.
   * - electron 
   * - positron
   * Only one of these will be correct; that can be checked by measure.p_scatter and measure.e_scatter.
   * I will make it easier to get the correct prediction.
   */

  
  //I rename the particles just for ease of use. I might take this out later for explicitness.
  fourvec electron=data.electron,
    positron=data.positron,
    hadronic=data.hadronic,
    electronic=data.electronic;
  double mag=data.mag; //Magnitude of the hadronic vector?
  double electronTheta=getTheta(p.electron); //Angle off of z-axis
  double positronTheta=getTheta(p.positron);
  
  //DEBUG
  //Find the mag difference between the electron & -hadron vector.
  //cout << "electronic:hadronic " << getTMag(electronic) << ":" << getTMag(hadronic) << endl;
  //cout << "delta: " << abs(getTMag(electronic)-getTMag(hadronic)) << endl;
  //END DEBUG

  if(mag>1.0){
    _prediction->Fill(electronTheta,positronTheta);
  }

  double p_mag = getMag(electronic);
  double e_theta=getTheta(electronic,p.electron);
  double p_theta=getTheta(electronic,p.positron);

  _p_theta->Fill(e_theta);
  _e_theta->Fill(p_theta);
  
  //Create position vectors
  fourvec real_e = getBeamcalPosition(data.electron);
  fourvec real_p = getBeamcalPosition(data.positron);
  fourvec pred_e = getBeamcalPosition(p.electron);
  fourvec pred_p = getBeamcalPosition(p.positron);

  //The following is a for a hit miss table to test efficiancy.
  //These booleans are true if the particle had hit.
  bool actual_electron=get_hitStatus(real_e)<3;
  bool actual_positron=get_hitStatus(real_p)<3;
  bool predicted_electron=get_hitStatus(pred_e)<3;
  bool predicted_positron=get_hitStatus(pred_p)<3;
  ++total;

  if(data.p_scatter){
    if(actual_positron&&predicted_positron)ph_th++;
    else if( actual_positron && !predicted_positron)pm_th++;
    else if(!actual_positron &&  predicted_positron)ph_tm++;
    else if(!actual_positron && !predicted_positron)pm_tm++;
    p_scatter++;
  }else if(data.e_scatter){
    if(actual_electron&&predicted_electron)ph_th++;
    else if( actual_electron && !predicted_electron)pm_th++;
    else if(!actual_electron &&  predicted_electron)ph_tm++;
    else if(!actual_electron && !predicted_electron)pm_tm++;	
    e_scatter++;
  }

}

void Prediction::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void Prediction::end(){ 
  Will::META meta = Will::getMETA();
  cout << interest << endl;
  _rootfile->Write();
  cout << "p_scatter:e_scatter = " << p_scatter << ":"<<e_scatter<<endl;
  cout << "Scat:noScat " << meta.SCATTERS << ":" << meta.NOSCATTERS << endl;
  cout << "MSC : " << meta.MSC << endl;
  int sum = ph_th+ph_tm+pm_th+pm_tm;
  double hh=1.0*ph_th/sum;
  double hm=1.0*ph_tm/sum;
  double mh=1.0*pm_th/sum;
  double mm=1.0*pm_tm/sum;
  
  cout << "          | Truth Hit | Truth Miss" << endl;
  cout << setprecision(4) << "Pred Hit  | " << hh << "   |  " << hm << endl;
  cout << setprecision(4) << "Pred Miss | " << mh << "   |  " << mm << endl;

}
