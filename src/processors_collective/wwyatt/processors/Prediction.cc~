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
  _p_theta = new TH1F("p_theta", "Theta between positron and hadronic system", 360, 0, .1);
  _e_theta = new TH1F("e_theta", "Theta between positron and hadronic system", 360, 0, .1);
  _vector = new TH1F("vector", "Vector", 200, 0.0, 0.05);
  zmom=new TH1F("zmom", "Z-Momentum Distribution", 300, 0, 500);
  tmom=new TH1F("tmom", "T-Momentum Distribution", 300, 0, 10);

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
  getJane(col);

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

  //DEBUG plotting the hadronic transverse momentum sum.
  tmom->Fill(getTMag(data.hadronic_nopseudo));

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
  //fourvec nelectron = fourvec(-data.electronic.x,-data.electronic.y,data.electronic.z,data.electronic.e);
  //data.hadronic=nelectron;
  double mag=data.mag; //Magnitude of the hadronic vector?
  double electronTheta=getTheta(p.electron); //Angle off of z-axis
  double positronTheta=getTheta(p.positron);
  
  //DEBUG
  //cout << "will : " << hadronic.e << " : " << hadronic.z << " : " << p.electron.z << endl;
  //fourvec pred_e_lab = transform_to_lab(p.electron);
  //cout << "will : " << pred_e_lab.x << ":" << pred_e_lab.y << ": " << pred_e_lab.e << endl;
  //cout << "before : " << p.electron.e <<endl;
  //END DEBUG

  if(mag>1.0){
    _prediction->Fill(electronTheta,positronTheta);
  }

  double p_mag = getMag(electronic);
  double e_theta=getTheta(electron,p.electron);
  double p_theta=getTheta(positron,p.positron);
  //cout << "BE " << e_theta << endl;
  _p_theta->Fill(e_theta);
  _e_theta->Fill(p_theta);
  
  //Create position vectors
  fourvec real_e = getBeamcalPosition(data.electron,  1);
  fourvec real_p = getBeamcalPosition(data.positron, -1);
  fourvec pred_e = getBeamcalPosition(p.electron,  1);
  fourvec pred_p = getBeamcalPosition(p.positron, -1);

  
  //cout << "Will real e pos : " << getTMag(real_e) << endl;
  //cout << "Will pred e cmp : " << pred_e.x << "\t" << pred_e.y << endl;
  //cout << "Will pred e pos : " << getTMag(pred_e) << endl;
  
  //The following is a for a hit miss table to test efficiancy.
  //These booleans are true if the particle had hit.
  bool actual_electron=get_hitStatus(real_e)<3;
  bool actual_positron=get_hitStatus(real_p)<3;
  bool predicted_electron=get_hitStatus(pred_e)<3;
  bool predicted_positron=get_hitStatus(pred_p)<3;
  ++total;
  Will::META meta = Will::getMETA();
  
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

  map<int, string> key;
  key[0]="pred_e";key[1]="real_e";key[2]="pred_p";key[3]="real_p";
  vector<fourvec> jane = {meta.pred_e, meta.real_e, meta.pred_p, meta.real_p};
  vector<fourvec> will = {pred_e, real_e, pred_p, real_p};
  /*for(unsigned short int i=0; i < jane.size(); ++i){
    if(jane[i].x-will[i].x > 0.0 || jane[i].y-will[i].y > 0.0 || jane[i].z-will[i].z > 0.0){
      cout <<endl << "=== "<< key[i] <<" ===" << endl;
      cout << "J " << jane[i].x << " : " << jane[i].y << " : " << jane[i].z << endl;
      cout << "W " << will[i].x << " : " << will[i].y << " : " << will[i].z << endl;
    }
    }*/

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
 
  cout << "TOTAL: " << meta.total << endl; 
  double hh2=(double)meta.hh/(double)meta.total;
  double hm2=(double)meta.hm/(double)meta.total;
  double mh2=(double)meta.mh/(double)meta.total;
  double mm2=(double)meta.mm/(double)meta.total;

  cout << "HH: " << hh2 << endl;
  cout << "HM: " << hm2 << endl;
  cout << "MH: " << mh2 << endl;
  cout << "MM: " << mm2 << endl;
}
