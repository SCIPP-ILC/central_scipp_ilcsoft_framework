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
static TH1F* amom;
static TH1F* bmom;
static TH1F* cmom;

static vector<Bundle> results;
static vector<fourvec> actual;
static vector<fourvec> predicted;
static vector<double> spread_e;
static vector<double> spread_p;
static vector<pair<double,double>> spread;

static int p_scatter=0;
static int e_scatter=0;
static int ph_th = 0;
static int ph_tm = 0;
static int pm_th = 0;
static int pm_tm = 0;
static int total = 0;
static double acc = 0.0;
static int _acc=0;
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
  zmom=new TH1F("zmom", "Hadronic system energy", 500, 400, 530);
  tmom=new TH1F("tmom", "Theta Distribution", 500, 0, .006);  
  amom=new TH1F("amom", "Distribution of Theta Energy Above 0", 500, 0,.1);
  bmom=new TH1F("bmom", "Distribution of Theta Energy Above 480", 500, 0,.1);
  cmom=new TH1F("cmom", "Distribution of Theta Energy Above 493", 500, 0,.1);
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
  //getJane(col);

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


  //If there was no scatter, then there is nothing to see.
  //Also excludes small magnitude events.
  if(!data.scattered || data.mag<=1) return;
  //if(!data.scattered) return;
  
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

  Bundle bundle;
  bundle.system_energy=data.positron.e+data.electron.e; //No hadronic energy
  if(data.p_scatter){
    bundle.actual=positron;
    bundle.predicted=p.positron;
  }else if(data.e_scatter){
    bundle.actual=electron;
    bundle.predicted=p.electron;
  }
  results.push_back(bundle);
}

void Prediction::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void Prediction::end(){ 
  Will::META meta = Will::getMETA();
  /*
    
    cout << "Scatter Ratios:"<<endl;
    cout << "(positrons:electrons)\t= " << p_scatter << ":"<<e_scatter<<endl;
    cout << "Misc data: " << meta.MSC << endl;
  */
  //General Analysis
  cout << "Predicted Z-Direction Errors:" << meta.err_direction << endl;

  cout << "(scattered:not-scattered)\t= " << meta.SCATTERS << ":" << meta.NOSCATTERS << endl;  
  cout << "Bundles Collected:" << results.size() << endl;
  cout << endl;
  cout << "Energy Above " << 0 << endl;
  printHMGrid(results);
  cout << endl;
  cout << "Energy Above " << 480 << endl;
  printHMGrid(results, 480);
  cout << endl;
  cout << "Energy Above " << 493 << endl;
  printHMGrid(results, 493);
  for(Bundle bundle: results){
    if(bundle.system_energy > 0)
      amom->Fill(getTheta(bundle.actual, bundle.predicted));
    if(bundle.system_energy > 480)
      bmom->Fill(getTheta(bundle.actual, bundle.predicted));
    if(bundle.system_energy > 493)
      cmom->Fill(getTheta(bundle.actual, bundle.predicted));
  }

  //HM Grids
  //This for loop will find out when the algorithm becomes very accurate
  //The sweet spot is somewhere between 480 and 493.8 GeV
  /*
  int len=10;
  double start = 480;
  double end = 492;
  double step=(end-start)/len;
  for(int i=0; i < len; ++i){
    double current_cut=start+step*i;
    cout << "Energy above "<< current_cut << endl;
    printHMGrid(results, current_cut);
  }
  */

  //Get energy Distribution
  for(Bundle bundle:results){
    zmom->Fill(bundle.system_energy);
  }
  
  _rootfile->Write();
}
