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
using namespace TwoPhoton;


Prediction Prediction;


static TFile* _rootfile;
static TH2F* _prediction;
static TH2F* _observe;
static TH1F* _vector;
static TH1F* _p_theta;
static TH1F* _e_theta;
static TH1F* zmom;
static TH1F* tmom;
static TH1F* amom;
static TH1F* bmom;
static TH1F* cmom;
static TH1F* dmom;
static TH1F* _alpha;
static TH1F* _beta;
static TH1F* _az;
static TH1F* _ez;

static vector<Result> positron_results;
static vector<Result> electron_results;
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

static int photons=0;
static int nphotons=0;

static int n_events=0;
static int i_events=0;

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
  _observe = new TH2F("energy", "B vs System Energy", 1000, 0.0, 1.5, 1000, 0.0, 525);
  _p_theta = new TH1F("p_theta", "Theta between positron and hadronic system", 360, 0, .1);
  _e_theta = new TH1F("e_theta", "Theta between positron and hadronic system", 360, 0, .1);
  _vector = new TH1F("vector", "Vector", 200, 0.0, 0.05);
  zmom=new TH1F("zmom", "System energy", 500, 450, 550);
  tmom=new TH1F("tmom", "Eletron system energy", 500, 450, 550);
  amom=new TH1F("amom", "Distribution of Actual Positron Theta", 500, 0,3.15);
  bmom=new TH1F("bmom", "Distribution of Acutal Electron Theta", 500, 0,3.15);
  cmom=new TH1F("cmom", "Distribution of Predicted Positron Theta", 500, 0,3.15);
  dmom=new TH1F("dmom", "Distribution of Predicted Electron Theta", 500, 0,3.15);
  _alpha=new TH1F("alpha", "Alpha Distribution", 500, 0,550);
  _beta=new TH1F("beta", "Beta Distribution", 500, 0,550);
  _az=new TH1F("az", "Posistron Z-Momentum", 500, -250,250);
  _ez=new TH1F("ez", "Electron Z-Momentum", 500, -250,250);
  _nEvt = 0 ;
}



void Prediction::processRunHeader( LCRunHeader* run) { 
  //    _nRun++ ;
} 

void Prediction::processEvent( LCEvent * evt ) { 
  // DISCLAIMER: THERE IS MORE DOCUMENTATION IN THE HEADER FILE (Will.h).
  LCCollection* col = evt->getCollection( _colName ) ;
  if( col == NULL )return;
  n_events++;
  //Run janes code for comparison.
  //getJane(col);

  TwoPhoton::bundle data = TwoPhoton::getHadronicSystem(col);  
  /* "getHadronicSystem(LCCollection)" Takes the collection, calculates the hadronic system.
   * A bundle is a struct with several fourvectors in it:
   * - electron - Highest energy final state electron
   * - positron - Highest energy final state positron
   * - hadronic - Calculated from all non final state high energy electron/positron
   * - electronic
   *   |=> Either the electron or the positron above WITH transverse momentum (e.g. scatter)
   *   |=> Transverse momentum means this is the particle that scattered.
   */


  //If there was no scatter, then there is nothing to see.
  //Also excludes small magnitude events.
  if(data.mag<=1) return;
  //if(!data.scattered) return;
  i_events++;

  prediction p(data); //Store prediction vector as variable 'p';
  /* "prediction" will calculate and return two prediction vectors.
   * - electron 
   * - positron
   * Only one of these will be correct; that can be checked by bundle.p_scatter and bundle.e_scatter.
   * I will make it easier to get the correct prediction.
   */

  /*Debugging wb events straight don't even work.*/
  _alpha->Fill(p.alpha);_beta->Fill(p.beta);
  _az->Fill(p.positron.z);_ez->Fill(p.electron.z);

  //Collecting Data for Analysis
  Result positron_result;
  Result electron_result;
  double tot_energy=data.positron.e+data.electron.e+data.hadronic_nopseudo.e;
  positron_result.system_energy=tot_energy;
  electron_result.system_energy=tot_energy;

  positron_result.actual=data.positron;
  positron_result.predicted=p.positron;
  electron_result.actual=data.electron;
  electron_result.predicted=p.electron;
  //cout << getPhi(data.electron) << "  :  " << getPhi(data.positron) << endl;

  if(tot_energy > 494.0){
    /*    double phi_ap=getPhi(data.positron), phi_ae=getPhi(data.electron), phi_pp=getPhi(p.positron), phi_pe=getPhi(p.electron);
    if(phi_ap==phi_ap) amom->Fill(getPhi(data.positron));
    if(phi_ae==phi_ae) bmom->Fill(getPhi(data.electron));
    if(phi_pp==phi_pp) cmom->Fill(getPhi(p.positron));
    if(phi_pe==phi_pe) dmom->Fill(getPhi(p.electron));
    */
    amom->Fill(getTheta(data.positron));
    bmom->Fill(getTheta(data.electron));
    cmom->Fill(getTheta(p.positron));
    dmom->Fill(getTheta(p.electron));
  }
  positron_results.push_back(positron_result);
  electron_results.push_back(electron_result);
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
  //cout << "total photons: " << photons << endl;
  //cout << "Photons with pi-0 parent: " << nphotons << endl;
  cout << "==== " << n_events << " Events ====" << endl;
  cout << "==== " << i_events << " Counted Events (mag>1) ====" << endl;
  cout << "Predicted Z-Direction Errors:" << meta.err_direction << endl;
  //cout << "(scattered:not-scattered)\t= " << meta.SCATTERS << ":" << meta.NOSCATTERS << endl;  
  cout << "Electron Results Collected:" << electron_results.size() << endl;
  cout << "Positron Results Collected:" << positron_results.size() << endl;
  cout << endl;
  cout << "Electron HM Grid: " << endl;
  printHMGrid(electron_results);
  cout << endl;
  cout << "Positron  HM Grid: " << endl;
  printHMGrid(positron_results);

  double cut=494;
  /*  for(Result result: electron_results){
    tmom->Fill(result.system_energy);

    if(result.system_energy > cut)
      cmom->Fill(getTheta(result.actual, result.predicted));
  }
  for(Result result: positron_results){
    zmom->Fill(result.system_energy);
    if(result.system_energy > cut)
      dmom->Fill(getTheta(result.actual, result.predicted));
  }
  cout << acos(-1) << endl;
  */
  
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
  
  _rootfile->Write();
}
