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

static vector<Result> results;
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
  zmom=new TH1F("zmom", "Hadronic system energy", 500, 450, 550);
  tmom=new TH1F("tmom", "B=P/E for Hadronic System", 500, 0, 1.5);
  amom=new TH1F("amom", "Distribution of Theta Energy Above 0", 500, 0,.1);
  bmom=new TH1F("bmom", "Distribution of Theta Energy Above 480", 500, 0,.1);
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
  MCParticle* max_photon=NULL;
  for(int i=0; i < col->getNumberOfElements(); ++i){
    MCParticle* particle=dynamic_cast<MCParticle*>(col->getElementAt(i));
    int pid=particle->getPDG();
    if(pid==22){
      photons++;
      vector<MCParticle *> parents=particle->getParents();
      for(auto parent:parents){
	if(parent->getPDG()==111){
	  nphotons++;
	  continue;
	}
	if(max_photon==NULL || max_photon->getEnergy()<particle->getEnergy())
	  max_photon=particle;
      }
    }
  }
    

  prediction p(data); //Store prediction vector as variable 'p';
  /* "prediction" will calculate and return two prediction vectors.
   * - electron 
   * - positron
   * Only one of these will be correct; that can be checked by bundle.p_scatter and bundle.e_scatter.
   * I will make it easier to get the correct prediction.
   */

  double mag=data.mag; //Magnitude of the hadronic vector?
  double electronTheta=getTheta(p.electron); //Angle off of z-axis
  double positronTheta=getTheta(p.positron);

  Result result;
  result.system_energy=data.positron.e+data.electron.e+data.hadronic_nopseudo.e; //No hadronic energy
  if(data.p_scatter){
    result.actual=data.positron;
    result.predicted=p.positron;
  }else if(data.e_scatter){
    result.actual=data.electron;
    result.predicted=p.electron;
  }
  results.push_back(result);
  if(max_photon!=NULL){
    double tot_energy= data.electron.e+data.positron.e+data.hadronic.e;
    double b=getMag(data.hadronic_nopseudo)/data.hadronic_nopseudo.e;
      tmom->Fill(b);

    _observe->Fill(b, tot_energy);
    zmom->Fill(tot_energy);

  }
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
  cout << "total photons: " << photons << endl;
  cout << "Photons with pi-0 parent: " << nphotons << endl;
  cout << "Predicted Z-Direction Errors:" << meta.err_direction << endl;

  cout << "(scattered:not-scattered)\t= " << meta.SCATTERS << ":" << meta.NOSCATTERS << endl;  
  cout << "Results Collected:" << results.size() << endl;
  cout << endl;
  cout << "Energy Above " << 0 << endl;
  printHMGrid(results);
  for(Result result: results){
      amom->Fill(getTheta(result.actual, result.predicted));
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
  for(Result result:results){

  }
  
  _rootfile->Write();
}
