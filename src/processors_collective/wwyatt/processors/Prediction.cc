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
    // usually a good idea to
    //printParameters() ;
    _prediction = new TH2F("predict", "Predicted Angle of Scatter, Correct vs Incorrect Kinematics", 1000, 0.0, 0.01, 1000, 0.0, 0.01);
    _p_theta = new TH1F("p_theta", "Theta between positron and hadronic system", 360, 0, 3.5);
    _e_theta = new TH1F("e_theta", "Theta between positron and hadronic system", 360, 0, 3.5);
    _vector = new TH1F("vector", "Vector", 200, 0.0, 0.05);
    zmom=new TH1F("zmom", "Z-Momentum Distribution", 300, 0, 500);
    tmom=new TH1F("tmom", "T-Momentum Distribution", 300, 0, 10);

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
    for(int i=0; i < col->getNumberOfElements(); ++i){
      MCParticle* particle=dynamic_cast<MCParticle*>(col->getElementAt(i));
      if(particle->getGeneratorStatus()!=1)continue;
      zmom->Fill(abs(particle->getMomentum()[2]));
      tmom->Fill(getTMag(particle->getMomentum()));
    }
    map<int,MCParticle*> max=maxParticle(col, {-11, 11});
    Will::measure data = Will::getMeasure(col);
    if(data.scattered){
      fourvec
	electron=data.electron,
	positron=data.positron,
	hadronic=data.hadronic,
	electronic=data.electronic;
      //      cout << "Z- " << positron.z << endl;
      //zmom->Fill(abs(positron.z));
      //tmom->Fill(positron.t);
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

      //The bool is true if the particle hits.
      
      bool predicted_electron=get_hitStatus(p.electron)<3;
      bool predicted_positron=get_hitStatus(p.positron)<3;
      bool actual_positron=get_hitStatus(positron)<3;
      bool actual_electron=get_hitStatus(electron)<3;
      
      /*
      bool predicted_electron=get_hitStatus(p.electron)<3;
      bool predicted_positron=get_hitStatus(p.positron)<3;
      bool actual_positron=get_hitStatus(getFourVector(max[-11]))<3;
      bool actual_electron=get_hitStatus(getFourVector(max[11]))<3;
      */

      //cout << "p:e = " << positron.T << ":" << electron.T<< endl;
      ++total;
      
      if(data.p_scatter){
	if(actual_positron&&predicted_positron)ph_th++;
	else if( actual_positron && !predicted_positron)ph_tm++;
	else if(!actual_positron &&  predicted_positron)pm_th++;
	else if(!actual_positron && !predicted_positron)pm_tm++;
	else cout << "err2" << endl;
	p_scatter++;
      }else if(data.e_scatter){
	if(actual_electron&&predicted_electron)ph_th++;
	else if( actual_electron && !predicted_electron)ph_tm++;
	else if(!actual_electron &&  predicted_electron)pm_th++;
	else if(!actual_electron && !predicted_electron)pm_tm++;	
	else cout << "err3" << endl;
	e_scatter++;
      }else{
	cout << "err" << endl;
      }
      
    }
}

void Prediction::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void Prediction::end(){ 
  Will::META meta = Will::getMETA();
  cout << interest << endl;
  _rootfile->Write();
  cout << "# Events: " << _nEvt << endl;
  cout << "p_scatter:e_scatter = " << p_scatter << ":"<<e_scatter<<endl;
  cout << "ph_th : " << (1.0*ph_th) << endl;
  cout << "ph_tm : " << (1.0*ph_tm) << endl;
  cout << "pm_th : " << (1.0*pm_th) << endl;
  cout << "pm_tm : " << (1.0*pm_tm) << endl;
  cout << "hit-miss-out-in: " <<meta.STATS[1]<< "-"<<meta.STATS[2]<<"-"<<meta.STATS[3]<<"-"<<meta.STATS[4]<<endl;
  cout << "Scat:noScat " << meta.SCATTERS << ":" << meta.NOSCATTERS << endl;
  cout << "MSC : " << meta.MSC << endl;
}
