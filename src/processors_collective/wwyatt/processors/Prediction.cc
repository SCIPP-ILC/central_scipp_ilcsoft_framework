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
    /*double tot_mom[]={0, 0};
    double compEn_e=0, compEn_p=0;
    double* mom   = new double[4]();//4vec
    double* mom_e = new double[4]();
    double* mom_p = new double[4]();
    double tmom, theta, good_t, bad_t, mag=0.0, eT, pT;
    bool scatter;
    //    MCParticle* electronic
    double* hadronic = new double[4]();
    double* electronic = new double[4]();

    int nElements = col->getNumberOfElements();
    scatter = false;*/

    Will::measure data = Will::getMeasure(col);

    /*mom_e=legacy(measure.electron);
    mom_p=legacy(measure.positron);
    hadronic=legacy(measure.hadronic);
    electronic=legacy(measure.electronic);
    mag=measure.mag;
    eT=measure.electron.T;
    pT=measure.positron.T;
    scatter=measure.scattered;*/
    
    if(data.scattered == true){
      fourvec
	electron=data.electron,
	positron=data.positron,
	hadronic=data.hadronic,
	electronic=data.electronic;
      double mag=data.mag;
      double alpha=(500-hadronic.E -hadronic.z);
      double  beta=(500-hadronic.E +hadronic.z);
      //Legacy
      fourvec predict;
      predict.x=-hadronic.x;
      predict.y=-hadronic.y;
      //END
      prediction p(-hadronic.x,-hadronic.y);
      p.electron.z = -(pow(electron.T, 2)-pow(alpha, 2))/(2*alpha);
      p.positron.z = (pow(positron.T, 2)-pow(beta, 2))/(2*beta);
      
      //      predict[3] = (pow(pT, 2)-pow(beta, 2))/(2*beta);
      //      predict[2] = -(pow(eT, 2)-pow(alpha, 2))/(2*alpha);
      //correct prediction (correct if positron deflection)

      double electronTheta=getTheta(p.electron);
      double positronTheta=getTheta(p.positron);
      if(mag>1.0){
	_prediction->Fill(electronTheta,positronTheta);
      }
      
      //Find a more accurate naming convention for these.
      
      //      double dot_c = electronic[0]*predict[0] + electronic[1]*predict[1] + electronic[2]*predict[3]; //Correct dot
      //      double dot_i = electronic[0]*predict[0] + electronic[1]*predict[1] + electronic[2]*predict[2]; //Incorrect dot

      double p_mag = getMag(electronic);

      //      double p_mag_c = sqrt(pow(predict[0], 2)+pow(predict[1], 2)+pow(predict[3], 2)); //Correct mag
      //      double p_mag_i = sqrt(pow(predict[0], 2)+pow(predict[1], 2)+pow(predict[2], 2)); //Incorrect mag

      //Log for the values above.
      /*
      if(false){
	cout << "dot_i " << dot_i << endl; //200
	cout << "dot_c " << dot_c << endl; //200
	cout << "p-mag " << p_mag << endl; //200
	cout << "p-mag_i " << p_mag_i << endl; //206
	cout << "p-mag_c " << p_mag_c << endl; //206
      }
      */
      double e_theta=getTheta(electronic,p.electron);
      double p_theta=getTheta(electronic,p.positron);
      //double theta_c = acos(dot_c/(p_mag*p_mag_c)); //Correct prediction
      //double theta_i = acos(dot_i/(p_mag*p_mag_i)); //Incorrect prediction

      /*
      //Log for anomoly #Tyco
      if(false){
	double electronUnit=acos(electronic[2]/getMag(electronic));
	//double predEUnit=acos(electronic[2]/getMag(electronic));
	cout << "Electron " << electronUnit << endl;
	cout << "Positron " << 3.14159-acos(mom_p[2]/getMag(mom_p)) << endl;
	cout << "Prediction (electron scatter) " << acos(predict[2]/p_mag_i) << endl;
	cout << "Prediction (positron scatter) " << 3.14159-acos(predict[3]/p_mag_c) << endl;
      }
      */
      _p_theta->Fill(e_theta);
      _e_theta->Fill(p_theta);


      //theta = acos(dot/(e_mag*p_mag)); 
      //cout << "Prediction Efficiency :" <<  theta << endl;
      //      double e_mag = getMag(electronic);
      //      double e_mag = sqrt(pow(electronic[0], 2)+pow(electronic[1], 2)+pow(electronic[2], 2)); 
      //cout << dot_c << endl; //dot_c 59066.6


    }
}

void Prediction::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void Prediction::end(){ 
    cout << interest << endl;
    _rootfile->Write();
}
