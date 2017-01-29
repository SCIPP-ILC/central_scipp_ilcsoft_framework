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

#include "first.h"
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


first first;

static TFile* _rootfile;
static TH2F* _hitmap1;
static TH2F* _hitmap2;
static TH2F* _hitmap3;
static TH2F* _hitmap4;

static TH2F* _e_hitmap;
static TH2F* _p_hitmap;

static TH1F* _hitmiss;

static int hit = 0;
static int miss = 0;


//static TH1F* _energy;
//static TH1F* _phi1;
//static TH1F* _phi2;
//static TH1F* _phi3;
static TH1F * _cosE;
static TH1F * _cosP;
static TH1F * _cos;

class Bundle{
public:
  const static bool VERBOSE = false;
  //Setting PDG constants to particle names (DID NOT VERIFY)
  const static int PHOTON = 22;
  const static int POSITRON = 11;
  const static int ELECTRON = -11;

  //Setting errors static so I can spell check them easily.
  const string ERROR_NOT_BAHBAH_PARTICLE = "Not a Bahbah particle User-Error: Trying to add a particle that is not accounted for.";
  const string ERROR_MAX_PHOTONS = "Max Photon User-Error: Trying to add photon but it is already full.";
  const string ERROR_ALREADY_PARTICLE = "Particle Already There User-Error: Trying to add a paticle but the space is not NULL.";

 
  MCParticle* photonA = NULL;
  MCParticle* photonB = NULL;
  MCParticle* positron = NULL;
  MCParticle* electron = NULL;
  MCParticle* photons[2] = {NULL, NULL};

  Bundle(){
    positron = NULL;
    electron = NULL;
  }

  void err(string _input);

  void graphHitStatus(double x, double y);

  //Returns the total number of particles added. MAX is 4 right now.
  int getCount();

  MCParticle* getElectron();
  MCParticle* getPositron();

  //-- Physics Section --\\
  //returns the angle from the x-y axis
  double getPhi(MCParticle* _input);

  //returns the angle from the z axis
  double getTheta(MCParticle* _input,bool);

  //Checks to see if particle is already there, then adds it for processing. Once full it runs initialize.
  void addParticle(MCParticle* _input);

  //Preparing for updating class for n photons. Currently max is 2 photons.
  void addPhoton(MCParticle* _input);

  //Gets the norm of the vector and finds it's manitude.
  double getMagnitude(MCParticle* _input);

  double getDotProduct(MCParticle* A, MCParticle* B);

};

//Prototypes\\
void graphHitStatus(double x, double y);
void initialize(Bundle* b);



first::first() : Processor("first") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}


void first::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;
    cout << "Initialized" << endl;
    _rootfile = new TFile("Phi_Bhabha.root","RECREATE");
    //    _energy = new TH1F("energy", "Energy", 520.0,  0.0, 260.0);
    _hitmiss = new TH1F("hm", "Hit Miss Ratio", 5, -1, 2);
    _hitmap1 = new TH2F("pos1", "Position Distribution On Beamcal", 300.0, -150.0, 150.0, 300.0, -150.0, 150.0);
    _hitmap2 = new TH2F("pos2", "Position Distribution On Beamcal", 300.0, -150.0, 150.0, 300.0, -150.0, 150.0);
    _hitmap3 = new TH2F("pos3", "Position Distribution On Beamcal", 300.0, -150.0, 150.0, 300.0, -150.0, 150.0);
    _hitmap4 = new TH2F("pos4", "Position Distribution On Beamcal", 300.0, -150.0, 150.0, 300.0, -150.0, 150.0);
    //    static TH2F* _e_hitmap;
    _e_hitmap = new TH2F("hitmap", "Position Distribution On Beamcal", 300.0, -150.0, 150.0, 300.0, -150.0, 150.0);
    //    _hitmap4 = new TH2F("pos4", "Position Distribution On Beamcal", 300.0, -150.0, 150.0, 300.0, -150.0, 150.0);
    //    _theta1 = new TH1F("theta_s", "Angle Difference (Electron - Positron) From 0-2π", 500, 0, 6.28);
    //    _theta2 = new TH1F("theta_m", "Angle Difference (Electron - Positron) From 1-π", 500, 1, 3.14);
    //    _theta3 = new TH1F("theta_t", "Angle Difference (Electron - Positron) From 0-1rad", 500, 0, .1);
    _cosE = new TH1F("electron_corth", "Cosine Theta Distribution For Electrons", 350, 0, 6.28);
    _cosP = new TH1F("positron_costh", "Cosine Theta Distribution For Positrons", 350, 0, 6.28);
    _cos = new TH1F("costh", "Cosine Theta Distribution", 350, 0, 6.28);
    //printParameters() ;
    _nEvt = 0 ;

    //Setting up Plots\\
    cout << "Vector allocated." << endl;
}



void first::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 


void first::processEvent( LCEvent * evt ) { 
    LCCollection* col = evt->getCollection( _colName );
    
    Bundle* bundle = new Bundle; //Created a bundle to process the data for the bahbah event.

    int stat, id =0;
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );        
            id = hit->getPDG();
            stat = hit->getGeneratorStatus();
            if(stat==1){
	      //hit is an end particle:
	      bundle->addParticle(hit); //I add it to a Bahbah class to do the rest of the work.
	      if(bundle->getPositron() != NULL && bundle->getElectron() != NULL){
		initialize(bundle);
	      }
            }//end final state   
        }//end for
    }
    delete bundle;
    _nEvt ++ ;
}



void first::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void first::end(){ 
  cout << "\n Hits: " << hit << endl;
  cout << " Miss: " << miss << endl;
  _rootfile->Write();
}


//Implementation of the Bundle Class
void initialize(Bundle* b){  
  //CosineTheta Not doing center to mass frame.
  double cosE = cos(b->getTheta(b->getElectron(), true));
  double cosP = cos(b->getTheta(b->getPositron(), true));
  _cosE->Fill(cosE);
  _cosP->Fill(cosP);
  _cos->Fill(cosE);
  _cos->Fill(cosP);
}

void init_hitmap_general(Bundle* b){
  double p_px = b->positron->getMomentum()[0];
  double p_E = b->positron->getEnergy();
  double e_px = b->electron->getMomentum()[0];
  double e_E = b->electron->getEnergy();

  scipp_ilc::transform_to_lab(p_px, p_E, p_px, p_E);
  scipp_ilc::transform_to_lab(e_px, e_E, e_px, e_E);
  
  //Put in graphp
  //  graphHitStatus(e_px, b->electron->getMomentum()[1]);
  //  graphHitStatus(p_px, b->positron->getMomentum()[1]);

  _e_hitmap->Fill(p_px, b->positron->getMomentum()[1]);
  _e_hitmap->Fill(e_px, b->electron->getMomentum()[1]);
}

void graphHitStatus(double x, double y, TH1F* hm = NULL, TH2F* hm1 = NULL, TH2F* hm2 = NULL, TH2F* hm3 = NULL, TH2F* hm4 = NULL){
  switch(scipp_ilc::get_hitStatus(x, y)){
  case(1):
    hit++;
    hm->Fill(1);
    hm1->Fill(x, y);
    break;
  case(2):
    miss++;
    hm->Fill(0);
    hm2->Fill(x, y);
    break;
  case(3):
    miss++;
    hm->Fill(0);
    hm3->Fill(x, y);
    break;
  case(4):
    miss++;
    hm->Fill(0);
    hm4->Fill(x, y);
    break;
  }
}

//Function used to calculate anglee difference using norms and phi to get theta.
void init_angleBetweenEP_1(){
  /*  //  Trial 1: Getting angle between electron and positron.
      double p_theta = getTheta(positron);
      double e_theta = getTheta(electron);
      double del_theta = p_theta-e_theta; 


  //  _hitmap->Fill(positron->getMomentum()[0], positron->getMomentum()[1]);
  //  _hitmap->Fill(electron->getMomentum()[0], electron->getMomentum()[1]);
  //  _theta1->Fill(del_theta);
  //  _theta2->Fill(del_theta);
  //  _theta3->Fill(del_theta);
  */
}

//Function used to calculate the angle difference using dot products.
void init_angleBetweenEP_2(){
/*// Trial 2: Getting angle between elctron and positron using dot products.
double dot = getDotProduct(positron, electron);
double mag_A = getMagnitude(positron);
double mag_B = getMagnitude(electron);
double val = -dot/(mag_A*mag_B);
double del_theta = acos(val); */
}
void Bundle::err(string _input){
  cout << _input << endl;
}

//-- Physics Section --\\
//returns the angle from the x-y axis
double Bundle::getPhi(MCParticle* _input){
  return atan(_input->getMomentum()[0]/_input->getMomentum()[1]);
}
//returns the angle from the z axis
double Bundle::getTheta(MCParticle* _input, bool lab = false){
  //Get the norm in the x,y plane:
  const double m_x = (_input->getMomentum())[0];
  const double m_y = (_input->getMomentum())[1];
  const double m_z = _input->getMomentum()[2];
  double p_x = m_x;
  if(lab){
    double energy = _input->getEnergy();
    scipp_ilc::transform_to_lab(p_x, energy, p_x, energy);
  }
  double norm = sqrt(p_x*p_x + m_y*m_y);
  double theta = atan(norm/m_z);

  //The norm in the opposite side, and the z axis is the same side.
  return theta;

}

//-- CompSci Section --\\

MCParticle* Bundle::getElectron(){
  //  if(electron == NULL)cout << "\nElectron is NUll\n";
  return electron;
}

MCParticle* Bundle::getPositron(){
  //  if(electron == NULL)cout << "\nPositron is NUll\n";
  return positron;
}


//Checks to see if particle is already there, then adds it for processing. Once full it runs initialize.
void Bundle::addParticle(MCParticle* _input){
  //  if(positron == NULL) cout << "pNULL" << endl;
  //  if(electron == NULL) cout << "eNULL" << endl;
  switch(_input->getPDG()){
  case(PHOTON):
    addPhoton(_input);
    break;
  case(POSITRON):
    if(positron == NULL){
      positron = _input;
    }else if(VERBOSE) err(ERROR_ALREADY_PARTICLE);
    break;
  case(ELECTRON):
    if(electron == NULL){
      electron = _input;
    }else if(VERBOSE)err(ERROR_ALREADY_PARTICLE);
    break;
  default:
    if(VERBOSE)err(ERROR_NOT_BAHBAH_PARTICLE);
  }
}
    
//Preparing for updating class for n photons. Currently max is 2 photons.
void Bundle::addPhoton(MCParticle* _input){
  if(photonA == NULL){
    photonA = _input;
    photons[0] = photonA;
  }else if(photonB ==NULL){
    photonB = _input;
    photons[1] = photonB;
  }else if(VERBOSE)err(ERROR_MAX_PHOTONS);
}

//Gets the norm of the vector and finds it's manitude.
double Bundle::getMagnitude(MCParticle* _input){
  const double m_x = (_input->getMomentum())[0];
  const double m_y = (_input->getMomentum())[1];
  const double m_z = _input->getMomentum()[2];
  return sqrt(m_x*m_x + m_y*m_y + m_z*m_z);
}
double Bundle::getDotProduct(MCParticle* A, MCParticle* B){
  const double a_x = (A->getMomentum())[0];
  const double a_y = (A->getMomentum())[1];
  const double a_z = A->getMomentum()[2];
  const double b_x = (B->getMomentum())[0];
  const double b_y = (B->getMomentum())[1];
  const double b_z = B->getMomentum()[2];
  return (a_x*b_x + a_y*b_y + a_z*b_z);
}
