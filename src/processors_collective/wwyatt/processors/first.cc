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
static TH2F* _hh;
static TH2F* _mm;
static TH2F* _hm;
static TH1F* _cos;
static TH1F* _cosp;
static TH1F* _cose;
static bool inlab = false;
//th2f hitmap
//th1f histogram

first::first() : Processor("first") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
    
    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}


void first::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;
    
    cout << "Initialized " << endl;
    _rootfile = new TFile("Phi_Bhabha.root","RECREATE");
    //These initialization are here as reference until I finish rewriting all the code.
    _hh = new TH2F("hh", "Hit-Hit Distribution", 300.0, -150.0, 150.0, 300.0, -150.0, 150.0);
    _mm = new TH2F("mm", "Miss-Miss Distribution", 300.0, -150.0, 150.0, 300.0, -150.0, 150.0);
    _hm = new TH2F("hm", "Hit-Miss Distribution", 300.0, -150.0, 150.0, 300.0, -150.0, 150.0);
    _cos = new TH1F("Theta", "Theta Distribution", 300, 0, 3.4);
    _cose = new TH1F("ETheta", "Theta Electron Distribution", 300, 0,3.4);
    _cosp = new TH1F("PTheta", "Theta Positron Distribution", 300, 0,3.4);

    //    _energy = new TH1F("energy", "Energy", 520.0,  0.0, 260.0);
    /*    _hitmiss = new TH1F("hm", "Hit Miss Ratio", 5, -1, 2);
    _hitmap1 = new TH2F("pos1", "Position Distribution On Beamcal", 300.0, -150.0, 150.0, 300.0, -150.0, 150.0);




    //    _hitmap4 = new TH2F("pos4", "Position Distribution On Beamcal", 300.0, -150.0, 150.0, 300.0, -    //    _theta1 = new TH1F("theta_s", "Angle Difference (Electron - Positron) From 0-2π", 500, 0, 6.28
    //    _theta2 = new TH1F("theta_m", "Angle Difference (Electron - Positron) From 1-π", 500, 1, 3.14)
    //    _theta3 = new TH1F("theta_t", "Angle Difference (Electron - Positron) From 0-1rad", 500, 0, .1);
    _cosE = new TH1F("electron_corth", "Cosine Theta Distribution For Electrons", 350, 0, 6.28);

    */
    _nEvt = 0 ;

    //Setting up Plots\\
    cout << "Vector allocated." << endl;
}



void first::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 


void first::processEvent( LCEvent * evt ) { 
    LCCollection* col = evt->getCollection( _colName );
    //    cout << "NEW Event Started." << endl;
    info = 0;
    electron=NULL;
    positron=NULL;
    int stat, id =0;
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );        
            id = hit->getPDG();
            stat = hit->getGeneratorStatus();
            if(stat==1){
	      //hit is an end particle:
	      addParticle(hit); //The bundle object has the rest of the processing.
	      if(getPositron() != NULL && getElectron() != NULL){
		init_hitmap(inlab);
		double e = acos(getCosTheta(getElectron(), inlab));
		double p = acos(getCosTheta(getPositron(), inlab));
		_cose->Fill(e);
		_cosp->Fill(p);
		_cos->Fill(e+p);
		electron=NULL;
		positron=NULL;

	      }
            }//end final state   
        }//end for
    }
    _nEvt ++ ;
}



void first::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void first::end(){ 
  //Mostly useless output, just to check that stuff is being processed.
  print("Hits: " + to_string(electronStore.hit + positronStore.hit));
  print("Miss: " + to_string(electronStore.miss + positronStore.miss));
  print("hh: " + to_string(electronStore.hh.size()));
  print("hm: " + to_string(electronStore.hm.size()));
  print("mh: " + to_string(electronStore.mh.size()));
  print("mm: " + to_string(electronStore.mm.size()));

  //Making thoes hit-miss plots on the data. 
  plotHitMiss(_hh, {electronStore.hh, positronStore.hh});
  plotHitMiss(_mm, {electronStore.mm, positronStore.mm});
  plotHitMiss(_hm, {electronStore.hm, positronStore.hm, electronStore.mh, positronStore.mh});
  _rootfile->Write();
}

//void plotHitMiss(TH2F* p, vector<double*> *m){

//}

//Implementation of the Bundle Class
void first::initialize(){
  /*//CosineTheta Not doing center to mass frame.
  double cosE = cos(b->getTheta(b->getElectron(), true));
  double cosP = cos(b->getTheta(b->getPositron(), true));
  _cosE->Fill(cosE);
  _cosP->Fill(cosP);
  _cos->Fill(cosE);
  _cos->Fill(cosP);*/

}
void first::init_hitmap(bool lab){
  graphHitStatus(getMomentum(getElectron(), lab), ELECTRON);
  graphHitStatus(getMomentum(getPositron(), lab), POSITRON);
}

void first::graphHitStatus(const double*  momentum, int id, bool lab){
  double x = momentum[0];
  double y = momentum[1];
  switch(scipp_ilc::get_hitStatus(x, y)){
  case(1): //hit beamcal
    if(id == this->ELECTRON){
      info += 10;
    }else if(id==this->POSITRON){
      info += 100;
    }
    info += 1;
    break;
  case(2): //outside beamcal
    if(id == this->ELECTRON){
      info += 10;
    }else if(id==this->POSITRON){
      info += 100;
    }
    info += 1;
    break;
  case(3): //outgoing beampipe
    if(id == this->ELECTRON){
      info += 70;
    }else if(id==this->POSITRON){
      info += 700;
    }
    info += 1;
    break;
  case(4): //incoming beampipe
    if(id == this->ELECTRON){
      info += 70;
    }else if(id==this->POSITRON){
      info += 700;
    }
    info += 1;
    break;
  }
  //cout << "# " << info << endl;
  // 1=hit; 7=miss; positron>electron;
  if(info%10==2){
    // Electron-Positron
    // Hit-Hit
    if(info%100/10==1 && info%1000/100==1){ 
      electronStore.hit++;
      positronStore.hit++;
      electronStore.hh.push_back(getMomentum(getElectron(), lab));
      positronStore.hh.push_back(getMomentum(getPositron(), lab));
    }
    // Miss-Miss
    else if(info%100/10==7 && info%1000/100==7){
      electronStore.miss++;
      positronStore.miss++;
      electronStore.mm.push_back(getMomentum(getElectron(), lab));
      positronStore.mm.push_back(getMomentum(getPositron(), lab));
    }
    // Miss-Hit
    else if(info%100/10==7 && info%1000/100==1){
      electronStore.miss++;
      positronStore.hit++;
      electronStore.mh.push_back(getMomentum(getElectron(), lab));
      positronStore.mh.push_back(getMomentum(getPositron(), lab));
    }
    // Hit-Miss
    else if(info%100/10==1 && info%1000/100==7){
      electronStore.hit++;
      positronStore.miss++;
      electronStore.hm.push_back(getMomentum(getElectron(), lab));
      positronStore.hm.push_back(getMomentum(getPositron(), lab));
    }
  }
}

void first::print(string _input){
  cout << _input << endl;
}

//Used to plot one momentum vector on a TH2F plot.
/*void first::plotHitMiss(TH2F* p, vector<double*> momentum){
    plotHotMiss(p, {momentum});
    }*/

 //Used to plot multiple vectors on a plot.
void first::plotHitMiss(TH2F* p, initializer_list<vector<double*>> momentums){
    for(vector<double*> arr: momentums){
      for(double* m: arr){
	p->Fill(m[0], m[1]);
      }
    }
  }


//-- Physics Section --\\

//Function used to calculate the angle difference (colinearity) using dot products.
double first::get_colinearity(bool lab){
  //Getting angle between electron and positron using dot products.
  double dot   = getDotProduct(getPositron(), getElectron());
  double mag_B = getMagnitude( getElectron());
  double mag_A = getMagnitude( getPositron());
  double val = dot/(mag_A*mag_B);
  double del_theta = acos(val); 
  return del_theta;
}

//returns the angle from the x-y axis
 double first::getPhi(MCParticle* _input, bool lab){
   return atan(getMomentum(_input, lab)[0]/getMomentum(_input, lab)[1]);
 }

//returns the angle from the z axis
 double first::getCosTheta(MCParticle* _input, bool lab){
   double* m = getMomentum(_input, lab);
   return m[2]/getMagnitude(m);
 }

 double first::getTheta(MCParticle* _input, bool lab){
  //Get the norm in the x,y plane:
  double* m = getMomentum(_input, lab);
  const double m_x = m[0];
  const double m_y = m[1];
  const double m_z = m[2];
  double p_x = m_x;
  double norm = sqrt(p_x*p_x + m_y*m_y);
  double theta = atan(norm/m_z);

  //The norm in the opposite side, and the z axis is the same side.
  return theta;

}

//-- CompSci Section --\\

MCParticle* first::getElectron(){
  //  if(electron == NULL)cout << "\nElectron is NUll\n";
  return electron;
}

MCParticle* first::getPositron(){
  //  if(electron == NULL)cout << "\nPositron is NUll\n";
  return positron;
}


//Checks to see if particle is already there, then adds it for processing.
//This is overkill, I only needed the photon and electron.
void first::addParticle(MCParticle* _input){
  switch(_input->getPDG()){
  case(first::PHOTON):
    addPhoton(_input);
    break;
  case(first::POSITRON):
    if(positron == NULL){
      positron = _input;
    }else if(VERBOSE) print(ERROR_ALREADY_PARTICLE);
    break;
  case(first::ELECTRON):
    if(electron == NULL){
      electron = _input;
    }else if(VERBOSE)print(ERROR_ALREADY_PARTICLE);
    break;
  default:
    if(VERBOSE)print(ERROR_NOT_BAHBAH_PARTICLE);
  }
}
    
//Preparing for updating class for n photons. Currently max is 2 photons.
void first::addPhoton(MCParticle* _input){
  if(photonA == NULL){
    photonA = _input;
    photons[0] = photonA;
  }else if(photonB ==NULL){
    photonB = _input;
    photons[1] = photonB;
  }else if(VERBOSE)print(ERROR_MAX_PHOTONS);
}

//Gets the norm of the vector and finds it's manitude.
 double first::getMagnitude(MCParticle* _input, bool lab){
  double* m = getMomentum(_input, lab);
  return getMagnitude(m);
}

 double first::getMagnitude(double*m){
  return sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
}
double first::getDotProduct(MCParticle* A, MCParticle* B, bool lab){
  double* a = getMomentum(A, lab);
  double* b = getMomentum(B, lab);
  return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}
double*  first::getMomentum(MCParticle* _input, bool lab){
  double x = _input->getMomentum()[0];
  double y = _input->getMomentum()[1];
  double z = _input->getMomentum()[2];
  double e = _input->getEnergy();

  if(lab)scipp_ilc::transform_to_lab(e,x,e,x);
  return new double[3]{x,y,z};
}
