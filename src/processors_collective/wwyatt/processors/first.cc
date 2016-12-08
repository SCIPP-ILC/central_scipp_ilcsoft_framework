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
static TH2F* _hitmap;
static TH1F* _energy;
static TH1F* _phi;
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
    _rootfile = new TFile("will_test.root","RECREATE");
    _energy = new TH1F("energy", "Energy", 520.0,  0.0, 260.0);
    _hitmap = new TH2F("pos", "Position Distribution", 300.0, -150.0, 150.0, 300.0, -150.0, 150.0);
    _phi = new TH1F("phi", "Angle Difference", 300, .1, 1);

    //printParameters() ;
    _nEvt = 0 ;

    //Setting up Plots\\
    cout << "Vector allocated." << endl;
}



void first::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 

//Namespace for all my bahbah processing. NOTE: PUT INTO SEPERATE FILE.
namespace BB{
  //Setting PDG constants to particle names (DID NOT VERIFY)
  const static int PHOTON = 22;
  const static int POSITRON = 11;
  const static int ELECTRON = -11;

  //Setting errors static so I can spell check them easily.
  const static string ERROR_NOT_BAHBAH_PARTICLE = "Not a Bahbah particle User-Error: Trying to add a particle that is not accounted for.";
  const static string ERROR_MAX_PHOTONS = "Max Photon User-Error: Trying to add photon but it is already full.";
  const static string ERROR_ALREADY_PARTICLE = "Particle Already There User-Error: Trying to add a paticle but the space is not NULL.";

  //Making a bundle to store and process my event because the data is being given to me nicly.
  class Bundle{
  public:
    Bundle(){}//No constructor

    //Returns the total number of particles added. MAX is 4 right now.
    int getCount(){
      return count;
    }

    //Checks to see if particle is already there, then adds it for processing. Once full it runs initialize.
    void addParticle(MCParticle* _input){
      switch(_input->getPDG()){
      case(PHOTON):
	addPhoton(_input);
	break;
      case(POSITRON):
	if(positron == NULL){
	  positron = _input;
	}else err(ERROR_ALREADY_PARTICLE);
	break;
      case(ELECTRON):
	if(electron == NULL){
	  electron = _input;
	}else err(ERROR_ALREADY_PARTICLE);
	break;
      default:
	err(ERROR_NOT_BAHBAH_PARTICLE);
      } 
      if(++count == 4){
	initialize();
      }
    }
    
    //Where I do my math to output to graph.
    void initialize(){
      double p_phi = getPhi(positron);
      double e_phi = getPhi(electron);
      double del_phi = p_phi-e_phi;

      //Put in graph
      _phi->Fill(del_phi);
    }
    //returns the angle from the x-y axis
    double getTheta(MCParticle* _input){
      return atan(_input->getMomentum()[0]/_input->getMomentum()[1]);
    }
    //returns the angle from the z axis
    double getPhi(MCParticle* _input){
      //Get the norm in the x,y plane:
      const double m_x = (_input->getMomentum())[0];
      const double m_y = (_input->getMomentum())[1];
      const double m_z = _input->getMomentum()[2];
      double norm = sqrt(m_x*m_x + m_y*m_y);
      double phi = atan(norm/m_z);

      //The norm in the opposite side, and the z axis is the same side.
      return phi;
    }
  private:
    int count = 0;
    MCParticle* photonA = NULL;
    MCParticle* photonB = NULL;
    MCParticle* positron = NULL;
    MCParticle* electron = NULL;
    MCParticle* photons[2] = {NULL, NULL};

    //Preparing for updating class for n photons. Currently max is 2 photons.
    void addPhoton(MCParticle* _input){
      if(photonA == NULL){
	photonA = _input;
	photons[0] = photonA;
      }else if(photonB ==NULL){
	photonB = _input;
	photons[1] = photonB;
      }else err(ERROR_MAX_PHOTONS);
    }

    //Helper function, can be set to the namespace for general use. I could not find println(), so I made it.
    void err(string _input){
      cout << _input << endl;
    }
  };
}
void first::processEvent( LCEvent * evt ) { 
    LCCollection* col = evt->getCollection( _colName );
    
    BB::Bundle* bundle = new BB::Bundle; //Created a bundle to process the data for the bahbah event.

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

            }//end final state   
        }//end for
    }
    delete[] bundle;
    _nEvt ++ ;
}



void first::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void first::end(){ 
    _rootfile->Write();
}
