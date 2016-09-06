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
 * author Summer Zuber
 * August 7, 2016
 * This processor analyses the distribution of S, V and M 
 * observables of SUSY events 
 */

#include "Thrust.h"
#include "scipp_ilc_utilities.h"
#include <iostream>
#include <cmath>

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

Thrust Thrust;

static TFile* _rootfile;

static TH1F* _V_n_C;
static TH1F* _V_n_A;
static TH1F* _V_N_A;
static TH1F* _S_n_C;
static TH1F* _S_n_A;
static TH1F* _S_N_A;
static TH1F* _M_n_C;
static TH1F* _M_n_A;
static TH1F* _M_N_A;

Thrust::Thrust() : Processor("Thrust") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
}



void Thrust::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("Thrust.root","RECREATE");
    _V_n_C = new TH1F("V_n_C","Detected Vector",40,0,20);
    _V_n_A = new TH1F("V_n_A","Detectable Vector",40,0,20);
    _V_N_A = new TH1F("V_N_A","True Vector",40,0,20);

    _S_n_C = new TH1F("S_n_C","Detected Scalar",40,0,20);
    _S_n_A = new TH1F("S_n_A","Detectable Scalar",40,0,20);
    _S_N_A = new TH1F("S_N_A","True Scalar",40,0,20);
 
  
    _M_n_C = new TH1F("M_n_C","Detected Mass",40,0,20);
    _M_n_A = new TH1F("M_n_A","Detectable Mass",40,0,20);
    _M_N_A = new TH1F("M_N_A","True Mass",40,0,20);
    // usually a good idea to
    //printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;
}



void Thrust::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void Thrust::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    LCCollection* col = evt->getCollection( _colName ) ;
    cout << endl;
    cout << endl;
    cout << endl;
    cout << "event = " << _nEvt << endl;
    
    double vec[4][3];
    double scalars[4];
    double energy[4]; 

    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
        
        // For each particle in Event ...
        for(int particleIndex = 0; particleIndex < nElements ; particleIndex++){
           MCParticle* particle = dynamic_cast<MCParticle*>( col->getElementAt(particleIndex) );
            
           int id = particle->getPDG(); 
           int stat = particle->getGeneratorStatus();
           // If Particle is FINAL-STATE 
           if(stat==1){
                bool isDarkMatter = (id == 1000022);
                if(isDarkMatter) continue ;
                double E = particle->getEnergy();
                const double* P = particle->getMomentum();
                double Pmag = sqrt(P[0]*P[0]+P[1]*P[1]+P[2]*P[2]);
                double cos = P[2]/Pmag;
                double scalar = sqrt(P[0]*P[0]+P[1]*P[1]); 
                bool isNeutrino = (
                    id == 12 || id == -12 ||
                    id == 14 || id == -14 ||
                    id == 16 || id == -16 ||
                    id == 18 || id == -18);
                bool isForward = ( cos > 0.9 || cos < -0.9);               
                scalars[0]+=scalar; //true
                vec[0][0]+=P[0];
                vec[0][1]+=P[1];
                vec[0][2]+=P[2];
                energy[0]+=E;                        
                if(!isDarkMatter && !isNeutrino){ //detectable
                    scalars[2]+=scalar;
                    vec[2][0]+=P[0];
                    vec[2][1]+=P[1];
                    vec[2][2]+=P[2];
                    energy[2]+=E;
                    if(!isForward){
                        scalars[1]+=scalar; //detected
                        vec[1][0]+=P[0];
                        vec[1][1]+=P[1];
                        vec[1][2]+=P[2];
                        energy[1]+=E;      
                    }
                }
                 
           }//end final state
        }//end for

        //all
        double total_true_scalar = scalars[0];
        double total_detected_scalar = scalars[1];
        double total_detectable_scalar = scalars[2];

        double total_true_mass_squared = energy[0]+energy[0]-
            (vec[0][0]*vec[0][0]+vec[0][1]*vec[0][1]+
            vec[0][2]*vec[0][2]);
        double total_true_mass = sqrt(total_true_mass_squared);
        double total_detected_mass_squared = energy[1]*energy[1]-
            (vec[1][0]*vec[1][0]+vec[1][1]*vec[1][1]+
            vec[1][2]*vec[1][2]);
        double total_detected_mass = sqrt(total_detected_mass_squared);
        double total_detectable_mass_squared = energy[2]*energy[2]-
            (vec[2][0]*vec[2][0]+vec[2][1]*vec[2][1]+
            vec[2][2]*vec[2][2]);
            
        double total_true_vector = sqrt(vec[0][0]*vec[0][0]+vec[0][1]*vec[0][1]+vec[0][2]*vec[0][2]); 
        double total_detected_vector = sqrt(vec[1][0]*vec[1][0]+vec[1][1]*vec[1][1]+vec[1][2]*vec[1][2]); 
        double total_detectable_vector = sqrt(vec[2][0]*vec[2][0]+vec[2][1]*vec[2][1]+vec[2][2]*vec[2][2]); 
        _V_n_C->Fill(total_detected_vector);
        _V_n_A->Fill(total_detectable_vector);
        _V_N_A->Fill(total_true_vector);
    }
    _nEvt ++ ;
}



void Thrust::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void Thrust::end(){ 

    _rootfile->Write();
}
