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

#include "SusySVMAnalysis.h"
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

SusySVMAnalysis SusySVMAnalysis;

static TFile* _rootfile;
static TH2F* _hitmap;
static TH1F* _mass;
static TH1F* _scalar;
static TH1F* _vector;
static TH1F* _neutrinos;

SusySVMAnalysis::SusySVMAnalysis() : Processor("SusySVMAnalysis") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
}



void SusySVMAnalysis::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("SusySVMAnalysis.root","RECREATE");
    _V_n_C = new TH1F("V_n_C","Detected Vector",40,0,20);
    _V_n_A = new TH1F("V_n_A","Detectable Vector",40,0,20);
    _V_N_A = new TH1F("V_N_A","True Vector",40,0,20);

    _S_n_A = new TH1F("S_n_C","Detected Scalar",40,0,20);
    _S_n_A = new TH1F("S_n_A","Detectable Scalar",40,0,20);
    _S_N_A = new TH1F("S_N_A","True Scalar",40,0,20);
 
  
    _M_n_A = new TH1F("M_n_C","Detected Mass",40,0,20);
    _M_n_A = new TH1F("M_n_A","Detectable Mass",40,0,20);
    _M_N_A = new TH1F("M_N_A","True Mass",40,0,20);
    // usually a good idea to
    //printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;
}



void SusySVMAnalysis::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void SusySVMAnalysis::processEvent( LCEvent * evt ) { 
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
                double P = particle->getMomentum().magnitude();
                double pz = particle->getPZ();
                double px = particle->getPX();
                double py = particle->getPY();
                double cos = pz/P;
                double scalar = sqrt(px*px+py*py); 
                bool isNeutrino = (
                    id == 12 || id == -12 ||
                    id == 14 || id == -14 ||
                    id == 16 || id == -16 ||
                    id == 18 || id == -18);
                bool isForward = ( cos > 0.9 || cos < -0.9);               
                scalars[0]+=scalar; //true
                vectors[0][0]+=px;
                vectors[0][1]+=py;
                vectors[0][2]+=pz;
                energy[0]+=E;                        
                if(!isDarkMatter && !isNeutrino){ //detectable
                    scalars[2]+=scalar;
                    vectors[2][0]+=px;
                    vectors[2][1]+=py;
                    vectors[2][2]+=pz;
                    energy[2]+=E;
                    if(!isForward){
                        scalars[1]+=scalar; //detected
                        vectors[1][0]+=px;
                        vectors[1][1]+=py;
                        vectors[1][2]+=pz;
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
            (vectors[0][0]*vectors[0][0]+vectors[0][1]vectors[0][1]+
            vectors[0][2]*vectors[0][2]);
        double total_true_mass = sqrt(total_true_mass_squared);
        double total_detected_mass_squared = energy[1]*energy[1]-
            (vectors[1][0]*vectors[1][0]+vectors[1][1]*vectors[1][1]+
            vectors[1][2]*vectors[1][2]);
        double total_detected_mass = sqrt(total_detected_mass_squared);
        double total_detectable_mass_squared = energy[2]*energy[2]-
            (vectors[2][0]*vectors[2][0]+vectors[2][1]*vectors[2][1]+
            vectors[2][2]*vectors[2][2]);
            
        double total_true_vector = sqrt(vectors[0][0]*vectors[0][0]+vectors[0][1]*vectors[0][1]+vectors[0][2]*vectors[0][2]); 
        double total_detected_vector = sqrt(vectors[1][0]*vectors[1][0]+vectors[1][1]*vectors[1][1]+vectors[1][2]*vectors[1][2]); 
        double total_detectable_vector = sqrt(vectors[2][0]*vectors[2][0]+vectors[2][1]*vectors[2][1]+vectors[2][2]*vectors[2][2]); 
        _V_n_C->Fill(total_detected_vector);
        _V_n_A->Fill(total_detectable_vector);
        _V_N_A->Fill(total_true_vector);
    }
    _nEvt ++ ;
}



void SusySVMAnalysis::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void SusySVMAnalysis::end(){ 

    _rootfile->Write();
}
