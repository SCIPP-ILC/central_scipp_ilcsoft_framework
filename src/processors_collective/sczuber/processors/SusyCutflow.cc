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
 */

#include "SusyCutflow.h"
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

SusyCutflow SusyCutflow;

static TFile* _rootfile;
static TH2F* _hitmap;
static TH1F* _mass;
static TH1F* _scalar;
static TH1F* _vector;
static TH1F* _neutrinos;

SusyCutflow::SusyCutflow() : Processor("SusyCutflow") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
}



void SusyCutflow::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("SusyCutflow.root","RECREATE");
    //_hitmap = new TH2F("hitmap","Hit Distribution",300.0,-150.0,150.0,300.0,-150.0,150.0);
    //_scalar = new TH1F("scalar", "Transverse Momentum Scalar Magnitude", 2000.0, 0.0, 20.0);
    //_vector = new TH1F("vector", "Transverse Momentum Vector Magnitude", 2000.0, 0.0, 20.0);
    //_mass = new TH1F("mass", "Mass Parameter", 2000.0, 0.0, 20.0);
    //_neutrinos = new TH1F("neutrinos", "Neutrinos per Event", 10.0, 0.0,10.0); 
    // usually a good idea to
    //printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;

   // _no_def_count = 0;
   // _e_def_count = 0;
   // _p_def_count = 0;
}



void SusyCutflow::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void SusyCutflow::processEvent( LCEvent * evt ) { 
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
                scalars[0]+=scalar;
                vectors[0][0]+=px;
                vectors[0][1]+=py;
                vectors[0][2]+=pz;
                energy[0]+=E;                        
                if(!isDarkMatter && !isNeutrino){
                    scalars[2]+=scalar;
                    vectors[2][0]+=px;
                    vectors[2][1]+=py;
                    vectors[2][2]+=pz;
                    energy[2]+=E;
                    if(!isForward){
                        scalars[1]+=scalar;
                        vectors[1][0]+=px;
                        vectors[1][1]+=py;
                        vectors[1][2]+=pz;      
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
        cuts[0][0]+=1;
        cuts[1][0]+=1;
        cuts[2][0]+=1;
        if(total_true_scalar > 0.5){
            cuts[0][1]+=1;
            if(total_true_mass > 0.5){
                cuts[0][2]+=1;
                if(total_true_scalar > 1){
                    cuts[0][3]+=1;
                    if(total_true_mass > 1){
                        cuts[0][4]+=1;
                    }
                }
            }
        }

        if(total_detected_scalar > 0.5){
            cuts[1][1]+=1;
            if(total_detected_mass > 0.5){
                cuts[1][2]+=1;
                if(total_detected_scalar > 1){
                    cuts[1][3]+=1;
                    if(total_detected_mass > 1){
                        cuts[1][4]+=1;
                    }
                }
            }
        }
        
        if(total_detectable_scalar > 0.5){
            cuts[2][1]+=1;
            if(total_detectable_mass > 0.5){
                cuts[2][2]+=1;
                if(total_detectable_scalar > 1){
                    cuts[2][3]+=1;
                    if(total_detectable_mass > 1){
                        cuts[2][4]+=1;
                    }
                }
            }
        }
        

        cout << "TRUE" << endl;
        cout<< "cut_0 " << cuts[0][0] << endl;
        cout<< "cut_1 " << cuts[0][1] << endl;
        cout<< "cut_2 " << cuts[0][2] << endl;
        cout<< "cut_3 " << cuts[0][3] << endl;
        cout<< "cut_4 " << cuts[0][4] << endl; 
        cout << "DETECTABLE" << endl;
        cout<< "cut_0 " << cuts[2][0] << endl;
        cout<< "cut_1 " << cuts[2][1] << endl;
        cout<< "cut_2 " << cuts[2][2] << endl;
        cout<< "cut_3 " << cuts[2][3] << endl;
        cout<< "cut_4 " << cuts[2][4] << endl;
        cout << "DETECTED" << endl;
        cout<< "cut_0 " << cuts[1][0] << endl;
        cout<< "cut_1 " << cuts[1][1] << endl;
        cout<< "cut_2 " << cuts[1][2] << endl;
        cout<< "cut_3 " << cuts[1][3] << endl;
        cout<< "cut_4 " << cuts[1][4] << endl;
    }
    _nEvt ++ ;
}



void SusyCutflow::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void SusyCutflow::end(){ 

    _rootfile->Write();
}
