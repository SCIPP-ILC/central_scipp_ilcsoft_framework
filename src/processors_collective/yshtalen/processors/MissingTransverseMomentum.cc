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

#include "MissingTransverseMomentum.h"
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

MissingTransverseMomentum MissingTransverseMomentum;

static TFile* _rootfile;
static TH2F* _hitmap;
static TH1F* _mass;
static TH1F* _scalar;
static TH1F* _vector;


MissingTransverseMomentum::MissingTransverseMomentum() : Processor("MissingTransverseMomentum") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
}



void MissingTransverseMomentum::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("hitmapeWpB_able.root","RECREATE");
    _hitmap = new TH2F("hitmap","Hit Distribution",300.0,-150.0,150.0,300.0,-150.0,150.0);
    _scalar = new TH1F("scalar", "Transverse Momentum Scalar Magnitude", 2000.0, 0.0, 20.0);
    _vector = new TH1F("vector", "Transverse Momentum Vector Magnitude", 2000.0, 0.0, 20.0);
    _mass = new TH1F("mass", "Mass Parameter", 2000.0, 0.0, 20.0);

    // usually a good idea to
    //printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;

    _no_def_count = 0;
    _e_def_count = 0;
    _p_def_count = 0;
}



void MissingTransverseMomentum::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void MissingTransverseMomentum::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    LCCollection* col = evt->getCollection( _colName ) ;
    cout << endl;
    cout << endl;
    cout << endl;
    cout << "event = " << _nEvt << endl;
    
    double scatter_vec[] = {0, 0, 0};
    double mag = 0;
    double energy = 0;

    MCParticle* high_e;
    MCParticle* high_p;


    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
        
        
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
    
           int id = hit->getPDG(); 
           int stat = hit->getGeneratorStatus();
           
           if(stat==1){
                if(id==11){
                    high_e = hit;
                }
                if(id==-11){
                    high_p = hit;
                }
           }//end final state
        }//end for loop
        
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
    
           int id = hit->getPDG(); 
           int stat = hit->getGeneratorStatus();
           
           if(stat==1){

                //find high energy electron
                if(id==11){
                    if(hit->getEnergy()>high_e->getEnergy()){
                        high_e = hit;
                    }               
                }    
                //find high energy positron
                if(id==-11){
                    if(hit->getEnergy()>high_p->getEnergy()){
                        high_p = hit;
                    }               
                }
                    
           }//end final state
        }//end for loop
        
        //--------------HADRONIC SYSTEM------------------------------------------------------//
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
    
           int id = hit->getPDG(); 
           int stat = hit->getGeneratorStatus();
           
           if(stat==1){
                const double* mom = hit->getMomentum();
                double in_x = mom[0];
                double in_energy = hit->getEnergy();
                double out_energy, out_x;
                
                double theta = 0.007;
                double in_E = 250.0;
                
                double beta = sin(theta);
                double gamma = pow((1-pow(beta, 2)), -0.5);

                 //*
                 //*    |gamma         gamma*beta| |p|                       |p'|
                 //*    |gamma*beta         gamma|*|E| = TRANSOFORMED4vector |E'|
                 //*        
                 //*                                     *                                     *      

                out_x = in_x*gamma + gamma*beta*in_energy;
                out_energy= in_energy*gamma + gamma*beta*in_x;
          

                //include hadronic only
                if(hit!=high_e && hit!=high_p){
                    //exclude neutrinos
                    if(id!=12 && id!=14 && id!=16){        
                        if(abs(out_x)>0.0){
                            scatter_vec[0]+=out_x;               
                        }
                        if(abs(mom[1])>0.0){
                            scatter_vec[1]+=mom[1];               
                        }
                        if(abs(mom[2])>0.0){
                            scatter_vec[2]+=mom[2];               
                        }
                        
                        energy+=out_energy;
                        
                        mag+=sqrt(pow(out_x, 2)+pow(mom[1], 2));  
                    }                  
                } 
           }//end final state
        }//end for

        //all
        if(_nEvt<1600000){
            double mass = sqrt(pow(energy, 2)-pow(scatter_vec[0], 2)-pow(scatter_vec[1], 2)-pow(scatter_vec[2], 2));
            _mass->Fill(mass);
            //cout << "Mass parameter: " << mass << endl;

            //fill scalar 
            _scalar->Fill(mag);
            //cout << "Scalar Momentum: " << mag << endl;

            //fill vector
            double vector = sqrt(pow(scatter_vec[0], 2) + pow(scatter_vec[1], 2));
            _vector->Fill(vector);
            //cout << "Vector Momentum: " << vector << endl;
              
            //detected
            //detectable
            }
        }

    _nEvt ++ ;
}



void MissingTransverseMomentum::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void MissingTransverseMomentum::end(){ 

    _rootfile->Write();
}
