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

#include "Basic.h"
#include "scipp_ilc_utilities.h"
#include "scipp_ilc_globals.h"
#include "polar_coords.h"
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

Basic Basic;

static TFile* _rootfile;
static TH2F* _hitmap;
static TH1F* _mass;
static TH1F* _scalar;
static TH1F* _vector;
static TH1F* _neutrinos;

Basic::Basic() : Processor("Basic") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
}



void Basic::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("vec_eBpB.root","RECREATE");
    _hitmap = new TH2F("hitmap","Hit Distribution",600.0,-300.0,300.0,600.0,-300.0,300.0);
    _scalar = new TH1F("scalar", "Transverse Momentum Scalar Magnitude", 2000.0, 0.0, 20.0);
    _vector = new TH1F("vector", "Transverse Momentum Vector Magnitude", 2000.0, 0.0, 20.0);
    _mass = new TH1F("mass", "Mass Parameter", 2000.0, 0.0, 20.0);
    // usually a good idea to
    //printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;

}



void Basic::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void Basic::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...


    LCCollection* col = evt->getCollection( _colName ) ;
    cout << endl;
    cout << endl;
    cout << endl;
    cout << "event = " << _nEvt << endl;
    
    double scatter_vec[] = {0, 0, 0};
    double mag = 0;
    double pos[] = {0, 0, 0};
    double energy = 0;
    double theta;
    int id, stat;

    MCParticle* high_e;
    MCParticle* high_p;


    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
        
        //first, find an electron and positron in the event
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
    
           id = hit->getPDG(); 
           stat = hit->getGeneratorStatus();
           
           if(stat==1){
                if(id==11){
                    high_e = hit;
                }
                if(id==-11){
                    high_p = hit;
                }
           }//end final state
        }//end for loop
        
        //------------------------HIGH-ENERGY ANALYSIS-----------------------------------------
        //determine if these are the high energy electron and positron
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
    
           id = hit->getPDG(); 
           stat = hit->getGeneratorStatus();
           
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
        
        //-------------------------HADRONIC SYSTEM-----------------------------------------------
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
           
           id = hit->getPDG(); 
           stat = hit->getGeneratorStatus();
           
           if(stat==1){


                //determine stdhep position
                const double* mom = hit->getMomentum();
                
                //include hadronic only
                if(hit!=high_e && hit!=high_p){
                    //exclude neutrinos
                   // if(id!=12 && id!=14 && id!=16){        
                        if(abs(mom[0])>0.0){
                            scatter_vec[0]+=mom[0];               
                        }
                        if(abs(mom[1])>0.0){
                            scatter_vec[1]+=mom[1];               
                        }
                        if(abs(mom[2])>0.0){
                            scatter_vec[2]+=mom[2];               
                        }
                        
                        energy+=hit->getEnergy();
                        
                        double tmag = sqrt(pow(mom[0], 2)+pow(mom[1], 2));
                        mag+=tmag;  

                        theta = atan(tmag/abs(mom[2]));
                    //}
                    //else{neutrino_counter++;} 

                } 
           }//end final state
        }//end for

        //all
        if(_nEvt<1600000){
                double mass = sqrt(pow(energy, 2)-pow(scatter_vec[0], 2)-pow(scatter_vec[1], 2)-pow(scatter_vec[2], 2));
                _mass->Fill(mass);
                cout << "Mass parameter: " << mass << endl;

                //fill scalar 
                _scalar->Fill(mag);
                cout << "Scalar Momentum: " << mag << endl;

                //fill vector
                double vector = sqrt(pow(scatter_vec[0], 2) + pow(scatter_vec[1], 2));
                _vector->Fill(vector);
                cout << "Vector Momentum: " << vector << endl;
                  
                
        }
    }
    _nEvt ++ ;
}



void Basic::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void Basic::end(){ 

    _rootfile->Write();
}
