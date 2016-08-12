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

MissingTransverseMomentum MissingTransverseMomentum;

static TFile* _rootfile;
static TH2F* _hitmap;
static TH2F* _hitmap_Lorentz;
static TH2F* _hitmap_Lorentz_shift;
static TH1F* _mass;
static TH1F* _scalar;
static TH1F* _vector;
static TH1F* _neutrinos;

MissingTransverseMomentum::MissingTransverseMomentum() : Processor("MissingTransverseMomentum") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
}



void MissingTransverseMomentum::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("pos_eBpW.root","RECREATE");
    _hitmap = new TH2F("hitmap","Hit Distribution",600.0,-300.0,300.0,600.0,-300.0,300.0);
    _hitmap_Lorentz = new TH2F("hitmap_Lorentz","Hit Distribution",600.0,-300.0,300.0,600.0,-300.0,300.0);
    _hitmap_Lorentz_shift = new TH2F("hitmap_Lorentz_shift","Hit Distribution",600.0,-300.0,300.0,600.0,-300.0,300.0);
    _scalar = new TH1F("scalar", "Transverse Momentum Scalar Magnitude", 2000.0, 0.0, 20.0);
    _vector = new TH1F("vector", "Transverse Momentum Vector Magnitude", 2000.0, 0.0, 20.0);
    _mass = new TH1F("mass", "Mass Parameter", 2000.0, 0.0, 20.0);
    _neutrinos = new TH1F("neutrinos", "Neutrinos per Event", 10.0, 0.0,10.0); 
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
    double pos[] = {0, 0, 0};
    double energy = 0;
    double theta;
    int neutrino_counter=0;
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
        //then, determine if these are the high energy electron and positron
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
                cout << "Momentum stdhep: " << mom[0] << ", " << mom[1] << ", " << mom[2] << endl; 
                
                //create position vector by ratios from known z position and momentum
                pos[2] = scipp_ilc::_BeamCal_zmin;
                pos[0] = mom[0]*pos[2]/mom[2];
                pos[1] = mom[1]*pos[2]/mom[2];
                
                cout << "Position from stdhep: " << pos[0] << ", " << pos[1] << ", " << pos[2] << endl; 
                //fill stdhep hitmap
                double rad = sqrt(pow(pos[0], 2)+ pow(pos[1], 2));
                if(rad<140){
                    _hitmap->Fill(pos[0], pos[1]);
                }
                //collect parameters necessary for Lorentz transform

                double in_x = mom[0];
                double in_energy = hit->getEnergy();
                double out_energy, out_x;
                
                //apply transform
                scipp_ilc::transform_to_lab(in_x, in_energy, out_x, out_energy);
                cout << "Momentum after Lorentz transform: " << mom[0] << ", " << mom[1] << ", " << mom[2] << endl; 
                


                pos[0] = out_x*pos[2]/mom[2];
           
                cout << "Position after Lorentz: " << pos[0] << ", " << pos[1] << ", " << pos[2] << endl; 
                rad = sqrt(pow(pos[0], 2) + pow(pos[1], 2));
                if(rad<140){
                    _hitmap_Lorentz->Fill(pos[0], pos[1]);
                }
                scipp_ilc::z_to_beam_out(pos[0], pos[1], pos[2]);

                cout << "Position after coordinate change " << pos[0] << ", " << pos[1] << ", " << pos[2] << endl; 
                
                rad = sqrt(pow(pos[0], 2) + pow(pos[1], 2));
                if(rad<140){
                    _hitmap_Lorentz_shift->Fill(pos[0], pos[1]);
                }
                //include hadronic only
                if(hit!=high_e && hit!=high_p){
                    //exclude neutrinos
                   // if(id!=12 && id!=14 && id!=16){        
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
                        
                        double tmag = sqrt(pow(out_x, 2)+pow(mom[1], 2));
                        mag+=tmag;  

                        theta = atan(tmag/abs(mom[2]));
                    //}
                    //else{neutrino_counter++;} 

                    if(id==12 || id==14 || id==16){
                        neutrino_counter++;
                    }    
                } 
           }//end final state
        }//end for

        //all
        if(_nEvt<1600000){
            if(cos(theta)<0.9){
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
            _neutrinos->Fill(neutrino_counter);
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
