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
 * author Jane Shtalenkovae
 * August 5, 2016
 */

#include "Basic.h"
#include "scipp_ilc_utilities.h"
#include "scipp_ilc_globals.h"
#include "polar_coords.h"
#include <iostream>
#include <cmath>

#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>

#include <TFile.h>
#include <TH2D.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"



using namespace lcio;
using namespace marlin;
using namespace std;

Basic Basic;

static TFile* _rootfile;
//static TH2F* _hitmap;
static TH1F* _mass;
//static TH1F* _scalar;
static TH1F* _vector;

static TH1F* _xSum;
static TH1F* _ySum;

Basic::Basic() : Processor("Basic") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
}



void Basic::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("eBpB_check.root","RECREATE");
    _vector = new TH1F("vector", "Scatter Transverse Momentum Vector Magnitude", 2000.0, 0.0, 20.0);
    _mass = new TH1F("mass", "Deflected Particle sqrt(Q^2) = sqrt(E^2 - <del_p>^2)", 2000.0, 0.0, 3.0);
    
    // usually a good idea to
    //printParameters() ;

    _neutrino_counter = 0;
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

    
    double scatter_vec[] = {0, 0, 0};
    double mag = 0;
    double energy = 0;
    double theta;
    int id, stat;

    MCParticle* high_e;
    MCParticle* high_p;

    const double* mom;


    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
        cout << nElements << endl;
         
        //first, find last electron and positron in the event
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
                //find neutrinos 
                if(id==12 || id==14 || id==16){_neutrino_counter++;}
           }//end final state
        }//end for loop
       
        //find last electron and positron in the event
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
    
           id = hit->getPDG(); 
           stat = hit->getGeneratorStatus();
           
           if(stat==1){
                if(id==11){
                    if(hit->getEnergy()>high_e->getEnergy()){
                        high_e = hit;
                    }
                }
                if(id==-11){
                    if(hit->getEnergy()>high_p->getEnergy()){
                        high_p = hit;
                    }
                }
                //find neutrinos 
                //if(id==12 || id==14 || id==16){_neutrino_counter++;}
           }//end final state
        }//end for loop
        
        cout << "event = " << _nEvt << endl;
        //create sum vector
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
            MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
            

            mom = hit->getMomentum();
            
            const double* mom_e = high_e->getMomentum();
            const double* mom_p = high_p->getMomentum();
            
            double e_pt = sqrt(pow(mom_e[0], 2)+ pow(mom_e[1], 2));
            double p_pt = sqrt(pow(mom_p[0], 2)+ pow(mom_p[1], 2));

            double theta_e = atan(e_pt/mom_e[2]);
            double theta_p = atan(p_pt/mom_p[2]);

            if(theta_p > 0.1){
                if(theta_e < 0.01) {
                    cout << "stat: " << stat << " id: " << hit->getPDG() << "     mom: [" << mom[0] << ", " << mom[1] << ", " << mom[2] << "]     energy: " << hit->getEnergy() << endl;
                }   
            }

            //int stat_std = getISTHEP(col->getElementAt(hitIndex));
            //final state excluding high energy electron/positron 
            //if(stat==1){
                //if(hit!=high_e && hit!=high_p){

                    
                    if(abs(mom[0])>0){
                        scatter_vec[0]+=mom[0];
                    }
                    if(abs(mom[1])>0){    
                        scatter_vec[1]+=mom[1];
                    }
                    if(abs(mom[2])>0){
                        scatter_vec[2]+=mom[2];
                    }    
                //}//end electronic system    
           //}//end final state
        }
        //all
        if(_nEvt<1600000){
            
            double q_2 = pow((250.0-energy), 2) - pow(scatter_vec[0], 2) - pow(scatter_vec[1], 2) - pow((250.0-abs(scatter_vec[2])), 2);
            double mass = sqrt(-q_2);
            _mass->Fill(mass);

            //fill vector
            double vector = sqrt(pow(scatter_vec[0], 2) + pow(scatter_vec[1], 2));
            //cout << "VECTOR SUM: " << vector << endl;
            //cout << endl;
            //cout << endl;
            cout << endl;
            _vector->Fill(vector);
                
        }
    }//end collection
    _nEvt ++ ;
}//end process



void Basic::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void Basic::end(){
    _rootfile->Write();
}

