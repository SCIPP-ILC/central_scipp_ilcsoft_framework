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
//static TH2F* _hitmap_hi;
//static TH2F* _hitmap_had;
static TH1F* _deltaX;
static TH1F* _deltaY;
static TH1F* _deltaZ;
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

    _rootfile = new TFile("predicting_eWpW.root","RECREATE");
  //  _hitmap_hi = new TH2F("high_hit","Hit Distribution",600.0,-300.0,300.0,600.0,-300.0,300.0);
  //  _hitmap_had = new TH2F("hadronic_hit","Hit Distribution",600.0,-300.0,300.0,600.0,-300.0,300.0);
    _deltaX = new TH1F("x_diff", "Delta X-Pos, Predicting + High Energy Vectors", 2000.0, -10.0, 10.0);
    _deltaY = new TH1F("y_diff", "Delta Y-Pos, Predicting + High Energy Vectors", 2000.0, -10.0, 10.0);
    _deltaZ = new TH1F("z_diff", "Delta Z-Pos, Predicting + High Energy Vectors", 2000.0, -10.0, 10.0);
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
    _b_def_count = 0;
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
    
    //vector construction variables
    double scatter_vec[] = {0, 0, 0};
    double high_vec[] = {0, 0, 0};
    double mag = 0;
    double pos[] = {0, 0, 0};
    double energy = 0;

    //Lorentz transform parameters
    double theta, out_energy, out_x;

    //MCParticle identifiers
    int id, stat;

    //hit status with respect to the Beamcal for high energy particle and prediction vectors
    // 1 - hit BeamCal
    // 2 - out of Beamcal radius range
    // 3 - down beampipe hole
    int h_hit_status=0;
    int p_hit_status=0;

    //booleans
    int e_def=0;
    int p_def=0;
    int b_def=0;

    //high energy electron and positron objects
    MCParticle* high_e;
    MCParticle* high_p;

    const double* mom_e;
    const double* mom_p;

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
        
        //------------------------HIGH-ENERGY-------------------------------------------------
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
        
        //create high energy vector, track deflections
        mom_e = high_e->getMomentum();
        mom_p = high_p->getMomentum();
        cout << "[" << mom_e[0] << ", " << mom_e[1] << ", " << mom_e[2] << "]" << endl;
        cout << "[" << mom_p[0] << ", " << mom_p[1] << ", " << mom_p[2] << "]" << endl;

        e_def=0;
        p_def=0;
        b_def=0;
        
        if(abs(mom_e[0])!=0 || abs(mom_e[1])!=0){
                e_def=1; 
        }
        if(abs(mom_p[0])!=0 || abs(mom_p[1])!=0){
                p_def=1;       
        }
        if(e_def==1 && p_def==1){
            e_def=0;
            p_def=0;
            b_def=1;
        }
        if(b_def==1){
            _b_def_count++;
            high_vec[0]=mom_e[0]+mom_p[0];
            high_vec[1]=mom_e[1]+mom_p[1];
            high_vec[2]=mom_e[2]+mom_p[2];
        }
        else{
            if(e_def==1){
                _e_def_count++;
                high_vec[0]=mom_e[0];
                high_vec[1]=mom_e[1];
                high_vec[2]=mom_e[2];
            }
            else if(p_def==1){
                _p_def_count++;
                high_vec[0]=mom_p[0];
                high_vec[1]=mom_p[1];
                high_vec[2]=mom_p[2];
            }
            else{_no_def_count++;} 
        }



        //-------------------------HADRONIC-----------------------------------------------
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
           
           id = hit->getPDG(); 
           stat = hit->getGeneratorStatus();
           
           if(stat==1){


                //determine stdhep position
                const double* mom = hit->getMomentum();
                
                //create position vector by ratios from known z position and momentum
                pos[2] = scipp_ilc::_BeamCal_zmin;
                pos[0] = mom[0]*pos[2]/mom[2];
                pos[1] = mom[1]*pos[2]/mom[2];
               


                //collect parameters necessary for Lorentz transform
                double in_x = mom[0];
                double in_energy = hit->getEnergy();

                
                //apply transform
                scipp_ilc::transform_to_lab(in_x, in_energy, out_x, out_energy);
                
                //adjust x position
                pos[0] = out_x*pos[2]/mom[2];
                
                //shift origin to center of beamcal beampipe hole
                scipp_ilc::z_to_beam_out(pos[0], pos[1], pos[2]);

                double rad = sqrt(pow(pos[0], 2)+ pow(pos[1], 2));
                
                //include hadronic only
                if(hit!=high_e && hit!=high_p){
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

                    //track neutrinos
                  //  if(id==12 || id==14 || id==16){
                  //      neutrino_counter++;
                  //  }    
                }
           }//end final state
        }//end for

        //create prediction vector
        scatter_vec[0] = -scatter_vec[0];
        scatter_vec[1] = -scatter_vec[1];
        scatter_vec[2] = sqrt(pow(250.0-energy, 2) - pow(mag, 2));

        //fill delta maps

        //determine prediction vector hit status

        //determine deflected high energy particle hit status

       

        //--------------------------PLOTTING------------------------------------------------------------------
        if(_nEvt<100000){
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
                cout << endl;
                cout << endl;
                cout << endl;
                cout << "Prediction Vector: [" << scatter_vec[0] << ", " << scatter_vec[1] << ", " << scatter_vec[2] << "]"  << endl;
                cout << "High Energy Vector: [" << high_vec[0] << ", " << high_vec[1] << ", " << high_vec[2] << "]"  << endl;
                
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
