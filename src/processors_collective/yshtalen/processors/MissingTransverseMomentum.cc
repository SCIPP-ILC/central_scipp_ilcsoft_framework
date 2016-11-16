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
 * authored by Jane Shtalenkova
 * August 5, 2016
 * Prediction Algorithm
 * Uses the sum of the hadronic system particles to predict the high energy e+/e- vector
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
static TH1F* _thetaDiff;
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

    _rootfile = new TFile("predict_test.root","RECREATE");
  //  _hitmap_hi = new TH2F("high_hit","Hit Distribution",600.0,-300.0,300.0,600.0,-300.0,300.0);
  //  _hitmap_had = new TH2F("hadronic_hit","Hit Distribution",600.0,-300.0,300.0,600.0,-300.0,300.0);
    _scalar = new TH1F("scalar", "Transverse Momentum Scalar Magnitude", 2000.0, 0.0, 20.0);
    _vector = new TH1F("vector", "Transverse Momentum Vector Magnitude", 2000.0, 0.0, 20.0);
    _mass = new TH1F("mass", "Mass Parameter", 2000.0, 0.0, 20.0);
    _thetaDiff = new TH1F("theta", "Difference in Angle Between Predicting and Deflected Lepton Vectors", 2000.0, 0.0, 0.1);
    
    // usually a good idea to
    //printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;
    _eventMax = 0;

    _no_def_count = 0;
    _e_def_count = 0;
    _p_def_count = 0;
    _b_def_count = 0;

    _correct = 0;
    _false_pos = 0;
    _false_neg = 0;
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

    int e_index;
    int p_index;

    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
        cout << "this collection has " << nElements << " elements" << endl;
        cout << endl; 
        cout << endl; 
        cout << endl; 


//----------------------------------------HIGH-ENERGY-------------------------------------------------
        //first, find an electron and positron in the event
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
    
           id = hit->getPDG(); 
           stat = hit->getGeneratorStatus();
           
           if(stat==1){
                if(id==11){
                    high_e = hit;
                    e_index = hitIndex;
                    cout << "PRIMARY ELECTRON INDEX: " << e_index << endl;
                }
                if(id==-11){
                    high_p = hit;
                    p_index = hitIndex;
                    cout << "PRIMARY POSITRON INDEX: " << p_index << endl;
                }
           }//end final state
        }//end for loop
        

        
        //create high energy vector, track deflections
        mom_e = high_e->getMomentum();
        cout << "PRIMARY ELECTRON at index " << e_index << ": [" << mom_e[0] << ", " << mom_e[1] << ", " << mom_e[2] << "]" << endl;
        mom_p = high_p->getMomentum();
        cout << "PRIMARY POSITRON at index " << p_index << ": [" << mom_p[0] << ", " << mom_p[1] << ", " << mom_p[2] << "]" << endl;

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

        //both particles deflected
        if(b_def==1){
            cout << "Both particles deflected" << endl;
            _b_def_count++;
            high_vec[0]=mom_e[0]+mom_p[0];
            high_vec[1]=mom_e[1]+mom_p[1];
            high_vec[2]=mom_e[2]+mom_p[2];
        }
        
        else{
            //electron deflection
            if(e_def==1){
                cout << "Electron deflected" << endl;
                _e_def_count++;
                high_vec[0]=mom_e[0];
                high_vec[1]=mom_e[1];
                high_vec[2]=mom_e[2];
            }
            //positron deflection
            else if(p_def==1){
                cout << "Positron deflected" << endl;
                _p_def_count++;
                high_vec[0]=mom_p[0];
                high_vec[1]=mom_p[1];
                high_vec[2]=mom_p[2];
            }
            //no deflection
            else{
                cout << "No deflection" << endl;
                _no_def_count++;
            } 
        }



//------------------------------------------HADRONIC-----------------------------------------------
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
               

/******************  ***  LORENTZ TRANSFORM ***  ***********************************/
                //collect parameters necessary for Lorentz transform
                double in_x = mom[0];
                double in_energy = hit->getEnergy();

                
                //apply transform
                scipp_ilc::transform_to_lab(in_x, in_energy, out_x, out_energy);
                
                //adjust x position
                pos[0] = out_x*pos[2]/mom[2];
/***********************************************************************************/
                
                //shift origin to center of beamcal beampipe hole
                scipp_ilc::z_to_beam_out(pos[0], pos[1], pos[2]);

                double rad = sqrt(pow(pos[0], 2)+ pow(pos[1], 2));
                
                //include hadronic only
                if(hitIndex!=e_index && hitIndex!=p_index){
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

        //--------------------------PLOTTING------------------------------------------------------------------
        if(_nEvt<100000){
            //vector magnitudes and theta calculation
            double mag_scatter = sqrt(pow(scatter_vec[0], 2)+pow(scatter_vec[1], 2)+pow(scatter_vec[2], 2));
            double mag_high = sqrt(pow(high_vec[0], 2)+pow(high_vec[1], 2)+pow(high_vec[2], 2));
            double dot = scatter_vec[0]*high_vec[0]+scatter_vec[1]*high_vec[1]+scatter_vec[2]*high_vec[2];
            double theta = acos(dot/(mag_scatter*mag_high));

            _thetaDiff->Fill(theta);

            //scatter composite and high energy particles position on BeamCal face
            double scatter_pos_x = scatter_vec[0]*pos[2]/scatter_vec[2];
            double scatter_pos_y = scatter_vec[1]*pos[2]/scatter_vec[2];
            double high_pos_x = high_vec[0]*pos[2]/high_vec[2];
            double high_pos_y = high_vec[1]*pos[2]/high_vec[2];
            

/**************************   DETERMINE HIT STATUS   **********************************/
            //determine prediction and high energy vector hit status
            int hit_scatter = scipp_ilc::get_hitStatus(scatter_pos_x, scatter_pos_y);
            int hit_high = scipp_ilc::get_hitStatus(high_pos_x, high_pos_y);

            //determine vector hit status
            //exclude events outside the Beamcal
            if(e_def==1 || b_def==1 || p_def==1){
                if(hit_scatter!=2&&hit_high!=2){
                    _eventMax++;
                    if(hit_scatter==hit_high){
                        _correct++;
                        cout << "CORRECT" << endl;
                    }
                    else{
                        if(hit_scatter==1&&hit_high!=1){
                            _false_pos++;
                            cout << "FALSE POSITIVE" << endl;
                            }
                        if(hit_scatter!=1&&hit_high==1){
                            _false_neg++; 
                            cout << "FALSE NEGATIVE" << endl;
                            }
                    }
                }
            }


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
    cout << "Events in Beamcal radius: " << _eventMax << endl;
    double correct = (double)_correct / (double)_eventMax;
    cout << "Ratio correctly identified: " << correct  << endl; 
    double positive =  (double)_false_pos / (double)_eventMax;
    cout << "Ratio false positive: " << positive << endl; 
    double negative = (double)_false_neg / (double)_eventMax;
    cout << "Ratio false negative: " << negative << endl; 
    _rootfile->Write();
}

