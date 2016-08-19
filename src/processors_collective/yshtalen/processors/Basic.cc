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

    _rootfile = new TFile("eBpW_all.root","RECREATE");
    //_hitmap = new TH2F("hitmap","Hit Distribution",600.0,-300.0,300.0,600.0,-300.0,300.0);
    //_scalar = new TH1F("scalar", "Transverse Momentum Scalar Magnitude", 2000.0, 0.0, 20.0);
    _vector = new TH1F("vector", "Deflected Particle Momentum Magnitude, sqrt(pX^2+pY^2)", 2000.0, 0.0, 20.0);
    _mass = new TH1F("mass", "Deflected Particle sqrt(Q^2) = sqrt(E^2 - <del_p>^2)", 2000.0, 0.0, 3.0);
    
    _xSum = new TH1F("xSum","X-Momentum Event Total",600.0,-10.0,10.0);
    _ySum = new TH1F("ySum","Y-Momentum Event Total",600.0,-10.0,10.0);
    
    // usually a good idea to
    //printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;

    _neutrino_counter=0;
    _p_def_count=0;
    _e_def_count=0;
    _b_def_count=0;
    _zero_scatter_count=0;
    _low_scatter_count=0;
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
    double energy = 0;
    double theta;
    int id, stat;

    MCParticle* high_e;
    MCParticle* high_p;

    int p_def=0;
    int e_def=0;
    int b_def=0;

    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
        
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
        

        //------------------------HIGH-ENERGY ANALYSIS-----------------------------------------
        //determine if these are the high energy electron and positron
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
    
           id = hit->getPDG();
           stat = hit->getGeneratorStatus();

           const double* mom = hit->getMomentum();
           if(stat==1){
               cout << id << endl;;
               cout << endl;
               cout << endl;
               cout << endl;
                //if(hit!=high_e && hit!=high_p){
                    scatter_vec[0]+=mom[0];
                    scatter_vec[1]+=mom[1];
                    scatter_vec[2]+=mom[2];
                    energy+=hit->getEnergy();
                //}
           }//end final state
        }//end for loop
        //-------------------------HADRONIC SYSTEM-----------------------------------------------
        /*for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
           
           id = hit->getPDG(); 
           stat = hit->getGeneratorStatus();
           
           if(stat==1){


                //determine stdhep position
                const double* mom = hit->getMomentum();
                
                //include hadronic only
                if(hit!=high_e && hit!=high_p){
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
                }
                else{
                    if(id==11){
                        cout << "Electron mom: " << mom[0] << " " << mom[1] << " " << mom[2] << endl;
                        if(abs(mom[0])!=0 || abs(mom[1])!=0){
                            cout << "Found deflected electron" << endl;
                            e_def=1;    
                        }        
                    } 
                    else{
                        cout << "Positron mom: " << mom[0] << " " << mom[1] << " " << mom[2] << endl;
                        if(abs(mom[0])!=0 || abs(mom[1])!=0){
                            cout << "Found deflected positron" << endl;
                            p_def=1;
                            if(e_def==1){b_def=1; e_def=0; p_def =0;}    
                        }        
                    } 
                       
                } 
           }//end final state
        }//end for
        */
        const double* mom_e = high_e->getMomentum();
        const double* mom_p = high_p->getMomentum();

        //all
        if(_nEvt<200000){
            
            _xSum->Fill(scatter_vec[0]);
            _ySum->Fill(scatter_vec[1]);

            if(e_def==1 && p_def==1){_b_def_count++;}
            else if(e_def==1){_e_def_count++;}
            else if(p_def==1){_p_def_count++;}
            
            //double mass = sqrt(pow(energy, 2)-pow(scatter_vec[0], 2)-pow(scatter_vec[1], 2)-pow(scatter_vec[2], 2));
            double q_2 = pow((250.0-energy), 2) - pow(scatter_vec[0], 2) - pow(scatter_vec[1], 2) - pow((250.0-abs(scatter_vec[2])), 2);
            double mass = sqrt(-q_2);
            _mass->Fill(mass);
    //        cout << "Mass parameter: " << mass << endl;

            //fill scalar 
            //_scalar->Fill(mag);
            //cout << "Scalar Momentum: " << mag << endl;

            //fill vector
            double vector = sqrt(pow(scatter_vec[0], 2) + pow(scatter_vec[1], 2));
            if(vector==0.0){_zero_scatter_count++;}
            if(vector>0.0 && vector< 0.01){_low_scatter_count++;}
            _vector->Fill(vector);
      //      cout << "Vector Momentum: " << vector << endl;
                  
                
        }
    }
    _nEvt ++ ;
}



void Basic::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void Basic::end(){
    /*cout << "Number of NEUTRINOS: " << _neutrino_counter << endl;
    cout << "Number deflected ELECTRONS: " << _e_def_count << endl;
    cout << "Number deflected POSITRONS: " << _p_def_count << endl;
    cout << endl; 
    cout << endl; 
    cout << endl; 
    
    double e_def_ratio= _e_def_count/10000.0;
    double p_def_ratio= _p_def_count/10000.0;
    double b_def_ratio= _b_def_count/10000.0;
    double zero_scatter_ratio = _zero_scatter_count/10000.0;
    double low_scatter_ratio = _low_scatter_count/10000.0;


    cout << "            DEFLECTION RATIOS" << endl;
    cout << endl; 
    cout << "Ratio deflected electrons: " << e_def_ratio  << endl;
    cout << "Ratio deflected positrons: " << p_def_ratio << endl;
    cout << "Ratio deflected both: " << b_def_ratio << endl;
    cout << endl; 
    cout << endl; 
    cout << endl; 


    cout << "            SCATTER RATIO" << endl;
    cout << endl; 
    cout << "Ratio events where scatter sums to zero: " << zero_scatter_ratio << endl;
    cout << "Ratio events where 0 < scatter_sum < 0.01: " << low_scatter_ratio << endl;
    cout << endl; 
    cout << endl; 
    cout << endl; */
    _rootfile->Write();
}

