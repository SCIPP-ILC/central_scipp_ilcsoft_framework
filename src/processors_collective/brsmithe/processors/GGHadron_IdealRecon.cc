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

#include "GGHadron_IdealRecon.h"
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


GGHadron_IdealRecon GGHadron_IdealRecon;

static TFile* _rootfile;
static TH2F* _hitmap;
static TH1F* _vector;
static TH1F* _scalar;
static TH1F* _mass;
static TH1F* _theta;
static TH1F* _energy;
static TH1F* _pMom;
static TH1F* _eMom;
static TH1F* _xSum;
static TH1F* _ySum;
static TH1F* _zSum;
static TH1F* _pTot;

GGHadron_IdealRecon::GGHadron_IdealRecon() : Processor("GGHadron_IdealRecon") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
}



void GGHadron_IdealRecon::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("hitmapeBpW.root","RECREATE");
    _hitmap = new TH2F("hitmap","Hit Distribution",300.0,-150.0,150.0,300.0,-150.0,150.0);
    _scalar = new TH1F("scalar", "Transverse Momentum Scalar Magnitude", 2000.0, 0.0, 20.0);
    _vector = new TH1F("vector", "Transverse Momentum Vector Magnitude", 2000.0, 0.0, 20.0);
    _mass = new TH1F("mass", "Mass Parameter", 2000.0, 0.0, 20.0);  
    _theta = new TH1F("theta", "Angle Between Prediction Vector and High Energy Vector", 1000.0, 0.0, 0.01); 
    _energy = new TH1F("energy", "Total Final State Energy, eBpB", 50000.0, 300.0, 505.0); 

    _pMom = new TH1F("pMom", "Positron Momentum Vector Magnitudes, eBpB", 1000.0, 0.0, 4.0);
    _eMom = new TH1F("eMom", "Electron Momentum Vector Magnitudes, eBpB", 1000.0, 0.0, 4.0);
    _xSum = new TH1F("xSum", "X-Momentum Event Total, eBpW", 1000.0, -2.0, 2.0); 
    _ySum = new TH1F("ySum", "Y-Momentum Event Total, eBpW", 1000.0, -2.0, 2.0); 
    _zSum = new TH1F("zSum", "Z-Momentum Event Total", 1000.0, -20.0, 20.0); 
    
    _pTot = new TH1F("pTot", "Transverse Momentum Event Total", 1000.0, 0.0, 2.0); 
    // usually a good idea to
    //printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;

    _no_def_count = 0;
    _e_def_count = 0;
    _p_def_count = 0;
}



void GGHadron_IdealRecon::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void GGHadron_IdealRecon::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    LCCollection* col = evt->getCollection( _colName ) ;
    cout << endl;
    cout << endl;
    cout << endl;
    cout << "event = " << _nEvt << endl;
    
    double scatter_vec[] = {0, 0, 0};
    double high_vec[] = {0, 0, 0};
    double mag = 0;
    double energy = 0;
    double e_scatter = 0;

    int no_def = 0;
     
    const double* mom_e;
    const double* mom_p;

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

        mom_e =  high_e->getMomentum();
        mom_p =  high_p->getMomentum();
        cout << "high energy electron momentum: [" << mom_e[0] << ", " << mom_e[1] << ", " << mom_e[2] << "]" << endl;
        cout << "high energy positron momentum: [" << mom_p[0] << ", " << mom_p[1] << ", " << mom_p[2] << "]" << endl;
        cout << endl;

        if(mom_e[0]!=0 || mom_e[1]!=0){
            high_vec[0] = mom_e[0];
            high_vec[1] = mom_e[1];
            high_vec[2] = mom_e[2];
        }
        else if(mom_p[0]!=0 || mom_p[1]!=0){
            high_vec[0] = mom_p[0];
            high_vec[1] = mom_p[1];
            high_vec[2] = mom_p[2];
        }
        else{
            no_def=1;
        }    
            
       //create high energy transverse momentum vectors 
       double e_mom = sqrt(pow(mom_e[0], 2)+pow(mom_e[1], 2));
       double p_mom = sqrt(pow(mom_p[0], 2)+pow(mom_p[1], 2));
       _pMom->Fill(p_mom);
       _eMom->Fill(e_mom);

        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
          
           //only check final state 
           int stat = hit->getGeneratorStatus();
           if(stat==1){
            
               //get endpoint and fill hitmap
               const double* pos = hit->getEndpoint();
               _hitmap->Fill(pos[0],pos[1]);
              cout << pos[0] << " " << pos[1] << " " << pos[2] << endl; 
                //add to event total energy
                energy+= hit->getEnergy();

                //add to momentum vector component-wise for every final state scatter particle
                const double* mom = hit->getMomentum();
                int id = hit->getPDG();
                if(hit!=high_e && hit!=high_p){
                    //cout << "id: " << id << "     " << "momentum: [" << mom[0] << ", " << mom[1] << ", " << mom[2] << "]" << endl;

                    //add particle mom magnitude and energy to event totals
                    mag += sqrt(pow(mom[0], 2)+pow(mom[1], 2));
                    e_scatter+=hit->getEnergy();

                    if(abs(mom[0]) > 0.0){
                        scatter_vec[0] += (double) mom[0];
                    }
                    if(abs(mom[1]) > 0.0){
                        scatter_vec[1] += (double) mom[1];
                    }
                    if(abs(mom[2]) > 0.0){
                        scatter_vec[2] += (double) mom[2];
                    }    
                }
            }//end final state
        }//end for loop
   

        
    }

//-----------------------------HADRONIC SYSTEM----------------------------------------------------------------------
    if(_nEvt<1600000){
        //fill mass
        double mass = sqrt(pow(e_scatter, 2)-pow(scatter_vec[0], 2)-pow(scatter_vec[1], 2)-pow(scatter_vec[2], 2));
        _mass->Fill(mass);
        //cout << "Mass parameter: " << mass << endl;
        
        //fill scalar 
        _scalar->Fill(mag);
        //cout << "Scalar Momentum: " << mag << endl;
        
        //fill vector
        double vector = sqrt(pow(scatter_vec[0], 2) + pow(scatter_vec[1], 2));
        _vector->Fill(vector);
        //cout << "Vector Momentum: " << vector << endl;
        //
    }
//-------------------------------------------------------------------------------------------------------------------
    if(no_def==0){
        //reflect vector through the origin
        scatter_vec[0] = -scatter_vec[0];
        scatter_vec[1] = -scatter_vec[1];
        scatter_vec[2] = -scatter_vec[2];

        //calculate longitudinal momentum prediction
        scatter_vec[2] = sqrt(pow(250.0 - e_scatter, 2) - pow(mag, 2)); 
        
        double smag = sqrt(pow(scatter_vec[0], 2)+pow(scatter_vec[1], 2)+pow(scatter_vec[2], 2));
        double hmag = sqrt(pow(high_vec[0], 2)+pow(high_vec[1], 2)+pow(high_vec[2], 2));

        double cos = (scatter_vec[0]*high_vec[0] + scatter_vec[1]*high_vec[1] + scatter_vec[2]*high_vec[2])/(smag*hmag);
        double theta = acos(cos);

        cout << "Prediction Vector: " << "[" << scatter_vec[0] << ", " << scatter_vec[1] << ", " << scatter_vec[2] << "]" << endl;
        cout << "High Energy Vector: " << "[" << high_vec[0] << ", " << high_vec[1] << ", " << high_vec[2] << "]" << endl;
        cout << "Error Angle: " << theta << endl;
        
        if(high_vec[2]>0){
            _e_def_count++;
        }
        else{
            _p_def_count++;
        }

        _theta->Fill(theta);
    }
    else{
        _no_def_count++;
    }    

//--------------------------------------TOTAL SYSTEM---------------------------------------------------------------------
    scatter_vec[0]= scatter_vec[0] + mom_e[0] + mom_p[0];
    scatter_vec[1]= scatter_vec[1] + mom_e[1] + mom_p[1];
    scatter_vec[2]= scatter_vec[2] + mom_e[2] + mom_p[2];

    //fill energy
    _energy->Fill(energy);
    double total_trans = sqrt(pow(scatter_vec[0], 2) + pow(scatter_vec[1], 2));
    _pTot->Fill(total_trans);

    //fill x and y mom check sums
    _xSum->Fill(scatter_vec[0]);
    _ySum->Fill(scatter_vec[1]);
    _zSum->Fill(scatter_vec[2]);
    _nEvt ++ ;
}



void GGHadron_IdealRecon::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void GGHadron_IdealRecon::end(){ 

    _rootfile->Write();
    cout << "Ratio of events with no deflection: " << _no_def_count/_nEvt << endl;
    cout << "Electron deflected: " << _e_def_count/_nEvt << endl;
    cout << "Positron deflected: " << _p_def_count/_nEvt << endl;
}
