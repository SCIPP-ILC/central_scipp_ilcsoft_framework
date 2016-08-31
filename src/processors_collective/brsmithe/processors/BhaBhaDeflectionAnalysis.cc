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

#include "BhaBhaDeflectionAnalysis.h"
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

BhaBhaDeflectionAnalysis BhaBhaDeflectionAnalysis;

static TFile* _rootfile;
//static TH2F* _hitmap;
static TH1F* _mass;
//static TH1F* _scalar;
static TH1F* _vector;

static TH1F* _xSum;
static TH1F* _ySum;

BhaBhaDeflectionAnalysis::BhaBhaDeflectionAnalysis() : Processor("BhaBhaDeflectionAnalysis") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
}



void BhaBhaDeflectionAnalysis::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("eBpW_dump.root","RECREATE");
    _vector = new TH1F("vector", "Deflected Particle Momentum Magnitude, sqrt(pX^2+pY^2)", 2000.0, 0.0, 20.0);
    _mass = new TH1F("mass", "Deflected Particle sqrt(Q^2) = sqrt(E^2 - <del_p>^2)", 2000.0, 0.0, 3.0);
    
    // usually a good idea to
    //printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;

}



void BhaBhaDeflectionAnalysis::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 


void BhaBhaDeflectionAnalysis::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...


    LCCollection* col = evt->getCollection( _colName ) ;

    //Lorentz transform parameters
    double theta, out_energy, out_x;
    
    double scatter_vec[] = {0, 0, 0};
    double mag = 0;
    double pos[] = {0, 0, 0};
    double energy = 0;
    double theta;
    int id, stat;
    const bool doTransform = true;

    MCParticle* high_e;
    MCParticle* high_p;

    const double* mom;


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
                //if(id==12 || id==14 || id==16){_neutrino_counter++;}
           }//end final state
        }//end for loop
        
        //create sum vector
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
    
        cout << "event = " << _nEvt << endl;
        
        mom = hit->getMomentum();
        
        const double* mom_e = high_e->getMomentum();
        const double* mom_p = high_p->getMomentum();
        
	
        if(stat==1){
	  
	  if (doTransform){
	    //create position vector by ratios from known z pos and momentum
	    pos[2] = scipp_ilc::_Beamcal_zmin;
	    pos[1] = mom[1]*pos[2]/mom[2];
	    pos[0] = mom[0]*pos[2]/mom[2];
	    
	    //collect parameters for Lorentz transform
	    double in_x = mom[0];
	    double in_energy = hit->getEnergy();
	    
	    //apply the transform
	    scipp_ilc::transform_to_lab(in_x, in_energy, out_x, out_energy);

	    //adjust the x position
	    pos[0] = out_x*pos[2]/mom[2];
		
	    //shift origin to the center of the beamcal beampipe hole
	    scipp_ilc::z_to_beam_out(pos[0], pos[1], pos[2]);
	    
	  }
	  else{
	    out_x = mom[0];
	    out_energy = hit->getEnergy();
	  }

          if(hit!=high_e && hit!=high_p){
	    scatter_vec[0]+=out_x;    
	    scatter_vec[1]+=mom[1];    
	    scatter_vec[2]+=mom[2];    
	  }
	  
	  //adjust the energy to post transform
	  energy+=out_energy;

	  //find the new angle
	  double tmag = sqrt(pow(out_x, 2)+pow(mom[1], 2));
	  mag+=tmag;
	  theta = atan(tmag/abs(mom[2]));
	  
	  //Debuggin'
	  cout << "The particle is at (" << pos[0] << "," << pos[1] << "," << pos[2] << ")." << endl;

        }
        }
        //all
        if(_nEvt<200000){
            
            double q_2 = pow((250.0-energy), 2) - pow(scatter_vec[0], 2) - pow(scatter_vec[1], 2) - pow((250.0-abs(scatter_vec[2])), 2);
            double mass = sqrt(-q_2);
            _mass->Fill(mass);

            //fill vector
            double vector = sqrt(pow(scatter_vec[0], 2) + pow(scatter_vec[1], 2));
            _vector->Fill(vector);
                
        }
    }//end collection
    _nEvt ++ ;
}//end process



void BhaBhaDeflectionAnalysis::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void BhaBhaDeflectionAnalysis::end(){
    _rootfile->Write();
    cout << "That's all folks!" << endl;
}

