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
static TH2F* _hitmap;
static TH1F* _Zpos;

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

    _rootfile = new TFile("BhaBhaAnalysis.root","RECREATE");
    _vector = new TH1F("vector", "Deflected Particle Momentum Magnitude, sqrt(pX^2+pY^2)", 2000.0, 0.0, 20.0);
    _mass = new TH1F("mass", "Deflected Particle sqrt(Q^2) = sqrt(E^2 - <del_p>^2)", 2000.0, 0.0, 3.0);
    _hitmap = new TH2F("hitmap", "Hitmap of Final State Particles", 2000.0, -400.0, 400.0, 2000.0, -400.0, 400.0);
    _Zpos = new TH1F("Zpos", "Locations on the Z axis", 2000.0, 0.0, 4000.0);

    // usually a good idea to
    //printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;

    _nHitHit = 0;
    _nEHitPMiss = 0;
    _nPHitEMiss = 0;
    _nMissMiss = 0;
    _nTotal = 0;
    
    _BeamCalz = scipp_ilc::_BeamCal_zmin;
    _BeamCalr = scipp_ilc::_BeamCal_outgoing_pipe_radius;
}



void BhaBhaDeflectionAnalysis::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 


void BhaBhaDeflectionAnalysis::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...


    LCCollection* col = evt->getCollection( _colName ) ;

    //Lorentz transform parameters
    double ptheta, pout_energy, pout_x;
    double etheta, eout_energy, eout_x;
    
    double scatter_vec[] = {0, 0, 0};
    double pmag = 0;
    double emag = 0;
    double Epos[] = {0, 0, 0};
    double Ppos[] = {0, 0, 0};
    double Eenergy = 0;
    double Penergy = 0;
    int id, stat;
    const bool doTransform = true;
    
    double ein_x, ein_energy;
    double pin_x, pin_energy;

    const double* mom_e;
    const double* mom_p;

    MCParticle* high_e;
    MCParticle* high_p;
    
    bool Phit_status;
    bool Ehit_status;


    // this will only be entered if the collection is available
    if( col != NULL ){
      
   
        int nElements = col->getNumberOfElements()  ;
        double debug;
	
        //first, find last electron and positron in the event. We'll use these as benchmarks to compare others to.
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
    
           id = hit->getPDG(); 
           stat = hit->getGeneratorStatus();
	   debug = 0;
	   

           if(stat==1){
	     switch(id){
	     case 11:
	       high_e = hit;
	       break;
	     case -11:
	       high_p = hit;
	       break;
	     default: break;
	     }
             
           }//end final state
	   
        }//end for loop
	
	//I declare this here just to make the next group look a little nicer.
	double hitEnergy;
	
	//Now, we loop through the particles again, and at each we compare them to their energy to whatever is the most energetic of its type. 
	//If it is more energetic than the current 'most energetic,' it becomes the new 'most energetic'

	for (int hitIndex=0; hitIndex < nElements ; hitIndex++){
	  MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
	  

	  id = hit->getPDG();
	  stat = hit->getGeneratorStatus();

	  // stat==1 means that this is a final state particle
	  if(stat==1){
	    
	    //Getting this hit's energy to compare it to the next positron/electron hit energy
	    hitEnergy = hit->getEnergy();
	    
	    switch(id){
	    case 11:
	      // If this hit is more energetic than the most energetic, then the hit becomes the most energetic
	      // Most energetic elecron
	      if(hitEnergy > high_e->getEnergy()){ 
		high_e = hit;
	      }
	      break;
	    case -11:
	      // Most energetic Positron
	      if(hitEnergy > high_p->getEnergy()){ 
		high_p = hit;
	      }
	      break;
	    default:
	      break;
	    }// End switch

	  }// End status if
	}// end for loop
        
        
        
	mom_e = high_e->getMomentum();
	mom_p = high_p->getMomentum();
        
	
        if(stat == 1){
	  //create position vector by ratios from known z pos and momentum
	    
	  Epos[2] = _BeamCalz;
	  Epos[1] = mom_e[1]*Epos[2]/mom_e[2];
	  Epos[0] = mom_e[0]*Epos[2]/mom_e[2];
	  
	  Ppos[2] = _BeamCalz;
	  Ppos[1] = mom_p[1]*Ppos[2]/mom_p[2];
	  Ppos[0] = mom_p[0]*Ppos[2]/mom_p[2];
	  
	  // Making sure the transofmation is working
	  //cout << "Initial Position" << "(" << pos[0] << "," << pos[1] << "," <<  pos[2] << ")." << endl;
	  
	  //collect parameters for Lorentz transform
	  ein_x = mom_e[0];
	  pin_x = mom_p[0];
	  ein_energy = high_e->getEnergy();
	  pin_energy = high_p->getEnergy();
	  
	  //cout << "EEnergy in: " << ein_energy << ", x momentum: " << ein_x << "." << endl;
	  //cout << "PEnergy in: " << pin_energy << ", x momentum: " << pin_x << "." << endl;
	  
	  //apply the transform
	  scipp_ilc::transform_to_lab(ein_x, ein_energy, eout_x, eout_energy);
	  scipp_ilc::transform_to_lab(pin_x, pin_energy, pout_x, pout_energy);
	  
	  //adjust the x position
	  Epos[0] = eout_x*Epos[2]/mom_e[2];
	  Ppos[0] = pout_x*Ppos[2]/mom_p[2];
	  
	  //shift origin to the center of the beamcal beampipe hole
	  scipp_ilc::z_to_beam_out(Epos[0], Epos[1], Epos[2]);
	  scipp_ilc::z_to_beam_out(Ppos[0], Ppos[1], Ppos[2]);
	  
	  
	  //adjust the energy to post transform
	  Eenergy=eout_energy;
	  Penergy=pout_energy;

	  //find the new angle
	  double tmag_e = sqrt(pow(eout_x, 2)+pow(mom_e[1], 2));
	  double tmag_p = sqrt(pow(pout_x, 2)+pow(mom_p[1], 2));
	  emag+=tmag_e;
	  pmag+=tmag_p;
	  etheta = atan(tmag_e/abs(mom_e[2]));
	  ptheta = atan(tmag_p/abs(mom_p[2]));
	  
	  //Filling in the endpoints (Z) and hitmap (x,y) plots.
	  _hitmap->Fill(Epos[0],Epos[1]);
	  _Zpos->Fill(Epos[2]);
	  _hitmap->Fill(Ppos[0],Ppos[1]);
	  _Zpos->Fill(Ppos[2]);

	  //TENTATIVE. Getting the hit status of the e- and e+
	  
	  int Estatus = scipp_ilc::get_hitStatus(Epos[0],Epos[1]);
	  int Pstatus = scipp_ilc::get_hitStatus(Ppos[0],Ppos[1]);
	  
	  if(Estatus==1){Ehit_status=true;}else{Ehit_status=false;}
	  if(Pstatus==1){Phit_status=true;}else{Phit_status=false;}
	  
	  //Increasing the counters
	  if(Phit_status){
	    if(Ehit_status){_nHitHit++;
	    }else{_nPHitEMiss++;}
	  }else{
	    if(Ehit_status){_nEHitPMiss++;
	    }else{_nMissMiss++;}
	  }
	  _nTotal++;

	  
	  
        }
	
    }//end collection
    _nEvt ++ ;
}//end process



void BhaBhaDeflectionAnalysis::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void BhaBhaDeflectionAnalysis::end(){
    _rootfile->Write();
    if (_nTotal > 0){
      double Phh = 100.0 * _nHitHit / _nTotal;
      double Pphem = 100.0 * _nPHitEMiss / _nTotal;
      double Pehpm = 100.0 * _nEHitPMiss / _nTotal;
      double Pmm = 100.0 * _nMissMiss / _nTotal;

      cout << "Percent hit-hit: %" << Phh << endl;
      cout << "Percent P-hit E-miss: %" << Pphem << endl;
      cout << "Percent E-hit P-miss: %" << Pehpm << endl;
      cout << "Percent miss-miss: %" <<Pmm << endl;
    }
    cout << "That's all folks!" << endl;
}

