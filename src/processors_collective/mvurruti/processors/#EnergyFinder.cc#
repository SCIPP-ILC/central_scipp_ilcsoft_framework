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

#include "EnergyFinder.h"
#include "scipp_ilc_utilities.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/MCParticle.h>

#include <TFile.h>
#include <TH2D.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#include "polar_coords.h"

using namespace lcio;
using namespace marlin;
using namespace std;


EnergyFinder EnergyFinder;

static TFile* _rootfile;
static TH2F* _hitmap;
static TH1D* _energy_e;
static TH1D* _energy_p;
static TH2D* _energy_rad;

EnergyFinder::EnergyFinder() : Processor("EnergyFinder") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}

void EnergyFinder::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile(_root_file_name.c_str(),"RECREATE");
    _hitmap = new TH2F("hitmap","Hit Distribution",300.0,-150.0,150.0,300.0,-150.0,150.0);
    _energy_e = new TH1D("_energy_e","Electron Energies",125,0.0,260.0);
    _energy_p = new TH1D("_energy_p","Positron Energies",125,0.0,260.0);
    _energy_rad = new TH2D("energy_rad","Energy versus Radius",125,0.0,2.6,125,0.0,260.0);

    // usually a good idea to
    //printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;

}







void EnergyFinder::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void EnergyFinder::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    LCCollection* col = evt->getCollection( _colName ) ;

    // this will only be entered if the collection is available
    double e_MaxEnerg = 0.0;
    double p_MaxEnerg = 0.0;

    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;

        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
	   // SimCalorimeterHit* hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt(hitIndex) );
	   // The line above casts the hits in the event as Sim Calorimeter Hits. 
	   // But your xml defines this processor as working with MCParticles. 
	   // So you have to cast these to MCParticles as below:
	   MCParticle* particle = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );	   

	   //These next lines pose another problem. MCParticle are 4-vectors (momentum, energy) of all the particles created after the beams interact.
	   //MCPartile do NOT have a position and so there is no getPosition() function defined for them. 
           //const float* pos = hit->getPosition();
           //_hitmap->Fill(pos[0],pos[1]);
	   //Instead, you can use the getEnergy() or getMomentum() functions. 
           const double* mom = particle->getMomentum();
	    double energ = particle->getEnergy();
	    double ID = particle->getPDG();
	//    double rad = sqrt(pow(mom[0],2) + pow(mom[1],2));
	    double radius,phi;
	    scipp_ilc::cartesian_to_polar(mom[0],mom[1],radius,phi);


	
	if( ID == 11) {
	   _energy_e->Fill(energ);
	     if(energ>e_MaxEnerg){
		e_MaxEnerg = energ;
}
}
	if( ID == -11) {
	   _energy_p->Fill(energ);
	     if(energ>p_MaxEnerg) {
		p_MaxEnerg = energ;
}
}
	
	_energy_rad->Fill(radius,energ);


	   //And let's go ahead and print the momentum just to make sure this is working. Notice that once I define a momentum array (as I did above),
	   //I can refer to each component of this array using index notation (i.e. mom[0] is the x-momentum)
	

	 if(hitIndex%10) {	
		
		cout << "Energy: " << energ << endl;  

		cout << "Momentum: [" << mom[0] << ", " << mom[1] << ", " << mom[2] << "]" << endl; 	
       
	}





	 } 
    }
cout << "Finished Event " << _nEvt << endl;
    _nEvt ++ ;
}



void EnergyFinder::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void EnergyFinder::end(){ 
    _rootfile->Write();
}
