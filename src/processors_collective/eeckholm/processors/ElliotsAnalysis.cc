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

#include "ElliotsAnalysis.h"
#include "scipp_ilc_utilities.h"
#include <iostream>

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


ElliotsAnalysis ElliotsAnalysis;

static TFile* _rootfile;
static TH2F* _hitmap;


ElliotsAnalysis::ElliotsAnalysis() : Processor("ElliotsAnalysis") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}



void ElliotsAnalysis::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("Elliothitmap.root","RECREATE");
    _hitmap = new TH2F("Elliothitmap","Hit Distribution",300.0,-150.0,150.0,300.0,-150.0,150.0);

    // usually a good idea to
    //printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;

}



void ElliotsAnalysis::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 

void ElliotsAnalysis::printParticleProperties(SimCalorimeterHit* hit){

  int highestParticleEnergy = 0;
  int type = 0;
  double energy = 0;
  float charge = 0;
  float px = 0, py = 0, pz = 0;

  MCParticle* particle; 

    for (int i = 0; i < hit->getNMCContributions(); i++){

      particle = hit->getParticleCont(i);
      type = particle->getPDG();
     
      px = particle ->getMomentum()[0];
      py = particle->getMomentum()[1];
      pz = particle->getMomentum()[2];
      charge = particle->getCharge();

      if (particle->getEnergy() > highestParticleEnergy){
	highestParticleEnergy = particle->getEnergy();
      }
    }

    printf("");

}

void ElliotsAnalysis::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    LCCollection* col = evt->getCollection( _colName ) ;
    double highestEnergy = 0;
    double lowestEnergy = 10000;
    double hParticleEnergy = 0;
    double lParticleEnergy = 0;

    float hPosX = 0, hPosY = 0, hPosZ = 0;
    float lPosX = 0, lPosY= 0, lPosZ = 0;

    SimCalorimeterHit* maxHit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt(0));
    SimCalorimeterHit* minHit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt(0));

   
    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
	
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           SimCalorimeterHit* hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt(hitIndex) );
	   
	   
	   //find hit with highest energy
	   if (hit->getEnergy() > highestEnergy){
	     highestEnergy = hit->getEnergy();
	     maxHit = hit;
	     hPosX = hit->getPosition()[0];
	     hPosY = hit->getPosition()[1];
	     hPosZ = hit->getPosition()[2];
	   }

	   //find hit with lowest energy
	   if (hit->getEnergy() < lowestEnergy){
             lowestEnergy = hit->getEnergy();
	     minHit = hit;
             lPosX = hit->getPosition()[0];
             lPosY = hit->getPosition()[1];
             lPosZ = hit->getPosition()[2];
           }

	   

           const float* pos = hit->getPosition();
           _hitmap->Fill(pos[0],pos[1]);
        } 
    }
    for (int i = 0; i < maxHit->getNMCContributions(); i++){
      hParticleEnergy += maxHit->getEnergyCont(i);

    }

    // printMaxParticle(SimCalorimetHit* hit)
    for (int i = 0; i < minHit->getNMCContributions(); i++){
      //  lParticleEnergy += minHit->getEnergyCont(i);
      MCParticle*  lParticle = minHit->getParticleCont(i);
      printf("Particle Energy: %0.25f\n", lParticle->getEnergy());
    }
    
    double highestParticleEnergy =0;
    for (int i = 0; i < maxHit->getNMCContributions(); i++){
                                                                                                              
      MCParticle*  hParticle = maxHit->getParticleCont(i);
      if (hParticle->getEnergy() > highestParticleEnergy){
	highestParticleEnergy =hParticle->getEnergy();
	
      }
      
    }

    printf("Highest Energy: %0.15f\n", highestEnergy);
    printf("Highest Energy Position:( %f,%f,%f) \n", hPosX, hPosY, hPosZ);
    printf("Sum of %d Particle Energies: %f\n",maxHit->getNMCContributions(), hParticleEnergy);
    printf("Highest Particle Energy: %0.25f\n", highestParticleEnergy);

    printf("Lowest Energy: %0.25f\n", lowestEnergy);
    printf("Lowest Energy Position:( %f,%f,%f) \n", lPosX, lPosY, lPosZ);

    
    _nEvt ++ ;
}



void ElliotsAnalysis::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void ElliotsAnalysis::end(){ 
    _rootfile->Write();
}
