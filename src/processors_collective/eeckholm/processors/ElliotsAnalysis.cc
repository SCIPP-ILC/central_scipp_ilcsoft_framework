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
#include <math.h>

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

double* ElliotsAnalysis::calculateMoments(LCCollection* col, double barycenters[4]){
  double moments[2];
  double pmom = 0, nmom = 0;
  double pnum = 0, pdenom = 0, nnum = 0, ndenom = 0;
  
  if( col != NULL ){
    int nElements = col->getNumberOfElements()  ;

    for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
      SimCalorimeterHit* hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt(hitIndex) );
     
      double prad = std::sqrt((std::pow(barycenters[0],2)) + (std::pow(barycenters[1],2)));
      double nrad = std::sqrt((std::pow(barycenters[2],2)) + (std::pow(barycenters[4],2)));
      
	pnum += prad * hit->getEnergy();
        pdenom += hit->getEnergy();

        nnum += nrad * hit->getEnergy();
	ndenom += hit->getEnergy();
    }
  }
  
  pmom = pnum / pdenom;
  nmom = nnum / ndenom;
  
  moments[0] = pmom;
  moments[1] = nmom;

  return moments;
  
}

double*  ElliotsAnalysis::calculateBarycenter( LCCollection* col ){
  
  double pbarycenterPosX = 0, pbarycenterPosY = 0, nbarycenterPosX = 0, nbarycenterPosY = 0;
  double pnumX = 0, pnumY = 0, pdenomX = 0, pdenomY = 0, nnumX = 0, nnumY = 0, ndenomX = 0, ndenomY = 0;
  double pEnergy = 0, nEnergy = 0;
    
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;

        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
            SimCalorimeterHit* hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt(hitIndex) );
	    double currentEnergy = hit->getEnergy();
	    double currentPosX = hit->getPosition()[0];
	    double currentPosY = hit->getPosition()[1];
	    double currentPosZ = hit->getPosition()[2];

	    currentPosX = currentPosX - std::abs(currentPosZ * 0.007);

	    if (currentPosZ < 0){
	      //calculate numerator and denominator of barycenter x value                                                                              
	      nnumX += currentPosX * currentEnergy;
	      ndenomX += currentEnergy;

	      //calculate numerator and denominator of barycenter y value                                                                              
	      nnumY += currentPosY * currentEnergy;
	      ndenomY += currentEnergy;
	    }
	    else {

	      //calculate numerator and denominator of barycenter x value
	      pnumX += currentPosX * currentEnergy;
	      pdenomX += currentEnergy;

	      //calculate numerator and denominator of barycenter y value
	      pnumY += currentPosY * currentEnergy;
	      pdenomY += currentEnergy;
	    }
	 }
    }
    
    pEnergy = pdenomX;
    pbarycenterPosX = pnumX / pdenomX;
    pbarycenterPosY = pnumY / pdenomY;

    nEnergy = ndenomX;
    nbarycenterPosX = nnumX / ndenomX;
    nbarycenterPosY = nnumY / ndenomY;

    double* barycenters; 
    barycenters[0] = pbarycenterPosX;
    barycenters[1] = pbarycenterPosY;
    barycenters[2] = nbarycenterPosX;
    barycenters[3] = nbarycenterPosY;
    
    printf("\n\nPositive Barycenter Position: (%f, %f) with Energy: %f\n\n", barycenters[0],barycenters[1], pEnergy);
    printf("\n\nNegative Barycenter Position: (%f, %f) with Energy: %f\n\n", nbarycenterPosX, nbarycenterPosY, nEnergy);
    
    return barycenters;
} 

void ElliotsAnalysis::printParticleProperties(SimCalorimeterHit* hit){

  
    int type = 0;
    double energy = 0;
    float charge = 0;
    float px = 0, py = 0, pz = 0;

    MCParticle* currentParticle; 
    MCParticle* highestEnergyParticle = hit->getParticleCont(0);
    
    

    for (int i = 0; i < hit->getNMCContributions(); i++){
      currentParticle = hit->getParticleCont(i);
      

      if (currentParticle->getEnergy() > highestEnergyParticle->getEnergy()){
        highestEnergyParticle = currentParticle;
      }

    }
    energy = highestEnergyParticle->getEnergy();
    type =  highestEnergyParticle->getPDG(); 
    px =  highestEnergyParticle->getMomentum()[0];
    py =  highestEnergyParticle->getMomentum()[1];
    pz =  highestEnergyParticle->getMomentum()[2];
    charge =  highestEnergyParticle->getCharge();
    

    printf("\nHighest energy Particle in hit: %0.5f\n", energy);
    printf("Type: %d\n",type);
    printf("Momentum: (%0.2f,%0.2f, %0.2f)\n", px,py,pz);
    printf("Charge: %0.2f\n", charge);
  
  
 

}

void ElliotsAnalysis::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    LCCollection* col = evt->getCollection( _colName ) ;
    
    double* barycenters;
    barycenters[0] = calculateBarycenter(col)[0];
    barycenters[1] = calculateBarycenter(col)[1];
    barycenters[2] = calculateBarycenter(col)[2];
    barycenters[3] = calculateBarycenter(col)[3];

    printf("Postive: %f %f Negative %f %f", barycenters[0], barycenters[1], barycenters[2], barycenters[3]);

    //    double *moments = calculateMoments(col, barycenters);

    // printf("\nPostive moment: %f \n", moments[0]);
    // printf("\nNegative moment: %f \n", moments[1]);



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
   
    printParticleProperties(maxHit);
    
    _nEvt ++ ;
}



void ElliotsAnalysis::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void ElliotsAnalysis::end(){ 
    _rootfile->Write();
}
