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
#include "include/Thrust.h"
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


int numEvents = 396;
int modEvents = 99;


double* pX = new double[modEvents];
double* pY = new double[modEvents];
double* nX = new double[modEvents];
double* nY = new double[modEvents];
double* peventBarycenterX = new double[modEvents];
double* peventBarycenterY = new double[modEvents];
double* neventBarycenterX = new double[modEvents];
double* neventBarycenterY = new double[modEvents];
double* pEnergyDep = new double[modEvents];
double* nEnergyDep = new double[modEvents];
double* pLR = new double[modEvents];
double* nLR = new double[modEvents];
double* pTD = new double[modEvents];
double* nTD = new double[modEvents];
double* pmeanDepth = new double[modEvents];
double* nmeanDepth = new double[modEvents];
double* prmoment = new double[modEvents];
double* nrmoment = new double[modEvents]; 
double* pinvrmoment = new double[modEvents];
double* ninvrmoment = new double[modEvents];
double* pThrustValue = new double[modEvents];
double* nThrustValue = new double[modEvents];

double* pScenBarycenterX = new double[numEvents / modEvents];
double* pScenBarycenterY = new double[numEvents / modEvents];
double* nScenBarycenterX = new double[numEvents / modEvents];
double* nScenBarycenterY = new double[modEvents / modEvents];
double* pScenEnergyDep = new double[numEvents / modEvents];
double* nScenEnergyDep = new double[numEvents / modEvents];
double* pScenLR = new double[numEvents / modEvents];
double* nScenLR = new double[numEvents / modEvents];
double* pScenTD = new double[numEvents / modEvents];
double* nScenTD = new double[numEvents / modEvents];
double* pScenMeanDepth = new double[numEvents / modEvents];
double* nScenMeanDepth = new double[numEvents / modEvents];
double* pScenrmoment = new double[numEvents / modEvents];
double* nScenrmoment = new double[numEvents / modEvents];
double* pSceninvrmoment = new double[numEvents / modEvents];
double* nSceninvrmoment = new double[numEvents / modEvents];




//double* firstBarycenters = new double[4];

int currentEvent = 0;
int currentScen = 0;



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

//pEnerygDep, nEnergyDep, pLR, nLR, pTD, nTD, pmeanDepth, nmeanDepth, prmoment, nrmoment, pinvrmoment, ninvrmoment, pX, pY, nX, nY   
double* ElliotsAnalysis::calculateObservables(LCCollection* col, double barycenters[4]){

  double ptotalEnergy = 0, ntotalEnergy = 0;
  double pLR = 0, nLR = 0, pTD = 0, nTD = 0;
  double pmeanDepth = 0, nmeanDepth = 0;
  double prmoment = 0, nrmoment = 0;
  double pinvrmoment = 0, ninvrmoment = 0;
  double pX = 0, pY = 0, nX = 0, nY = 0;

  //LR
  double pnum_LR = 0, pdenom_LR = 0, nnum_LR = 0, ndenom_LR = 0;
  //TD
  double pnum_TD = 0, pdenom_TD = 0, nnum_TD = 0, ndenom_TD = 0; 
  //meanDepth
  double pnum_meanDepth = 0, nnum_meanDepth = 0;
  //r-moment
  double pnum_moment = 0, pdenom_moment = 0, nnum_moment = 0, ndenom_moment = 0;
  //invr-moment
  double pnum_invrmoment = 0, nnum_invrmoment = 0;
  
  double prad = 0, nrad = 0;

  //Thrust
  double pThrustValue = 0, pThrustAxis = 0;
  double nThrustValue = 0, nThrustAxis = 0;
 
  if( col != NULL ){
    int nElements = col->getNumberOfElements()  ;

    for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
      SimCalorimeterHit* hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt(hitIndex) );
      double currentEnergy = hit->getEnergy();
      double currentPosX = hit->getPosition()[0];
      double currentPosY = hit->getPosition()[1];
      double currentPosZ = hit->getPosition()[2];
      currentPosX = currentPosX - std::abs(hit->getPosition()[2] * 0.007);

     
      if (hit->getPosition()[2] > 0){
	currentPosX = currentPosX - barycenters[0];				
	currentPosY = currentPosY - barycenters[1];
 
	//Positions
	pX = currentPosX;
	pY = currentPosY;

	//total  energy deposit
	ptotalEnergy += currentEnergy;
	
	//LR
	pnum_LR += currentPosX * currentEnergy;
	pdenom_LR += std::abs(currentPosX) * currentEnergy;

	//TD
	pnum_TD += currentPosY * currentEnergy;
	pdenom_TD += std::abs(currentPosY) * currentEnergy;

	//Mean Depth
	pnum_meanDepth += currentPosZ * currentEnergy;

	//r-moment
	prad = std::sqrt((std::pow(currentPosX,2) + (std::pow(currentPosY,2))));
	pnum_moment += prad * currentEnergy;
	pdenom_moment += currentEnergy;

	//invr-moment
	pnum_invrmoment += (1 / prad) * currentEnergy;

	//Thrust
	pThrustValue = 5;//processThrust(currentEnergy*currentPosX, currentEnergy*currentPosY, currentPosZ);

      }else {

	currentPosX = currentPosX - barycenters[2];
        currentPosY = currentPosY - barycenters[3];


	//Positions                                                                                                                        
	nX = currentPosX;
	nY = currentPosY;

	//total energy deposit
	ntotalEnergy += currentEnergy;

	//LR                                                                                                                               
        nnum_LR+= currentPosX * currentEnergy;
	ndenom_LR += std::abs(currentPosX) * currentEnergy;

	//TD                                                                                                                               
	nnum_TD+= currentPosY * currentEnergy;
	ndenom_TD += std::abs(currentPosY) * currentEnergy;

	//Mean Depth                                                                                                                       
	nnum_meanDepth += currentPosZ * currentEnergy;

	//r-moment
	nrad = std::sqrt((std::pow(currentPosX,2) + (std::pow(currentPosY,2))));
	nnum_moment += nrad * currentEnergy;
        ndenom_moment += currentEnergy;

	//invr-moment                                                                                                                      
	nnum_invrmoment += (1 / nrad) * currentEnergy;
      
	//Thrust                                                                                                                            
        nThrustValue = 5;//Thrust::processThrust(currentEnergy*currentPosX, currentEnergy*currentPosY, currentPosZ);
      }


    }
  }

  double* obs = new double[18];

  //Energy Deposit
  obs[0] = ptotalEnergy;
  obs[1] = ntotalEnergy;
  
  //LR
  pLR = pnum_LR / pdenom_LR;
  nLR = nnum_LR / ndenom_LR;
  obs[2] = pLR;
  obs[3] = nLR;

  //TD
  pTD = pnum_TD / pdenom_TD;
  nTD = nnum_TD / ndenom_TD;
  obs[4] = pTD;
  obs[5] = nTD;

  //mean Depth
  pmeanDepth = pnum_meanDepth / ptotalEnergy;
  nmeanDepth = nnum_meanDepth / ntotalEnergy;
  obs[6] = pmeanDepth;
  obs[7] = nmeanDepth;

  //r-moment
  prmoment = pnum_moment / pdenom_moment;
  nrmoment = nnum_moment / ndenom_moment;
  obs[8] = prmoment;
  obs[9] = nrmoment;

  //invr-moment
  pinvrmoment = pnum_invrmoment / ptotalEnergy;
  ninvrmoment = nnum_invrmoment / ntotalEnergy;
  obs[10] = pinvrmoment;
  obs[11] = ninvrmoment;

  //Positons
  obs[12] = pX;
  obs[13] = pY;    
  obs[14] = nX;
  obs[15] = nY;

  //Thrust
  obs[16] = pThrustValue;
  obs[17] = nThrustValue;

  return obs;
  
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

    double* barycenters = new double[4]; 
    barycenters[0] = pbarycenterPosX;
    barycenters[1] = pbarycenterPosY;
    barycenters[2] = nbarycenterPosX;
    barycenters[3] = nbarycenterPosY;
    
    //    printf("\n\nPositive Barycenter Position: (%f, %f) with Energy: %f\n\n", barycenters[0],barycenters[1], pEnergy);
    // printf("\n\nNegative Barycenter Position: (%f, %f) with Energy: %f\n\n", barycenters[2],barycenters[3], nEnergy);
    
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

double ElliotsAnalysis::findAvgObs(double* obs){
  double sum = 0;
  double avgObs = 0;

  int num = 0;

  for (int i = 0;i < modEvents ; i++ ){

    sum += obs[i];

    num++;

  }

  avgObs = sum / num;

  return avgObs;

}



void ElliotsAnalysis::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    LCCollection* col = evt->getCollection( _colName ) ;
    

    

        
    double* barycenters = calculateBarycenter(col);    
    double* obs;



   
    double firstBarycenters[] = { -15.036040, -8.101348, -15.051319,  8.064840};

    obs = calculateObservables(col, firstBarycenters);


   
    peventBarycenterX[currentEvent] = barycenters[0];
    peventBarycenterY[currentEvent] = barycenters[1];
    neventBarycenterX[currentEvent] = barycenters[2];
    neventBarycenterY[currentEvent] = barycenters[3];
    pEnergyDep[currentEvent] = obs[0];
    nEnergyDep[currentEvent] = obs[1];
    pLR[currentEvent] = obs[2];
    nLR[currentEvent] = obs[3];
    pTD[currentEvent] = obs[4];
    nTD[currentEvent] = obs[5];
    pmeanDepth[currentEvent] = obs[6];
    nmeanDepth[currentEvent] = obs[7];
    prmoment[currentEvent] = obs[8];
    nrmoment[currentEvent] = obs[9];
    pinvrmoment[currentEvent] = obs[10];
    ninvrmoment[currentEvent] = obs[11];
    pX[currentEvent] = obs[12];
    pY[currentEvent] = obs[13];
    nX[currentEvent] = obs[14];
    nY[currentEvent] = obs[15];
    pThrustValue[currentEvent] = obs[16];
    nThrustValue[currentEvent] = obs[17];
    

    /*
    printf("\nBARYCENTER Postive:( %f,%f) Negative (%f,%f)", barycenters[0], barycenters[1], barycenters[2], barycenters[3]);
    printf("ENERYG DEPOSIT Postive: %f  Negative: %f", obs[0], obs[1]);
    */
    //    printf("THRUST VALUE Postive: %f  Negative: %f", obs[16], obs[17]);
    
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
   
    // printParticleProperties(maxHit);
  
    currentEvent++;
    
    if (currentEvent % modEvents == 0){


     
      double pAvgBarycenterX = findAvgObs(peventBarycenterX);
      double pAvgBarycenterY = findAvgObs(peventBarycenterY);
      double nAvgBarycenterX = findAvgObs(neventBarycenterX);
      double nAvgBarycenterY = findAvgObs(neventBarycenterY);
      double pAvgEnergyDep = findAvgObs(pEnergyDep);
      double nAvgEnergyDep = findAvgObs(nEnergyDep);
      double pAvgLR = findAvgObs(pLR);
      double nAvgLR = findAvgObs(nLR);
      double pAvgTD = findAvgObs(pTD);
      double nAvgTD = findAvgObs(nTD);
      double pAvgmeanDepth = findAvgObs(pmeanDepth);
      double nAvgmeanDepth = findAvgObs(nmeanDepth);
      double pAvgrmoment = findAvgObs(prmoment);
      double nAvgrmoment = findAvgObs(nrmoment);
      double pAvginvrmoment = findAvgObs(pinvrmoment);
      double nAvginvrmoment = findAvgObs(ninvrmoment);
      double pAvgX = findAvgObs(pX);
      double pAvgY = findAvgObs(pY);
      double nAvgX = findAvgObs(nX);
      double nAvgY = findAvgObs(nY);

      
      pScenBarycenterX[currentScen] = pAvgBarycenterX;
      pScenBarycenterY[currentScen] = pAvgBarycenterY;
      nScenBarycenterX[currentScen] = nAvgBarycenterX;
      nScenBarycenterY[currentScen] = nAvgBarycenterY;
      pScenEnergyDep[currentScen] = pAvgEnergyDep;
      nScenEnergyDep[currentScen] = nAvgEnergyDep;
      pScenLR[currentScen] = pAvgLR;
      nScenLR[currentScen] = nAvgLR;
      pScenTD[currentScen] = pAvgTD;
      nScenTD[currentScen] = nAvgTD;
      pScenMeanDepth[currentScen] = pAvgmeanDepth;
      nScenMeanDepth[currentScen] = nAvgmeanDepth;
      pScenrmoment[currentScen] = pAvgrmoment;
      nScenrmoment[currentScen] = nAvgrmoment;
      pSceninvrmoment[currentScen] = pAvginvrmoment;
      nSceninvrmoment[currentScen] = nAvginvrmoment;

      printf("\n=====================Scenario %d========================== \n", currentScen + 1);

      printf("\nAVERAGE BARYCENTER: Postive: (%f,%f) Negative: (%f,%f)", pAvgBarycenterX, pAvgBarycenterY,nAvgBarycenterX, nAvgBarycenterY);
      printf("\nAVERAGE Energy Deposit: Postive: %f Negative: %f", pAvgEnergyDep, nAvgEnergyDep);

      printf("\nAVERAGE R-MOMENT: Postive: %f Negative: %f", pAvgrmoment, nAvgrmoment);

      printf("\nAVERAGE Mean Depth: Postive: %f Negative: %f", pAvgmeanDepth, nAvgmeanDepth);

      printf("\nAVERAGE LR: Postive: %f Negative: %f", pAvgLR, nAvgLR);

      printf("\nAVERAGE TD: Postive: %f Negative: %f\n", pAvgTD, nAvgTD);


      currentScen++;
      std::fill_n(peventBarycenterX, modEvents, 0);
      std::fill_n(peventBarycenterY, modEvents, 0);
      std::fill_n(neventBarycenterX, modEvents, 0);
      std::fill_n(neventBarycenterY, modEvents, 0);
      std::fill_n(pEnergyDep, modEvents, 0);
      std::fill_n(nEnergyDep, modEvents, 0);
      std::fill_n(pLR, modEvents, 0);
      std::fill_n(nLR, modEvents, 0);
      std::fill_n(pTD, modEvents, 0);
      std::fill_n(nTD, modEvents, 0);
      std::fill_n(pmeanDepth, modEvents, 0);
      std::fill_n(nmeanDepth, modEvents, 0);
      std::fill_n(prmoment, modEvents, 0);
      std::fill_n(nrmoment, modEvents, 0);
      std::fill_n(pinvrmoment, modEvents, 0);
      std::fill_n(ninvrmoment, modEvents, 0);
      std::fill_n(pX, modEvents, 0);
      std::fill_n(pY, modEvents, 0);
      std::fill_n(nX, modEvents, 0);
      std::fill_n(nY, modEvents, 0);
  
        currentEvent = 0;

      /*
      if (currentScen = (numEvents / modEvents)){
	currentScen = 0;
	}*/

	if (currentScen == (numEvents / modEvents)){
	  printf("\nFirst BARYCENTER: Postive: (%f,%f) Negative: (%f,%f)\n", pScenBarycenterX[0],  pScenBarycenterY[0],nScenBarycenterX[0],
		 nScenBarycenterY[0]);


	}

    }
    // if (currentScen = (numEvents / modEvents)){
    // currentScen = 0;
    
      //   printf("\nScenario: %d R-Moment: Postive: %f Negative: %f\n", currentScen, pScenrmoment, nScenrmoment);
      // }


   


    

    _nEvt++ ;
}



void ElliotsAnalysis::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void ElliotsAnalysis::end(){


 
    _rootfile->Write();
}
