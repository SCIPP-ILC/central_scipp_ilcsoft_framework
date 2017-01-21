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

#include "Prediction.h"
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


Prediction Prediction;

static TFile* _rootfile;
static TH2F* _hitmap;
static TH1F* _vector;

Prediction::Prediction() : Processor("Prediction") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}



void Prediction::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("eBpB_vector.root","RECREATE");
    // usually a good idea to
    //printParameters() ;
    _nEvt = 0 ;

}



void Prediction::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void Prediction::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    LCCollection* col = evt->getCollection( _colName ) ;

   /* vector<int*> particles;
    int stat, id =0;
    double tot_mom[]={0, 0};
    double compEn_e=0;
    double compEn_p=0;
    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;

        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
            MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
            MyParticle particle = new MyParticle(hit);

            id = particle.getPDG();
            stat = particle.getStat();
            if(stat==1){
                //only final state particles in vector
                particles.push_back(&particle); 
                cout << "Particle " << hitIndex << " with ID: " << id << endl;
                //find highest energy electron and positron
                if(id==11){
                    if(particle.getEnergy() > compEn_e){compEn_e=particle.getEnergy();}    
                }
                else if(id==-11){
                    if(particle.getEnergy() > compEn_p){compEn_p=particle.getEnergy();}    
                }
                //set all other particle types to hadronic system
                else{
                    particle.setHadronic(true);
                    particle.setDetectable(true);
                }
            }//end final state   
        }//end for

        for(MyParticle* particle : particles){
            if(id==11 && particle->getEnergy()==compEn_e){particle->setElectronic(true);}    
            else if(id==-11 && particle->getEnergy()==compEn_p){particle->setElectronic(true);}
            else{
                if(abs(id)==12||abs(id)==14||abs(id)==16||abs(id)==18){particle->setDetectable(false);}    
            }
        }

        for(MyParticle* particle : particles){
            if(particle->isElectronic){
                cout << "ELECTRONIC:: id: " << particle->id() << " Energy: " << particle->energy() << endl;
            }    
            if(particle->isHadronic){
                cout << "HADRONIC:: id: " << particle->id() << " Energy: " << particle->energy() << endl;
            }
        }
         
    }*/

    _nEvt ++ ;
}



void Prediction::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void Prediction::end(){ 
    _rootfile->Write();
}
