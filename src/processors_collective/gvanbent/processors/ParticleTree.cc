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

#include "ParticleTree.h"
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

ParticleTree ParticleTree;

static TFile* _rootfile;

ParticleTree::ParticleTree() : Processor("ParticleTree") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
}

void ParticleTree::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("eBpW_dump.root","RECREATE");
    //_hitmap = new TH2F("hitmap","Hit Distribution",600.0,-300.0,300.0,600.0,-300.0,300.0);
    //_scalar = new TH1F("scalar", "Transverse Momentum Scalar Magnitude", 2000.0, 0.0, 20.0);
    
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



void ParticleTree::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 

void ParticleTree::getChildren(MCParticle* hit, int gen){
    int i=1;
    for(MCParticle* kid : hit->getDaughters()){
        for(int i=0 ; i<=gen; i++){
            cout << "    ";
        }    

        cout << "Child " << i << "     id: " << kid->getPDG() << "     stat: " << kid->getGeneratorStatus();
        cout << "     energy:" << kid->getEnergy() << "     children: " << kid->getDaughters().size() << endl;
        i++;
        int this_gen = gen;
        if(kid->getGeneratorStatus()==2){gen++; getChildren(kid, gen);}
        gen=this_gen;
    }
    cout << endl;
}

void ParticleTree::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...


    LCCollection* col = evt->getCollection( _colName ) ;

    int id, stat;
    const double* mom;


    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
        
        cout << "event = " << _nEvt << endl;
        

        //particle printer
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );

           id = hit->getPDG();
           stat = hit->getGeneratorStatus();
           mom = hit->getMomentum(); 
           cout << "n: " << hitIndex << "  id: " << id << "  stat: " << stat;
           cout << "  mom: [" << mom[0] << ", " << mom[1] << ", " << mom[2] << "]  en: " << hit->getEnergy();   
           cout << " children: " << hit->getDaughters().size() << endl;
           /*for(MCParticle* parent : hit->getParents()){
                cout << "  parent [id, energy]: [" << parent->getPDG() << ", " << parent->getEnergy() << "]";  
           }*/
           if(hit->getParents().size()==0){
                cout << endl;
                getChildren(hit, 1);                   
           }
           cout << endl;

        }//end for loop

    }//end collection
    _nEvt ++ ;
}//end process



void ParticleTree::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void ParticleTree::end(){
    _rootfile->Write();
}

