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
//static TH2F* _hitmap;
static TH1F* _mass;
//static TH1F* _scalar;
static TH1F* _vector;

static TH1F* _xSum;
static TH1F* _ySum;

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



void ParticleTree::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 

void ParticleTree::getChildren(MCParticle* hit, int gen){
    int i=1;
    for(MCParticle* kid : hit->getDaughters()){
        for(int i=0 ; i<=gen; i++){
            cout << "    ";
        }    

        cout << "Child " << i << "     id: " << kid->getPDG() << "     status: " << kid->getGeneratorStatus();
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

    
    double scatter_vec[] = {0, 0, 0};
    double mag = 0;
    double energy = 0;
    double theta;
    int id, stat;

    MCParticle* high_e;
    MCParticle* high_p;

    const double* mom;

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
        
        //create sum vector
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
    
           id = hit->getPDG();
           stat = hit->getGeneratorStatus();
           mom = hit->getMomentum();
           if(stat==1){
                    scatter_vec[0]+=mom[0];
                    scatter_vec[1]+=mom[1];
                    scatter_vec[2]+=mom[2];
                    energy+=hit->getEnergy();
           }//end final state
        }//end for loop
        
        cout << "event = " << _nEvt << endl;
        
        if(abs(scatter_vec[0])>1.9){
            _b_def_count++;
            cout << endl;
            cout << endl;
            cout << endl;

            //particle printer
            for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
               MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
        
               id = hit->getPDG();
               stat = hit->getGeneratorStatus();
               mom = hit->getMomentum(); 
               cout << "n: " << hitIndex << "  id: " << id << "  stat: " << stat;
               cout << "  momentum: [" << mom[0] << ", " << mom[1] << ", " << mom[2] << "]  energy: " << hit->getEnergy();   
               for(MCParticle* parent : hit->getParents()){
                    cout << "  parent [id, energy]: [" << parent->getPDG() << ", " << parent->getEnergy() << "]";  
               }
               cout << endl;

            }//end for loop
            cout << "XSum: " << scatter_vec[0] << endl;
            cout << "YSum: " << scatter_vec[1] << endl;
        }//end delta cut
        const double* mom_e = high_e->getMomentum();
        const double* mom_p = high_p->getMomentum();

        //all
        if(_nEvt<200000){
            
            _xSum->Fill(scatter_vec[0]);
            _ySum->Fill(scatter_vec[1]);

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



void ParticleTree::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void ParticleTree::end(){
    _rootfile->Write();
    cout << "Number of events with abs(x-mom sum) > 1.9: " << _b_def_count << endl;
    double rat = _b_def_count/_nEvt;
    cout << "Ratio of events with abs(x-mom sum) > 1.9: " << rat << endl;
}

