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

#include "Observables.h"
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
using namespace scipp_ilc;

Observables Observables;

static TFile* _rootfile;
static TH1F* _V;
static TH1F* _S;
static TH1F* _Vzoom;
static TH1F* _Szoom;


Observables::Observables() : Processor("Observables") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}



void Observables::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile=new TFile("WW_had.root", "RECREATE");
    _V = new TH1F("V", "Magnitude of Transverse Momentum Vector Sum, Hadronic System", 100, 0.0, 20.0);
    _S = new TH1F("S", "Sum of Transverse Momentum Magnitudes, Hadronic System", 100, 0.0, 20.0);

    _Vzoom = new TH1F("V", "Magnitude of Transverse Momentum Vector Sum, Hadronic System", 100, 0.0, 2.0);
    _Szoom = new TH1F("S", "Sum of Transverse Momentum Magnitudes, Hadronic System", 100, 0.0, 2.0);

    
    // usually a good idea to
    //printParameters() ;
    _nEvt = 0 ;

}



void Observables::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void Observables::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    LCCollection* col = evt->getCollection( _colName ) ;

    int stat, id = 0;
    double E_e, E_p = 0;
    //double tot[]={0, 0, 0};
    double tot_mom[]={0, 0};
    double tmag = 0;
    double S = 0;
    double V = 0;
    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
    cout << endl;
    cout << endl;
    cout << "***************************EVENT: " << _nEvt << "****************************" << endl;

    vector<MCParticle*> system;

        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
        
            id = hit->getPDG();
            stat = hit->getGeneratorStatus();
            if(stat==1){ 
                system.push_back(hit);
                if(id==11){E_e = (E_e < hit->getEnergy()) ? hit->getEnergy() : E_e;}
                if(id==-11){E_p = (E_p < hit->getEnergy()) ? hit->getEnergy() : E_p;}
                cout << "Particle " << hitIndex << " with ID: " << id;
                cout << " status: " << stat;

                /*double mom[4];
                mom[0] = hit->getMomentum()[0]; 
                mom[1] = hit->getMomentum()[1]; 
                mom[2] = hit->getMomentum()[2];
                mom[3] = hit->getEnergy();

                cout << " momentum [" << mom[0] << ", " << mom[1] << ", " << mom[2] << "] with energy: " << mom[3] << endl;
                
                double mag = sqrt(pow(mom[0], 2)+pow(mom[1], 2));
                tmag+=mag;
                tot_mom[0]+=mom[0];
                tot_mom[1]+=mom[1];*/
                
            }//end final state
        }//end for
        
        for(MCParticle* hit : system){
            if(hit->getEnergy()!=E_e && hit->getEnergy()!=E_p){
                    
                double mom[4];
                mom[0] = hit->getMomentum()[0]; 
                mom[1] = hit->getMomentum()[1]; 
                mom[2] = hit->getMomentum()[2];
                mom[3] = hit->getEnergy();

                cout << " momentum [" << mom[0] << ", " << mom[1] << ", " << mom[2] << "] with energy: " << mom[3] << endl;
                double mag = sqrt(pow(mom[0], 2)+pow(mom[1], 2));
                tmag+=mag;
                tot_mom[0]+=mom[0];
                tot_mom[1]+=mom[1];
            }    
        }

        V = sqrt(pow(tot_mom[0], 2)+pow(tot_mom[1], 2));
        _V->Fill(V); 
        _Vzoom->Fill(V); 
        cout << "V: " << V << endl;
        S = tmag;
        _S->Fill(S); 
        _Szoom->Fill(S); 
        cout << "S: " << S << endl;
    }

    _nEvt ++ ;
}



void Observables::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void Observables::end(){ 
    _rootfile->Write();
}
