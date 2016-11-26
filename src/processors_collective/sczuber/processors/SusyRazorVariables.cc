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
 * author Summer Zuber
 * November 23, 2016
 * This is an initial attempt to use calculate and use Razor Variables with our degenerate Susy Events 
 */

#include "SusyRazorVariables.h"
#include "scipp_ilc_utilities.h"
#include <iostream>
#include <cmath>

#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/MCParticle.h>

#include <TFile.h>
#include <TH2D.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include <list>


using namespace lcio;
using namespace marlin;
using namespace std;

SusyRazorVariables SusyRazorVariables;

static TFile* _rootfile;

static TH1F* _FirstRazorPlot;

SusyRazorVariables::SusyRazorVariables() : Processor("SusyRazorVariables") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
}



void SusyRazorVariables::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("SusyRazorVariables.root","RECREATE");
    _FirstRazorPlot = new TH1F("FirstRazorPlot","My First Razor Plot", 40,0,20); 
    // usually a good idea to
    //printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;
}



void SusyRazorVariables::processRunHeader( LCRunHeader* run) { 
    //    _nRun++ ;
} 



void SusyRazorVariables::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    LCCollection* col = evt->getCollection( _colName ) ;
    cout << endl;
    cout << endl;
    cout << endl;
    cout << "event = " << _nEvt << endl;

    //creat tau list? 

    // might not need these:
    double vec[4][3];  // momentum 3 vectors of: tau 1, tau 2, lsp 1, lsp 2
    double scalars[4]; // magnitude of 3 vectors of same categories 
    double energy[4];  // energy of same categories 

    //particle identifiers
    int id, stat; 

    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;

        // For each particle in Event ...
        for(int particleIndex = 0; particleIndex < nElements ; particleIndex++){
            MCParticle* particle = dynamic_cast<MCParticle*>( col->getElementAt(particleIndex) ); 
            try{ 
                id = particle->getPDG(); 
                stat = particle->getGeneratorStatus();
            }
            catch(const std::exception& e){
                cout << "exception caught with message " << e.what() << "\n";
            }
           

            if (id==15 || id == -15){
                cout << particle << " " << particle->getPDG() << endl;
                for(MCParticle* parent : particle->getParents()){
                    cout << "parent: parent, id" << parent << " " << parent->getPDG() << endl;
                    if(parent->getPDG() == 1000015 || parent->getPDG() == -1000015){
                        //tauLIST.ADD(PARTICLE)
                    }
                }

            }

            // transform momentum-energy four vector to R frame:
            // beta = (E of tau 1 - E of tau 2)/(pz of tau 1 - pz of tau 2)
            if(id == 1000022){
                cout << particle <<"  "<< particle->getPDG() << endl; 
            }
            // _V_n_C->Fill(total_detected_vector);


        }//end for
    }//ind if col
    _nEvt ++;
    cout << "event "<< _nEvt <<" finished " << endl;
}//end process
void SusyRazorVariables::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void SusyRazorVariables::end(){ 

    _rootfile->Write();
}
