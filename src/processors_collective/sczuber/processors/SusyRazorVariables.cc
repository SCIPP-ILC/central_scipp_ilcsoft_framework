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
    std::list<MCParticle*> tauList;
    //instead of list, these will be the two tau's for the event:
    MCParticle* tau1;
    MCParticle* tau2;
    // might not need these:
    double vec[4][3];  // momentum 3 vectors of: tau 1, tau 2, lsp 1, lsp 2
    double scalars[4]; // magnitude of 3 vectors of same categories 
    double energy[4];  // energy of same categories 

    //particle identifiers
    int id, stat; 

    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;

        int i = 0; // this is to count the 2 tau leptons 
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
                        //tauList.Add(particle);
                        if(i==0){tau1 = particle;}
                        if(i==1){tau2 = particle;}
                        i++;
                    }
                }

            } 
            if(id == 1000022){
                cout << particle << "  " << particle->getPDG() << endl; 
            }

        }//end for
        cout << "tau 1 " << tau1 << " " << tau1->getPDG() << endl;
        cout << "tau 2 " << tau2 << " " << tau2->getPDG() << endl;
        for(int particleIndex = 0; particleIndex < nElements ; particleIndex++){
            MCParticle* particle = dynamic_cast<MCParticle*>( col->getElementAt(particleIndex) );
            cout << endl;
            cout << endl;
            try{
                id = particle->getPDG();
                stat = particle->getGeneratorStatus();
            }
            catch(const std::exception& e){
                cout << "exception caught with message "<< e.what() << "\n";
            }
            double particle4Vector[4] = {particle->getEnergy(), particle->getMomentum()[0], 
                                         particle->getMomentum()[1], particle->getMomentum()[2]};
            //transform to R frame 
            double beta = (tau1->getEnergy() - tau2->getEnergy())/(tau1->getMomentum()[2] - tau2->getMomentum()[2]);
            cout << "BETA :"<< beta<<endl;
            double *R4Vector[4] = {Transform2RFrame( particle4Vector, beta )};
            cout << "Four Vec    " << endl;
            cout << particle4Vector[0] <<" "<<particle4Vector[1]<<" "<<particle4Vector[2]<<" "<<particle4Vector[3]<< endl;
            cout << "Transformed " << R4Vector[0]        <<" "<<R4Vector[1]       <<" "<<R4Vector[2]       <<" "<<R4Vector[3]<<" "<<R4Vector[4]<< endl; 
        } 
    }//ind if col

    _nEvt ++;

    cout << "event "<< _nEvt <<" finished " << endl;
}//end process

// function to transform into R frame 
double *SusyRazorVariables::Transform2RFrame(double in[4], double beta){
    cout << "----------------------------"<<endl;
    cout << "Running Transform Function!"<<endl;
    double beta2 = pow(beta,2);
    cout << "BETA SQUARED: "<<beta2<<endl; 
    double gamma = 1/(sqrt(1-pow(beta,2)));
    cout << "GAMMA: "<<gamma<<endl;
    cout << "In 4 Vec: " <<endl;
    cout << in[0] <<" "<<in[1]<<" "<<in[2]<<" "<<in[3]<<endl;
    double out[4] = {gamma*in[0]-gamma*beta*in[3], in[1], in[2], -gamma*beta*in[0]+gamma*in[3]}; 
    cout << "Out:"<<endl;
    cout << out[0]<<" "<<out[1]<<" "<<out[2]<<" "<<out[3]<<endl;
    cout << "----------------------------"<<endl; 
    return out;     
    
}

void SusyRazorVariables::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void SusyRazorVariables::end(){ 

    _rootfile->Write();
}
