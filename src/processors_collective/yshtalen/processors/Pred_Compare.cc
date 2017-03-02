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

#include "Pred_Compare.h"
#include "scipp_ilc_utilities.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/MCParticle.h>

#include <TFile.h>
#include <TH2D.h>
#include <MyParticle.h>
// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"



using namespace lcio;
using namespace marlin;
using namespace std;


Pred_Compare Pred_Compare;

static TFile* _rootfile;
static TH2F* _prediction;
static TH2F* _table;
static TH1F* _vector;

Pred_Compare::Pred_Compare() : Processor("Pred_Compare") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}



void Pred_Compare::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("BW_predicters_more.root","RECREATE");
    // usually a good idea to
    //printParameters() ;
    _prediction = new TH2F("predict", "Predicted Angle of Scatter, Correct vs Incorrect Kinematics", 1000, 0.0, 0.01, 1000, 0.0, 0.01);
    _table = new TH2F("table", "Reality vs Predicted Hit Status", 6, 0.0, 6.0, 6, 0.0, 6.0);
    _vector = new TH1F("vector", "Vector", 200, 0.0, 0.05);
    _nEvt = 0 ;

}



void Pred_Compare::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void Pred_Compare::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    LCCollection* col = evt->getCollection( _colName ) ;

    //vector<MyParticle*> particles;
    vector<MCParticle*> final_system;
    int stat, id =0;
    bool scatter;

    double compEn_e=0;
    double compEn_p=0;

    double mom[4];
    double real_e[4];
    double real_p[4];
    double pred_e[4];
    double pred_p[4];

    double eT, pT, mag;

    double hadronic[] = {0, 0, 0, 0};
    double electronic[] = {0, 0, 0, 0};

    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
        cout << endl;
        //cout << "************************EVENT: " << _nEvt << "*****************************" << endl;
        
        scatter = false;

//****************************************************INITIAL*****PASS******************************************************************************
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
            //cast to MCParticle
            MCParticle* particle = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );

            id = particle->getPDG();
            stat = particle->getGeneratorStatus();

            if(stat==1){
                //add to particle vector
                final_system.push_back(particle); 
                if(id==11){
                    if(particle->getEnergy() > compEn_e){compEn_e=particle->getEnergy();}    
                }
                else if(id==-11){
                    if(particle->getEnergy() > compEn_p){compEn_p=particle->getEnergy();}    
                }
                //set all other particle types to hadronic system
                else{
                    //particle->setHadronic(true);
                    //particle->setDetectable(true);
                }
            }//end final state   
        }//end for

//****************************************************SECOND*****PASS******************************************************************************
        for(MCParticle* particle : final_system){
            id = particle->getPDG();
            //BEAM ELECTRON
            if(particle->getEnergy()==compEn_e){
                real_e[0]=particle->getMomentum()[0];    
                real_e[1]=particle->getMomentum()[1];    
                real_e[2]=particle->getMomentum()[2];
                real_e[3]=particle->getEnergy();
                eT = sqrt(pow(real_e[0], 2)+pow(real_e[1], 2));
                //cout << "HEElectron: [" << real_e[0] << ", " << real_e[1] << ", " << real_e[2] << ", " << real_e[3] << "]" << endl; 
                if(abs(real_e[0])!=0||abs(real_e[1])!=0){
                    scatter = true;
                    //scipp_ilc::transform_to_lab(real_e[0], real_e[3], real_e[0], real_e[3]);
                    //cout << "HEElectron after transform: [" << real_e[0] << ", " << real_e[1] << ", " << real_e[2] << ", " << real_e[3] << "]" << endl; 
                    electronic[0]+=real_e[0];    
                    electronic[1]+=real_e[1];    
                    electronic[2]+=real_e[2];    
                    electronic[3]+=real_e[3];    
                }
            }//end unscattered
                
            //BEAM POSITRON
            else if(particle->getEnergy()==compEn_p){
                real_p[0]=particle->getMomentum()[0];    
                real_p[1]=particle->getMomentum()[1];    
                real_p[2]=particle->getMomentum()[2];
                real_p[3]=particle->getEnergy();
                pT = sqrt(pow(real_p[0], 2)+pow(real_p[1], 2));
                //cout << "Scattered: [" << real_p[0] << ", " << real_p[1] << ", " << real_p[2] << ", " << real_p[3] << "]" << endl; 
                if(abs(real_p[0])!=0||abs(real_p[1])!=0){
                    scatter = true;
                    //scipp_ilc::transform_to_lab(real_p[0], real_p[3], real_p[0], real_p[3]);
                    //cout << "Scattered after transform: [" << real_p[0] << ", " << real_p[1] << ", " << real_p[2] << ", " << real_p[3] << "]" << endl; 
                    electronic[0]+=real_p[0];    
                    electronic[1]+=real_p[1];    
                    electronic[2]+=real_p[2];    
                    electronic[3]+=real_p[3];    
                }    
            }//end scattered

            //HADRONIC SYSTEM
            else{
                mom[0]=particle->getMomentum()[0];    
                mom[1]=particle->getMomentum()[1];    
                mom[2]=particle->getMomentum()[2];
                mom[3]=particle->getEnergy();
                //scipp_ilc::transform_to_lab(mom[0], mom[3], mom[0], mom[3]);
                hadronic[0]+=mom[0];    
                hadronic[1]+=mom[1];    
                hadronic[2]+=mom[2];    
                hadronic[3]+=mom[3];    
                //cout << "Hadronic Particle ID: " << id <<" MOM [" << mom[0] << ", " << mom[1] << ", " << mom[2] << ", " << mom[3] << "]" << endl; 
                double tmag = sqrt(pow(mom[0], 2)+pow(mom[1], 2));
                mag+=tmag;
            }//end hadronic system    
        }//end for


        //create prediction vector
        double pred_e[4];
        double pred_p[4];
        pred_e[0] = -hadronic[0];
        pred_e[1] = -hadronic[1];
        pred_p[0] = -hadronic[0];
        pred_p[1] = -hadronic[1];

        double alpha = 500 - hadronic[3] - hadronic[2];
        double beta = 500 - hadronic[3] + hadronic[2];
        
        pred_e[2] = (pow(pT, 2)-pow(beta, 2))/(2*beta);
        pred_e[3] = 500 - hadronic[3] - pred_e[2] - hadronic[2];
        pred_p[2] = -(pow(eT, 2)-pow(alpha, 2))/(2*alpha);
        pred_p[3] = 500 - hadronic[3] + pred_p[2] + hadronic[2];
        

        
         
    }//end collection

    _nEvt ++ ;
}



void Pred_Compare::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void Pred_Compare::end(){ 
    cout << interest << endl;
    _rootfile->Write();
}
