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
#include <MyParticle.h>
// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"



using namespace lcio;
using namespace marlin;
using namespace std;


Prediction Prediction;

static TFile* _rootfile;
static TH1F* _prediction;

Prediction::Prediction() : Processor("Prediction") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}



void Prediction::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("predict_test.root","RECREATE");
    // usually a good idea to
    //printParameters() ;
    _prediction = new TH1F("predict", "Prediction", 200, 0.0, 0.05);
    _nEvt = 0 ;

}



void Prediction::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void Prediction::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    LCCollection* col = evt->getCollection( _colName ) ;

    //vector<MyParticle*> particles;
    vector<MCParticle*> final_system;
    int stat, id =0;
    double tot_mom[]={0, 0};
    double compEn_e=0;
    double compEn_p=0;

    double mom[4];
    double mom_e[4];
    double mom_p[4];

    double tmom, theta, mag;

    bool scatter;

    double hadronic[] = {0, 0, 0, 0};
    double electronic[] = {0, 0, 0, 0};

    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
        cout << endl;
        cout << endl;
        cout << endl;
        cout << endl;
        cout << "************************EVENT: " << _nEvt << "*****************************" << endl;
        cout << "****************************** " << _nEvt << "*****************************" << endl;
        
        scatter = false;

        //INITIAL CYCLE THROUGH COLLECTION
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
            //cast to MCParticle and create MyParticle object for each
            MCParticle* particle = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
            //cout << "MCParticle* cast, typeof: " << typeid(hit).name() << endl;
            //MyParticle particle;
            //particle.setSource(hit);

            id = particle->getPDG();
            stat = particle->getGeneratorStatus();


            /*cout << "Particle " << hitIndex << " with stat: " << stat << " with ID: " << id; 
            cout << " with mom " << particle->getMomentum()[0] << ", " << particle->getMomentum()[1] << ", " << particle->getMomentum()[2]; 
            cout << " with energy: " << particle->getEnergy() << endl;*/
            //cut on final state
            if(stat==1){
                //add to particle vector
                final_system.push_back(particle); 
                cout << "Particle " << hitIndex << " with ID: " << id; 
                cout << " with mom " << particle->getMomentum()[0] << ", " << particle->getMomentum()[1] << ", " << particle->getMomentum()[2]; 
                cout << " with energy: " << particle->getEnergy() << endl;

                //find highest energy electron and positron
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

        cout << endl;
        cout << endl;
        //SECOND PASS THROUGH FINAL STATE PARTICLES ONLY
        /*for(MyParticle* particle : particles){
            //assign electronic system e+/e-
            if(id==11 && particle->energy()==compEn_e){particle->setElectronic(true);}    
            else if(id==-11 && particle->energy()==compEn_p){particle->setElectronic(true);}
            //all rest of particles are hadronic
            else{
                //all other particles are hadronic 
                particle->setHadronic(true);
                //exclude neutrinos
                if(abs(id)==12||abs(id)==14||abs(id)==16||abs(id)==18){particle->setDetectable(false);}    
            }
        }

        for(MyParticle* particle : particles){
            if(particle->isElectronic()){
                cout << "ELECTRONIC:: id: " << particle->id() << " Energy: " << particle->energy() << endl;
            }    
            if(particle->isHadronic()){
                cout << "HADRONIC:: id: " << particle->id() << " Energy: " << particle->energy() << endl;
            }
        }*/
        //cout << "Highest energies: " << compEn_e << ", " << compEn_p << endl;

        for(MCParticle* particle : final_system){
            id = particle->getPDG();

            if(particle->getEnergy()==compEn_p){
                mom_e[0]=particle->getMomentum()[0];    
                mom_e[1]=particle->getMomentum()[1];    
                mom_e[2]=particle->getMomentum()[2];
                mom_e[3]=particle->getEnergy();
                //cout << "HEElectron: [" << mom_e[0] << ", " << mom_e[1] << ", " << mom_e[2] << ", " << mom_e[3] << "]" << endl; 
                if(abs(mom_e[0])!=0||abs(mom_e[1])!=0){
                    scatter = true;
                    //scipp_ilc::transform_to_lab(mom_e[0], mom_e[3], mom_e[0], mom_e[3]);
                    //cout << "HEElectron after transform: [" << mom_e[0] << ", " << mom_e[1] << ", " << mom_e[2] << ", " << mom_e[3] << "]" << endl; 
                    electronic[0]+=mom_e[0];    
                    electronic[1]+=mom_e[1];    
                    electronic[2]+=mom_e[2];    
                    electronic[3]+=mom_e[3];    
                }
                else{
                    //scipp_ilc::transform_to_lab(mom_e[0], mom_e[3], mom_e[0], mom_e[3]);
                    //cout << "HEElectron after transform: [" << mom_e[0] << ", " << mom_e[1] << ", " << mom_e[2] << ", " << mom_e[3] << "]" << endl; 
                }   
            }//end beam electron    
            else if(particle->getEnergy()==compEn_e){
                mom_p[0]=particle->getMomentum()[0];    
                mom_p[1]=particle->getMomentum()[1];    
                mom_p[2]=particle->getMomentum()[2];
                mom_p[3]=particle->getEnergy();
                //cout << "HEPositron: [" << mom_p[0] << ", " << mom_p[1] << ", " << mom_p[2] << ", " << mom_p[3] << "]" << endl; 
                if(abs(mom_p[0])!=0||abs(mom_p[1])!=0){
                    scatter = true;
                    //cout << "----------------------------SCATTERED POSITRON------------------------------" << endl;
                    //scipp_ilc::transform_to_lab(mom_p[0], mom_p[3], mom_p[0], mom_p[3]);
                    //cout << "HEPositron after transform: [" << mom_p[0] << ", " << mom_p[1] << ", " << mom_p[2] << ", " << mom_p[3] << "]" << endl; 
                    electronic[0]+=mom_p[0];    
                    electronic[1]+=mom_p[1];    
                    electronic[2]+=mom_p[2];    
                    electronic[3]+=mom_p[3];    
                }    
            }//end beam positron
            else{
                mom[0]=particle->getMomentum()[0];    
                mom[1]=particle->getMomentum()[1];    
                mom[2]=particle->getMomentum()[2];
                mom[3]=particle->getEnergy();
                scipp_ilc::transform_to_lab(mom[0], mom[3], mom[0], mom[3]);
                hadronic[0]+=mom[0];    
                hadronic[1]+=mom[1];    
                hadronic[2]+=mom[2];    
                hadronic[3]+=mom[3];    
                //cout << "Hadronic Particle ID: " << id <<" MOM [" << mom[0] << ", " << mom[1] << ", " << mom[2] << ", " << mom[3] << "]" << endl; 
                double tmag = sqrt(pow(mom[0], 2)+pow(mom[1], 2));
                mag+=tmag;
            }//end hadronic system    
        }//end for

        cout << endl;
        cout << endl;
        //cout << "Hadronic Vector: [" << hadronic[0] << ", " << hadronic[1] << ", " << hadronic[2] << ", " << hadronic[3] << "]" << endl; 
        
        if(scatter == true){
            //create prediction vector
            double predict[3];
            predict[0] = -hadronic[0];
            predict[1] = -hadronic[1];
            predict[2] = sqrt(pow(250.0-hadronic[3], 2) - pow(mag, 2));
            
            //cout << "Electronic Vector: [" << electronic[0] << ", " << electronic[1] << ", " << electronic[2] << "]" << endl;
            //cout << "Prediction Vector: [" << predict[0] << ", " << predict[1] << ", " << predict[2] << "]" << endl;

            double dot = electronic[0]*predict[0] + electronic[1]*predict[1] + electronic[2]*predict[2];
            double e_mag = sqrt(pow(electronic[0], 2)+pow(electronic[1], 2)+pow(electronic[2], 2)); 
            double p_mag = sqrt(pow(predict[0], 2)+pow(predict[1], 2)+pow(predict[2], 2)); 
            theta = acos(dot/(e_mag*p_mag));
            //cout << "Angle between vectors: " << theta << endl;           
            _prediction->Fill(theta);
            
            //if(theta>0.005){
                //cout << "****************************EVENT " << _nEvt << "*******************************************" << endl;
                //cout << "********************************************************************************************" << endl;
                for(MCParticle* particle : final_system){ 
                    double x, e;
                    //scipp_ilc::transform_to_lab(particle->getMomentum()[0], particle->getEnergy(), x, e);   
                    /*cout << "Particle with ID: " << particle->getPDG(); 
                    cout << " with mom " << x << ", " << particle->getMomentum()[1] << ", " << particle->getMomentum()[2]; 
                    cout << " with energy: " << e << endl;*/

                }
                cout << endl;
                cout << endl;
                cout << "Hadronic Vector: [" << hadronic[0] << ", " << hadronic[1] << ", " << hadronic[2] << ", " << hadronic[3] << "]" << endl; 
                double h_trans = sqrt(pow(hadronic[0], 2)+pow(hadronic[1], 2));
                cout << "Hadronic Transverse Momentum Magnitude: " << h_trans << endl;
                cout << "Electronic Vector: [" << electronic[0] << ", " << electronic[1] << ", " << electronic[2] << "]" << endl;
                double e_trans = sqrt(pow(electronic[0], 2)+pow(electronic[1], 2));
                cout << "Electronic Transverse Momentum Magnitude: " << e_trans << endl;
                cout << "Prediction Vector: [" << predict[0] << ", " << predict[1] << ", " << predict[2] << "]" << endl;
                cout << "Angle between vectors: " << theta << endl;           
            //}       
        }
         
    }//end collection

    _nEvt ++ ;
}



void Prediction::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void Prediction::end(){ 
    _rootfile->Write();
}
