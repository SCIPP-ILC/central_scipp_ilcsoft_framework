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
static TH1F* _prediction2;
static TH1F* _vector;
static TH1F* _z;
static TH1F* _eSum;
static TH2F* _rad;
static TH2F* _rad2;

Prediction::Prediction() : Processor("Prediction") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}



void Prediction::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("WW_CM_conserve.root","RECREATE");
    // usually a good idea to
    //printParameters() ;
    _prediction = new TH1F("predict", "Prediction Efficiency Angle, Transverse Mom Conserved", 500, 0.0, 0.005);
    _prediction2 = new TH1F("predict2", "Prediction Efficiency Angle, Transverse Mom Not Conserved", 500, 0.0, 0.005);
    _vector = new TH1F("vector", "Vector", 200, 0.0, 0.05);
    _z = new TH1F("z", "Z-mom Hadronic", 5000, -250.0, 250.0);
    _rad = new TH2F("rad", "Prediction Efficiency Due to Initial State Radiation", 500, 0.0, 0.01, 2000, 200.0, 500.0);
    _rad2 = new TH2F("rad2", "Prediction Efficiency Due to Initial State Radiation", 500, 0.0, 0.01, 2000, 200.0, 500.0);
    _eSum = new TH1F("eSum", "Sum of Electron/Positron Energies", 1000, 200.0, 500.0);
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

    double tmom, theta, thetau, mag, eT, pT;
    
    bool scatter;

    double hadronic[] = {0, 0, 0, 0};
    double electronic[] = {0, 0, 0, 0};
    double unscat[] = {0, 0, 0, 0};

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


            /*cout << "Particle " << hitIndex << " with stat: " << stat << " with ID: " << id; 
            cout << " with mom " << particle->getMomentum()[0] << ", " << particle->getMomentum()[1] << ", " << particle->getMomentum()[2]; 
            cout << " with energy: " << particle->getEnergy() << endl;*/
            //cut on final state
            if(stat==1){
                //add to particle vector
                final_system.push_back(particle); 
                /*cout << "Particle " << hitIndex << " with ID: " << id; 
                cout << " with mom " << particle->getMomentum()[0] << ", " << particle->getMomentum()[1] << ", " << particle->getMomentum()[2]; 
                cout << " with energy: " << particle->getEnergy() << endl;
                */
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
//****************************************************SECOND*****PASS******************************************************************************
        for(MCParticle* particle : final_system){
            id = particle->getPDG();
            //UNSCATTERED BEAM PARTICLE
            if(particle->getEnergy()==compEn_p){
                mom_e[0]=particle->getMomentum()[0];    
                mom_e[1]=particle->getMomentum()[1];    
                mom_e[2]=particle->getMomentum()[2];
                mom_e[3]=particle->getEnergy();
                eT = sqrt(pow(mom_e[0], 2)+pow(mom_e[1], 2));
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
                    unscat[0]=mom_e[0];
                    unscat[1]=mom_e[1];
                    unscat[2]=mom_e[2];
                    //scipp_ilc::transform_to_lab(mom_e[0], mom_e[3], mom_e[0], mom_e[3]);
                    //cout << "HEElectron after transform: [" << mom_e[0] << ", " << mom_e[1] << ", " << mom_e[2] << ", " << mom_e[3] << "]" << endl; 
                }   
            }//end unscattered    
            //SCATTERED BEAM PARTICLE
            else if(particle->getEnergy()==compEn_e){
                mom_p[0]=particle->getMomentum()[0];    
                mom_p[1]=particle->getMomentum()[1];    
                mom_p[2]=particle->getMomentum()[2];
                mom_p[3]=particle->getEnergy();
                pT = sqrt(pow(mom_p[0], 2)+pow(mom_p[1], 2));
                //cout << "Scattered: [" << mom_p[0] << ", " << mom_p[1] << ", " << mom_p[2] << ", " << mom_p[3] << "]" << endl; 
                if(abs(mom_p[0])!=0||abs(mom_p[1])!=0){
                    scatter = true;
                    //scipp_ilc::transform_to_lab(mom_p[0], mom_p[3], mom_p[0], mom_p[3]);
                    //cout << "Scattered after transform: [" << mom_p[0] << ", " << mom_p[1] << ", " << mom_p[2] << ", " << mom_p[3] << "]" << endl; 
                    electronic[0]+=mom_p[0];    
                    electronic[1]+=mom_p[1];    
                    electronic[2]+=mom_p[2];    
                    electronic[3]+=mom_p[3];    
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
        _z->Fill(hadronic[2]);
        cout << endl;
        cout << endl;
        //cout << "Hadronic Vector: [" << hadronic[0] << ", " << hadronic[1] << ", " << hadronic[2] << ", " << hadronic[3] << "]" << endl; 
        
        if(scatter == true){

            //create prediction vector
            double predict[3];
            predict[0] = -hadronic[0];
            predict[1] = -hadronic[1];
            double alpha = 500 - hadronic[3] - hadronic[2];
            double beta = 500 - hadronic[3] + hadronic[2];
            predict[2] = -(pow(eT, 2)-pow(alpha, 2))/(2*alpha);
            //predict[2] = (pow(eT, 2)-pow(beta, 2))/(2*beta);

            double upredict[3];
            upredict[0]=unscat[0];
            upredict[1]=unscat[1];
            upredict[2]=-hadronic[2]-electronic[2];




            //cout << "Electronic Vector: [" << electronic[0] << ", " << electronic[1] << ", " << electronic[2] << "]" << endl;
            //cout << "Prediction Vector: [" << predict[0] << ", " << predict[1] << ", " << predict[2] << "]" << endl;

            double dot = electronic[0]*predict[0] + electronic[1]*predict[1] + electronic[2]*predict[2];
            double e_mag = sqrt(pow(electronic[0], 2)+pow(electronic[1], 2)+pow(electronic[2], 2)); 
            double p_mag = sqrt(pow(predict[0], 2)+pow(predict[1], 2)+pow(predict[2], 2)); 
            theta = acos(dot/(e_mag*p_mag)); 
            _prediction2->Fill(theta);
            cout << endl;
            //cout << "Scattered Prediction Efficiency: " << theta << endl;
        
            double x = hadronic[0]+electronic[0];
            double pseudo_x = -x; 
            double y = hadronic[1]+electronic[1];
            double pseudo_y = -y;
            double vector = sqrt(pow(x, 2)+pow(y, 2));           
            //cout << "V: " << vector << endl;

            //if(vector<0.001){

                    cout << endl;
                    cout << endl;
                //if(theta>0.003){
                    cout << "************************EVENT: " << _nEvt << "*****************************" << endl;
                    cout << "SUM VECTOR MAGNITUDE, V = " << vector << endl;
                    cout << "Electronic Vector: [" << electronic[0] << ", " << electronic[1] << ", " << electronic[2] << "]" << endl;
                    cout << "Prediction Vector: [" << predict[0] << ", " << predict[1] << ", " << predict[2] << "]" << endl;
                    cout << "Theta Difference: " << theta << endl;

                    double eTot = 0;
                    double zMom = 0;
                    for(MCParticle* particle : final_system){
                        eTot+=particle->getEnergy();
                        zMom+=particle->getMomentum()[2];


                    }
                    cout << "Total Event Energy: " << eTot << endl;
                    cout << "Total Event Z-Momentum: " << zMom << endl;
                    hadronic[0]+=pseudo_x;
                    hadronic[1]+=pseudo_y;
                    x = hadronic[0]+electronic[0];
                    y = hadronic[1]+electronic[1];
                    vector = sqrt(pow(x, 2)+pow(y,2));
                    cout << "****** After Adjustment *****" << endl;    
                    _vector->Fill(vector);
                    cout << "SUM VECTOR MAGNITUDE, V = " << vector << endl;
                    //double En = 500.0-eTot;
                    //create prediction vector
                    predict[0] = -hadronic[0];
                    predict[1] = -hadronic[1];
                    //hadronic[3]+=En;
                    alpha = 500 - hadronic[3] - hadronic[2];
                    beta = 500 - hadronic[3] + hadronic[2];
                    predict[2] = -(pow(eT, 2)-pow(alpha, 2))/(2*alpha);
                    //predict[2] = (pow(eT, 2)-pow(beta, 2))/(2*beta);
                    
                    //get prediction efficiency angle
                    dot = electronic[0]*predict[0] + electronic[1]*predict[1] + electronic[2]*predict[2];
                    e_mag = sqrt(pow(electronic[0], 2)+pow(electronic[1], 2)+pow(electronic[2], 2)); 
                    p_mag = sqrt(pow(predict[0], 2)+pow(predict[1], 2)+pow(predict[2], 2)); 
                    theta = acos(dot/(e_mag*p_mag)); 
                    if(theta<0.000001){
                        interest++;
                        for(MCParticle* particle : final_system){ 
                            id = particle->getPDG();
                            cout << "Particle ID: " << id; 
                            cout << " with mom " << particle->getMomentum()[0] << ", " << particle->getMomentum()[1] << ", " << particle->getMomentum()[2]; 
                            cout << " with energy: " << particle->getEnergy() << endl;
                        }
                    }
                    cout << "Electronic Vector: [" << electronic[0] << ", " << electronic[1] << ", " << electronic[2] << "]" << endl;
                    cout << "Prediction Vector: [" << predict[0] << ", " << predict[1] << ", " << predict[2] << "]" << endl;
                    cout << "Theta Difference: " << theta << endl;
                    _prediction->Fill(theta);
                    double eSum = abs(mom_e[3]+mom_p[3]+hadronic[3]);
                    _eSum->Fill(eSum);
                    if(theta>0.001){
                        _rad->Fill(theta, eSum);
                    }
                    else{_rad2->Fill(theta, eSum);}
                //}
            //}

            cout << endl;
            //cout << "Unscattered Actual: " << unscat[0] << " " << unscat[1] << " " << unscat[2] << endl;
            //cout << "Unscattered Predict: " << upredict[0] << " " << upredict[1] << " " << upredict[2] << endl;
            double diff = abs(upredict[2]-unscat[2]);
            //cout << "Unscattered Prediction Efficiency: " << diff << endl;
            cout << endl;
            
            


            
        }
         
    }//end collection

    _nEvt ++ ;
}



void Prediction::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void Prediction::end(){ 
    cout << interest << endl;
    _rootfile->Write();
}
