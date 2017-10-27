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
 * author Y. Shtalenkova
 * April 5, 2016
 */

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

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
static TH2F* _hit;
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

    _rootfile = new TFile("BW_table.root","RECREATE");
    // usually a good idea to
    //printParameters() ;
    _prediction = new TH2F("predict", "Predicted Angle of Scatter, Correct vs Incorrect Kinematics", 1000, 0.0, 0.01, 1000, 0.0, 0.01);
    _table = new TH2F("table", "Reality vs Predicted Hit Status", 7, 0.0, 7.0, 7, 0.0, 7.0);
    _hit = new TH2F("hit", "Hit Status on BeamCal Face", 6000, -150.0, 150.0, 6000, -150.0, 150.0);
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
    vector<MCParticle*> had;
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

    int real, pred;

    srand(time(NULL));
    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
        //cout << endl;
        //cout << "************************EVENT: " << _nEvt << "*****************************" << endl;
        scatter = false;

//****************************************************INITIAL*****PASS******************************************************************************
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
            MCParticle* particle = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );

            id = particle->getPDG();
            stat = particle->getGeneratorStatus();

            if(stat==1){
                //push all final state particles to separate vector
                final_system.push_back(particle); 
                if(id==11){
                    if(particle->getEnergy() > compEn_e){compEn_e=particle->getEnergy();}    
                }
                else if(id==-11){
                    if(particle->getEnergy() > compEn_p){compEn_p=particle->getEnergy();}    
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
                
                //if deflected, add to electronic
                if(abs(real_e[0])!=0||abs(real_e[1])!=0){
                    scatter = true;
                    electronic[0]+=real_e[0];    
                    electronic[1]+=real_e[1];    
                    electronic[2]+=real_e[2];    
                    electronic[3]+=real_e[3];    
                }
            }
                
            //BEAM POSITRON
            else if(particle->getEnergy()==compEn_p){
                real_p[0]=particle->getMomentum()[0];    
                real_p[1]=particle->getMomentum()[1];    
                real_p[2]=particle->getMomentum()[2];
                real_p[3]=particle->getEnergy();
                pT = sqrt(pow(real_p[0], 2)+pow(real_p[1], 2));
                
                //if deflected, add to electronic
                if(abs(real_p[0])!=0||abs(real_p[1])!=0){
                    scatter = true;
                    electronic[0]+=real_p[0];    
                    electronic[1]+=real_p[1];    
                    electronic[2]+=real_p[2];    
                    electronic[3]+=real_p[3];    
                }    
            }

            //HADRONIC SYSTEM
            else{
                had.push_back(particle);
                mom[0]=particle->getMomentum()[0];    
                mom[1]=particle->getMomentum()[1];    
                mom[2]=particle->getMomentum()[2];
                mom[3]=particle->getEnergy();

                double tmag = sqrt(pow(mom[0], 2)+pow(mom[1], 2));
                mag+=tmag;
                
                //scipp_ilc::transform_to_lab(mom[0], mom[3], mom[0], mom[3]);
                hadronic[0]+=mom[0];    
                hadronic[1]+=mom[1];    
                hadronic[2]+=mom[2];    
                hadronic[3]+=mom[3];    
            }//end hadronic system    
        }//end for

        if(scatter==true){
            


	  //cout << "HADRONIC before correction: " << hadronic[0] << " " << hadronic[1] << endl;
            
            //create balancing particle
            double x = hadronic[0] + electronic[0];
            double y = hadronic[1] + electronic[1];

            double pseudo_x = -x;
            double pseudo_y = -y;

            double p_mag = sqrt(pow(pseudo_x, 2)+pow(pseudo_y, 2));
            mag+=p_mag;

            hadronic[0]=0;
            hadronic[1]=0;

            //refill hadronic system with cuts and balancing particle
            for(MCParticle* particle : had){
                mom[0]=particle->getMomentum()[0];    
                mom[1]=particle->getMomentum()[1];    
                mom[2]=particle->getMomentum()[2];
                mom[3]=particle->getEnergy();

                double hyp = sqrt(pow(mom[0], 2)+pow(mom[1], 2)+pow(mom[2], 2));
                double cos = mom[2]/hyp;
                if(abs(cos)<0.5){
                    //scipp_ilc::transform_to_lab(mom[0], mom[3], mom[0], mom[3]);
                    hadronic[0]+=mom[0];    
                    hadronic[1]+=mom[1];    
                    hadronic[2]+=mom[2];    
                    hadronic[3]+=mom[3];    
                }
            }

            //add balancing particle to hadronic system
            hadronic[0]+=pseudo_x;
            hadronic[1]+=pseudo_y;

            //cut on S
            if(mag>1.0){
                
                total++;
                //create prediction vectors
                pred_e[0] = -hadronic[0];
                pred_e[1] = -hadronic[1];
                pred_p[0] = -hadronic[0];
                pred_p[1] = -hadronic[1];

                double alpha = 500 - hadronic[3] - hadronic[2];
                double beta = 500 - hadronic[3] + hadronic[2];
                

                pred_e[2] = -(pow(eT, 2)-pow(alpha, 2))/(2*alpha);
                pred_e[3] = 500 - hadronic[3] - pred_e[2] - hadronic[2];
                pred_p[2] = (pow(pT, 2)-pow(beta, 2))/(2*beta);
                pred_p[3] = 500 - hadronic[3] + pred_p[2] + hadronic[2];
               
                //Lorentz transform - frome center of mass to lab frame 
                scipp_ilc::transform_to_lab(real_e[0], real_e[3], real_e[0], real_e[3]);
                scipp_ilc::transform_to_lab(real_p[0], real_p[3], real_p[0], real_p[3]);
                scipp_ilc::transform_to_lab(pred_e[0], pred_e[3], pred_e[0], pred_e[3]);
                scipp_ilc::transform_to_lab(pred_p[0], pred_p[3], pred_p[0], pred_p[3]);
               
                //create position vector
                double real_e_pos[3];
                double real_p_pos[3];
                double pred_e_pos[3];
                double pred_p_pos[3];



                //set z-positions as beamcal face
                real_e_pos[2] = 3265;
                real_p_pos[2] = -3265;
                pred_e_pos[2] = 3265;
                pred_p_pos[2] = -3265;

                //extrapolate transverse positions from mom vector
                real_e_pos[0] = real_e[0]*real_e_pos[2]/real_e[2];
                real_e_pos[1] = real_e[1]*real_e_pos[2]/real_e[2];
                real_p_pos[0] = real_p[0]*real_p_pos[2]/real_p[2];
                real_p_pos[1] = real_p[1]*real_p_pos[2]/real_p[2];
                
                pred_e_pos[0] = pred_e[0]*pred_e_pos[2]/pred_e[2];
                pred_e_pos[1] = pred_e[1]*pred_e_pos[2]/pred_e[2];
                pred_p_pos[0] = pred_p[0]*pred_p_pos[2]/pred_p[2];
                pred_p_pos[1] = pred_p[1]*pred_p_pos[2]/pred_p[2];
           
                 
                //transform to BeamCal frame
                real_e_pos[0] = real_e_pos[0] - real_e_pos[2]*0.007;
                real_p_pos[0] = real_p_pos[0] + real_p_pos[2]*0.007;
                pred_e_pos[0] = pred_e_pos[0] - pred_e_pos[2]*0.007;
                pred_p_pos[0] = pred_p_pos[0] + pred_p_pos[2]*0.007;
                

                //get hit status
                int re_hit = scipp_ilc::get_hitStatus(real_e_pos[0], real_e_pos[1], real_e_pos[2]);
                int rp_hit = scipp_ilc::get_hitStatus(real_p_pos[0], real_p_pos[1], real_p_pos[2]);
                int pe_hit = scipp_ilc::get_hitStatus(pred_e_pos[0], pred_e_pos[1], pred_e_pos[2]);
                int pp_hit = scipp_ilc::get_hitStatus(pred_p_pos[0], pred_p_pos[1], pred_p_pos[2]);

               /* 
                //eWpB 
                if(re_hit!=3 && re_hit!=4){real = 1;}
                else{real = 2;}
                if(pe_hit!=3 && pe_hit!=4){pred = 1;}
                else{pred = 2;}

                if(pred == 1 && real == 1){hh++;}
                else if(pred == 1 && real == 2){hm++;}
                else if(pred == 2 && real == 1){mh++;}
                else if(pred == 2 && real == 2){mm++;}
                

                cout << "REAL ELECTRON: " << real_e_pos[0] << " " << real_e_pos[1] << " " << real_e_pos[2] << endl;
                cout << "PRED ELECTRON: " << pred_e_pos[0] << " " << pred_e_pos[1] << " " << pred_e_pos[2] << endl;
             
             */
                //eBpW 
                if(rp_hit!=3 && rp_hit!=4){real = 1;}
                else{real = 2;}
                if(pp_hit!=3 && pp_hit!=4){pred = 1;}
                else{pred = 2;}

                if(pred == 1 && real == 1){hh++;}
                else if(pred == 1 && real == 2){hm++;}
                else if(pred == 2 && real == 1){mh++;}
                else if(pred ==2 && real == 2){mm++;}

                //cout << "REAL POSITRON: " << real_p_pos[0] << " " << real_p_pos[1] << " " << real_p_pos[2] << endl;
                //cout << "PRED POSITRON: " << pred_p_pos[0] << " " << pred_p_pos[1] << " " << pred_p_pos[2] << endl;
                
                real = 0;
                pred = 0;
            }


        }//end if scatter
    }//end collection

    _nEvt ++ ;
}



void Pred_Compare::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void Pred_Compare::end(){
    cout << "TOTAL: " << total << endl; 
    double hh2=(double)hh/(double)total;
    double hm2=(double)hm/(double)total;
    double mh2=(double)mh/(double)total;
    double mm2=(double)mm/(double)total;

    cout << "HH: " << hh2 << endl;
    cout << "HM: " << hm2 << endl;
    cout << "MH: " << mh2 << endl;
    cout << "MM: " << mm2 << endl;
    _rootfile->Write();
}
