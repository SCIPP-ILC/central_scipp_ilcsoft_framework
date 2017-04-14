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

    _rootfile = new TFile("WB_table.root","RECREATE");
    // usually a good idea to
    //printParameters() ;
    _prediction = new TH2F("predict", "Predicted Angle of Scatter, Correct vs Incorrect Kinematics", 1000, 0.0, 0.01, 1000, 0.0, 0.01);
    _table = new TH2F("table", "Reality vs Predicted Hit Status", 7, 0.0, 7.0, 7, 0.0, 7.0);
    _hit = new TH2F("hit", "Hit Status on BeamCal Face", 600, -30.0, 30.0, 600, -30.0, 30.0);
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

    int real, pred;

    srand(time(NULL));
    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
        cout << endl;
        cout << "************************EVENT: " << _nEvt << "*****************************" << endl;
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

                double tmag = sqrt(pow(mom[0], 2)+pow(mom[1], 2));
                mag+=tmag;
                
                double hyp = sqrt(pow(mom[0], 2)+pow(mom[1], 2)+pow(mom[2], 2));
                double cos = mom[2]/hyp;
                //if(cos<0.9){
                    //scipp_ilc::transform_to_lab(mom[0], mom[3], mom[0], mom[3]);
                    hadronic[0]+=mom[0];    
                    hadronic[1]+=mom[1];    
                    hadronic[2]+=mom[2];    
                    hadronic[3]+=mom[3];    
                //}
            }//end hadronic system    
        }//end for

        cout << "HADRONIC before correction: " << hadronic[0] << " " << hadronic[1] << endl;
        
        //create balancing particle
        double x = hadronic[0] + electronic[0];
        double y = hadronic[1] + electronic[1];

        double pseudo_x = -x;
        double pseudo_y = -y;

        //reset hadronic vector
        hadronic[0]=0;
        hadronic[1]=0;
        hadronic[2]=0;
        hadronic[3]=0;
        mag=0;
        
        //refill with detectable only  
        for(MCParticle* particle : final_system){
            id = particle->getPDG();
            //BEAM ELECTRON
            if(particle->getEnergy()==compEn_e){}
            else if(particle->getEnergy()==compEn_p){}
            else {
                if(abs(id)!=12&&abs(id)!=14&&abs(id)!=16&&abs(id)!=18){
                    mom[0]=particle->getMomentum()[0];    
                    mom[1]=particle->getMomentum()[1];    
                    mom[2]=particle->getMomentum()[2];
                    mom[3]=particle->getEnergy();
                    double tmag = sqrt(pow(mom[0], 2)+pow(mom[1], 2));
                    mag+=tmag;
                    
                    double hyp = sqrt(pow(mom[0], 2)+pow(mom[1], 2)+pow(mom[2], 2));
                    double cos = mom[2]/hyp;
                    //if(abs(cos)<0.9){
                        //scipp_ilc::transform_to_lab(mom[0], mom[3], mom[0], mom[3]);
                        hadronic[0]+=mom[0];    
                        hadronic[1]+=mom[1];    
                        hadronic[2]+=mom[2];    
                        hadronic[3]+=mom[3];    
                    //}
                }
            } 
        }

        cout << "HADRONIC: " << hadronic[0] << " " << hadronic[1] << endl;

        //add balancing particle to hadronic system
        hadronic[0]+=pseudo_x;
        hadronic[1]+=pseudo_y;

        //create prediction vectors
        pred_e[0] = -hadronic[0];
        pred_e[1] = -hadronic[1];
        pred_p[0] = -hadronic[0];
        pred_p[1] = -hadronic[1];

        double alpha = 500 - hadronic[3] - hadronic[2];
        double beta = 500 - hadronic[3] + hadronic[2];
        

        pred_e[2] = -(pow(pT, 2)-pow(alpha, 2))/(2*alpha);
        pred_e[3] = 500 - hadronic[3] - pred_e[2] - hadronic[2];
        pred_p[2] = (pow(eT, 2)-pow(beta, 2))/(2*beta);
        pred_p[3] = 500 - hadronic[3] + pred_p[2] + hadronic[2];
       
        //Lorentz transform
        scipp_ilc::transform_to_lab(real_e[0], real_e[3], real_e[0], real_e[3]);
        scipp_ilc::transform_to_lab(real_p[0], real_p[3], real_p[0], real_p[3]);
        scipp_ilc::transform_to_lab(pred_e[0], pred_e[3], pred_e[0], pred_e[3]);
        scipp_ilc::transform_to_lab(pred_p[0], pred_p[3], pred_p[0], pred_p[3]);
       
        //create position vector
        double real_e_pos[3];
        double real_p_pos[3];
        double pred_e_pos[3];
        double pred_p_pos[3];

        cout << "REAL ELECTRON: " << real_e[0] << " " << real_e[1] << " " << real_e[2] << endl;
        cout << "PRED ELECTRON: " << pred_e[0] << " " << pred_e[1] << " " << pred_e[2] << endl;

        real_e_pos[2] = scipp_ilc::_BeamCal_zmin;
        real_p_pos[2] = -scipp_ilc::_BeamCal_zmin;
        pred_e_pos[2] = scipp_ilc::_BeamCal_zmin;
        pred_p_pos[2] = -scipp_ilc::_BeamCal_zmin;

        real_e_pos[0] = real_e[0]*real_e_pos[2]/real_e[2];
        real_e_pos[1] = real_e[1]*real_e_pos[2]/real_e[2];
        real_p_pos[0] = real_p[0]*real_p_pos[2]/real_p[2];
        real_p_pos[1] = real_p[1]*real_p_pos[2]/real_p[2];
        pred_e_pos[0] = pred_e[0]*pred_e_pos[2]/pred_e[2];
        pred_e_pos[1] = pred_e[1]*pred_e_pos[2]/pred_e[2];
        pred_p_pos[0] = pred_p[0]*pred_p_pos[2]/pred_p[2];
        pred_p_pos[1] = pred_p[1]*pred_p_pos[2]/pred_p[2];
    
        //transform to BeamCal frame
        //scipp_ilc::z_to_beam_out(real_e_pos[0], real_e_pos[1], real_e_pos[2]);     
        //scipp_ilc::z_to_beam_out(real_p_pos[0], real_p_pos[1], real_p_pos[2]);     
        //scipp_ilc::z_to_beam_out(pred_e_pos[0], pred_e_pos[1], pred_e_pos[2]);     
        //scipp_ilc::z_to_beam_out(pred_p_pos[0], pred_p_pos[1], pred_p_pos[2]);    
        real_e_pos[0] = real_e_pos[0] - abs(real_e_pos[2]*scipp_ilc::_transform);
        real_p_pos[0] = real_p_pos[0] - abs(real_p_pos[2]*scipp_ilc::_transform);
        pred_e_pos[0] = pred_e_pos[0] - abs(pred_e_pos[2]*scipp_ilc::_transform);
        pred_p_pos[0] = pred_p_pos[0] - abs(pred_p_pos[2]*scipp_ilc::_transform);

        _hit->Fill(real_p_pos[0], real_p_pos[1]);
        //get hit status
        int re_hit = scipp_ilc::get_hitStatus(real_e_pos[0], real_p_pos[1]);
        int rp_hit = scipp_ilc::get_hitStatus(real_p_pos[0], real_p_pos[1]);
        int pe_hit = scipp_ilc::get_hitStatus(pred_e_pos[0], pred_p_pos[1]);
        int pp_hit = scipp_ilc::get_hitStatus(pred_p_pos[0], pred_p_pos[1]);


        //compare hit status
        if(re_hit==1 && rp_hit == 1){real=1;}
        else if (re_hit!=1 && rp_hit!=1){real=5;}
        else{real=3;}       
        
         
        if(pe_hit==1 && pp_hit == 1){pred=1;}
        else if (pe_hit!=1 && pp_hit!=1){pred=5;}
        else{pred=3;}  



        

        
        if(mag>1.0){ 
            total++;   
            if(pred==1 && real==1){hh_hh++;}
            if(pred==1 && real==3){hh_hm++;}
            if(pred==1 && real==5){hh_mm++;}
            if(pred==3 && real==1){hm_hh++;}
            if(pred==3 && real==3){hm_hm++;}
            if(pred==3 && real==5){hm_mm++;}
            if(pred==5 && real==1){mm_hh++;}
            if(pred==5 && real==3){mm_hm++;}
            if(pred==5 && real==5){mm_mm++;}
            _table->Fill(pred, real);     
        }


    }//end collection

    _nEvt ++ ;
}



void Pred_Compare::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void Pred_Compare::end(){
    cout << "TOTAL: " << total << endl; 
    double hh_hh_2=100.0*(double)hh_hh/(double)total;
    double hh_hm_2=100.0*(double)hh_hm/(double)total;
    double hh_mm_2=100.0*(double)hh_mm/(double)total;
    double hm_hh_2=100.0*(double)hm_hh/(double)total;
    double hm_hm_2=100.0*(double)hm_hm/(double)total;
    double hm_mm_2=100.0*(double)hm_mm/(double)total;
    double mm_hh_2=100.0*(double)mm_hh/(double)total;
    double mm_hm_2=100.0*(double)mm_hm/(double)total;
    double mm_mm_2=100.0*(double)mm_mm/(double)total;

    cout << "HH/HH: " << hh_hh_2 << endl;
    cout << "HH/HM: " << hh_hm_2 << endl;
    cout << "HH/MM: " << hh_mm_2 << endl;
    cout << "HM/HH: " << hm_hh_2 << endl;
    cout << "HM/HM: " << hm_hm_2 << endl;
    cout << "HM/MM: " << hm_mm_2 << endl;
    cout << "MM/HH: " << mm_hh_2 << endl;
    cout << "MM/HM: " << mm_hm_2 << endl;
    cout << "MM/MM: " << mm_mm_2 << endl;
    _rootfile->Write();
}
