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

#include "Tester.h"
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

Tester Tester;

static TFile* _rootfile;
static TH2F* _pos;
static TH2F* _angle;

Tester::Tester() : Processor("Tester") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}



void Tester::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("WW_ang.root","RECREATE");
    _pos = new TH2F("pos", "Distribution of Deflected Beam Particles at Position of BeamCal Face", 900, -150.0, 150.0, 900, -150.0, 150.0);
    _angle = new TH2F("ang", "Positron vs Electron Deflection Angle", 100, 0.0, 0.02, 100, 0.0, 0.02);

    
    // usually a good idea to
    //printParameters() ;
    _nEvt = 0 ;

}



void Tester::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void Tester::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...
    double tot[]={0, 0, 0};
    LCCollection* col = evt->getCollection( _colName ) ;

    int stat, id = 0;
    double theta_e, theta_p;
    double E_e, E_p = 0;
    double mom_e[4], mom_p[4];
    double pos_e[3], pos_p[3];

    int hit = 0;

    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
    cout << endl;
    cout << endl;
    cout << "***************************EVENT: " << _nEvt << "****************************" << endl;

    vector<MCParticle*> system;
        
        //create final state subsystem
        //determine beam particle energies for identification
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

                double mom[4];
                mom[0] = hit->getMomentum()[0]; 
                mom[1] = hit->getMomentum()[1]; 
                mom[2] = hit->getMomentum()[2];
                mom[3] = hit->getEnergy();

                cout << " momentum [" << mom[0] << ", " << mom[1] << ", " << mom[2] << "] with energy: " << mom[3] << endl;

            }//end final state
        }//end for

        //find beam particles
        for(MCParticle* particle : system){
            id = particle->getPDG();
            if(id==11&&particle->getEnergy()==E_e){
                    const double* mom = particle->getMomentum();

                    mom_e[0]=mom[0];
                    mom_e[1]=mom[1];
                    mom_e[2]=mom[2];
                    mom_e[3]=particle->getEnergy();
                    if(abs(mom_e[0])!=0||abs(mom_e[1])!=0){
                        hit++;
                        double r = sqrt(pow(mom_e[0], 2)+pow(mom_e[1], 2));
                        double mag = sqrt(pow(mom_e[0], 2)+pow(mom_e[1], 2)+pow(mom_e[2], 2));
                        theta_e = asin(r/mag);
                    }//end scatter
            }//end high energy electron    
            else if(id==-11&&particle->getEnergy()==E_p){
                    const double* mom = particle->getMomentum();
                    
                    mom_p[0] = mom[0]; 
                    mom_p[1] = mom[1]; 
                    mom_p[2] = mom[2];
                    mom_p[3] = particle->getEnergy();
                    if(abs(mom_p[0])!=0||abs(mom_p[1])!=0){
                        hit++;
                        double r = sqrt(pow(mom_p[0], 2)+pow(mom_p[1], 2));
                        double mag = sqrt(pow(mom_p[0], 2)+pow(mom_p[1], 2)+pow(mom_p[2], 2));
                        theta_p = asin(r/mag);
                    }//end scatter
            }//end high energy electron    
            
            else{
                double mom[4];
                mom[0] = particle->getMomentum()[0]; 
                mom[1] = particle->getMomentum()[1]; 
                mom[2] = particle->getMomentum()[2];
                mom[3] = particle->getEnergy();
                double mag = sqrt(pow(mom[0], 2)+pow(mom[1], 2));    
            }
        }//end for
        if(hit==2){
            _angle->Fill(theta_e, theta_p);
        }

        
        scipp_ilc::transform_to_lab(mom_e[0], mom_e[3], mom_e[0], mom_e[3]);
        scipp_ilc::transform_to_lab(mom_p[0], mom_p[3], mom_p[0], mom_p[3]);

        /*pos_e[0] = mom_e[0]*scipp_ilc::_BeamCal_zmin/mom_e[2];
        pos_e[1] = mom_e[1]*scipp_ilc::_BeamCal_zmin/mom_e[2];
        pos_e[2] = scipp_ilc::_BeamCal_zmin;
        pos_p[0] = mom_p[0]*scipp_ilc::_BeamCal_zmin/mom_p[2];
        pos_p[1] = mom_p[1]*scipp_ilc::_BeamCal_zmin/mom_p[2];
        pos_p[2] = scipp_ilc::_BeamCal_zmin;

        scipp_ilc::z_to_beam_out(pos_p[0], pos_p[1], pos_p[2]);
        scipp_ilc::z_to_beam_out(pos_e[0], pos_e[1], pos_e[2]);

        _pos->Fill(pos_p[0], pos_p[1]);
        _pos->Fill(pos_e[0], pos_e[1]);
        */
    }

    _nEvt ++ ;
}



void Tester::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void Tester::end(){ 
    _rootfile->Write();
}
