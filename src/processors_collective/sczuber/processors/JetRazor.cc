#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0
/* 
 * author Summer Zuber 
 * January 18, 2017 
 *
 * Including fastjet/ClusterSequence.hh to use fastjet algorithms to identify jets. Uses those jets to make Razor Variables.  
 */

#include "JetRazor.h"
#include "scipp_ilc_utilities.h"

#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/MCParticle.h>

#include <TFile.h>
#include <TH2D.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio> // trying to creat log file
//#include <map> 

#include <CLHEP/Vector/ThreeVector.h>
// #include <CLHEP/Random/RanluxEngine.h>

#include <EVENT/LCCollection.h>
#include <EVENT/LCIO.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/LCTOOLS.h>

#include "megajet.h"

using namespace lcio;
using namespace marlin;
using namespace fastjet;
using namespace std;

static TFile* _rootfile;

static TH1F* _R_T;
static TH1F* _R_DAB;
static TH1F* _R_DED;
static TH1F* _MR_T;
static TH1F* _MR_DAB;
static TH1F* _MR_DED; 

static TH1F* _NJ_T; // histogram of the number of jets using TRUE jets
static TH1F* _NJ_DAB; // histogram of the number of jets using DAB jets
static TH1F* _NJ_DED; // histogram of the number of jets using DED jets

JetRazor JetRazor;

JetRazor::JetRazor() : Processor("JetRazor") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") ); 
    registerProcessorParameter( "jetDetectability",
            "Detectability of the Thrust Axis/Value to be used:\n#\t0 : True \n#t1 : Detectable \n#t2 : Detected" ,
            _jetDetectability, 1 );
}


void JetRazor::init() { 
    streamlog_out(DEBUG)  << "   init called  " << std::endl ;

    if(_jetDetectability==0){_rootfile = new TFile("JetRazor_.39133._T8.0.root","RECREATE");
        _R_T = new TH1F("R_T", "R =MTR/MR",130,-3,10);
        _MR_T = new TH1F("MR_T","MR", 100, 0 ,10);
        _NJ_T = new TH1F("NJ_T","NJ",20, 0, 20); 
    }
    if(_jetDetectability==1){_rootfile = new TFile("JetRazor_.39133._DAB3.0.root","RECREATE");
        _MR_DAB = new TH1F("MR_DAB","MR", 100, 0 ,10); 
        _R_DAB = new TH1F("R_DAB", "R =MTR/MR",130,-3,10);
        _NJ_DAB = new TH1F("NJ_DAB","NJ",20,0,20);
    }
    if(_jetDetectability==2){_rootfile = new TFile("JetRazor_.39133._DED.root","RECREATE");
        _MR_DED = new TH1F("MR_DED","MR", 100, 0 ,10); 
        _R_DED = new TH1F("R_DED", "R =MTR/MR",130,-3,10);
    }
    if(_jetDetectability==0){ freopen( "JetRazor_.39133._T8.0.log", "w", stdout ); }
    if(_jetDetectability==1){ freopen( "JetRazor_.39133._DAB3.0.log", "w", stdout ); }
    if(_jetDetectability==2){ freopen( "JetRazor_.39133._DED8.0.log", "w", stdout ); }
    // irameters() ;

    // config ranlux 
    filename = "Ranlux.coonf";
    ifstream rndcfgfile( filename.c_str() );

    if (!rndcfgfile)
    {
        long int ss=1234;
        myrnd.setSeeds(&ss,4);
        myrnd.showStatus();
    }
    else
    {
        rndcfgfile.close();
        myrnd.restoreStatus(filename.c_str());
        myrnd.showStatus();
    } // if file not existusually a good idea to
    //printParameters() ;
    _nEvt = 0 ; 

    pname[11  ]="              electron";
    pname[12  ]="     electron neutrino";
    pname[-12 ]="anti electron neutrino";
    pname[13  ]="                  muon";
    pname[-13 ]="             anti muon";
    pname[14  ]="         muon neutrino";
    pname[-14 ]="    anti muon neutrino";
    pname[15  ]="                   tau";
    pname[16  ]="          tau neutrino";
    pname[-16 ]="     anti tau neutrino";
    pname[22  ]="                photon";
    pname[1000022]="              X10";
    pname[211 ]="                   pi+";
    pname[-211]="                   pi-";

}

void JetRazor::processRunHeader( LCRunHeader* run) { 
    //run->parameters().setValue("thrust",12300321);
    //    _nRun++ ;
} 

void JetRazor::processEvent( LCEvent * evt ) { 


    parp1 = false;
    parp0 = false;  
    // usually the working horse ...
    cout << "EVENT: " << _nEvt << endl; 
    _inParVec = evt->getCollection( _colName) ;
    //cout << "num of elements " << _inParVec->getNumberOfElements() << endl;
    if (!_parp.empty()) _parp.clear();
    int id, stat;

    for (int n=0;n<_inParVec->getNumberOfElements() ;n++)
    {
        MCParticle* aPart = dynamic_cast<MCParticle*>( _inParVec->getElementAt(n) );
        try{
            id = aPart->getPDG();
            stat = aPart->getGeneratorStatus();
        }
        catch(const std::exception& e){
            cout << "exception caught with message " << e.what() << "\n";
        }
        const double* parp = aPart->getMomentum();
        double parp_mag = sqrt(parp[0]*parp[0]+parp[1]*parp[1]+parp[2]*parp[2]);

        if(stat==1){
            cout << "id: " << id << " " << pname[id] <<" " << parp[0]<<" "<< parp[1]<<" "<< parp[2] <<endl; 
            double penergy = aPart->getEnergy(); 
            bool isDarkMatter = (id == 1000022);
            bool isNeutrino = (
                    id == 12 || id == -12 ||
                    id == 14 || id == -14 ||
                    id == 16 || id == -16 ||
                    id == 18 || id == -18);

            double cos = parp[2]/parp_mag;
            bool isForward = ( cos > 0.9 || cos < - 0.9);
            bool isDetectable = (!isDarkMatter && !isNeutrino);
            bool isDetected = (isDetectable &&  !isForward  );
            if(_jetDetectability == 0){
                if(!isDarkMatter){ 
                    _parp.push_back( PseudoJet(parp[0], parp[1], parp[2], penergy) );
                }
            }
            if(_jetDetectability == 1){
                if(isDetectable){ 
                    _parp.push_back( PseudoJet(parp[0], parp[1], parp[2], penergy) );
                }
            }

            if(_jetDetectability == 2){ 
                if(isDetected){  
                    _parp.push_back( PseudoJet(parp[0], parp[1], parp[2], penergy) ); 
                }
            }
        } // stat = 1
    } // for particle  

    JetDefinition jet_def(antikt_algorithm,   _JetRParameter   ); 
    // run the clustering, extract the jets
    ClusterSequence cs(_parp, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
    // print out some info
    //cout << "Clustered with " << jet_def.description() << endl;
    // print the jets
    //cout << " pt y phi" << endl;
    cout << "NUMBER OF JETS: "<< jets.size() << endl; 
    for (unsigned i = 0; i < jets.size(); i++) {
        cout << "jet " << i << ": "<< jets[i].perp() << " "<< jets[i].rap() << " energy: " <<jets[i].e()<< endl; 
        vector<PseudoJet> constituents = jets[i].constituents();
        for (unsigned j = 0; j < constituents.size(); j++) {
            //cout << " constituent " << j << "’s pt: "<< constituents[j].perp() << endl; 
            cout <<" constituent "<< j << " "<< constituents[j].px() << " "<<constituents[j].py() << " "<< constituents[j].pz()<< endl;
        }
    }

    vector<PseudoJet> megajets = findMegaJets(jets);
    //_NJ_T->Fill(jets.size());
    _NJ_DAB->Fill(jets.size());
    //_NJ_DED->Fill(jets.size());
    //_R_DAB->Fill(R);
    //_R_DED->Fill(R);

    // megajet algorithm: the combination of jets into two disjoint sets (megajets) is the combination that minimized the 
    // sum of the squared-invariant-mass values of the two megajets 

    double vec[2][3][4]; // jet 1, jet 2 : true detectable, detected : energy, momx, momy, momz
    double Rvec[3][4]; // true, detectable, detected : energy, px, py, pz 
    int d = _jetDetectability;
    double R;
    double beta2;

    //int id, stat;

    for (int n=0;n<_inParVec->getNumberOfElements() ;n++){

        MCParticle* aPart = dynamic_cast<MCParticle*>( _inParVec->getElementAt(n) );

        try{
            id = aPart->getPDG();
            stat = aPart->getGeneratorStatus();
        }
        catch(const std::exception& e){
            cout << "exception caught with message " << e.what() << "\n";
        }

        if(stat==1){
            const double* parp = aPart->getMomentum();
            double part4mom[4];

            part4mom[0] = aPart->getEnergy(); 
            part4mom[1] = parp[0];
            part4mom[2] = parp[1]; 
            part4mom[3] = parp[2]; 


            // need momentum and energy of entire jet 

            bool isDarkMatter = (id == 1000022);
            bool isNeutrino = (
                    id == 12 || id == -12 ||
                    id == 14 || id == -14 ||
                    id == 16 || id == -16 ||
                    id == 18 || id == -18);
            double cos = parp[2]/(sqrt(parp[0]*parp[0]+parp[1]*parp[1]+parp[2]*parp[2]));
            bool isForward = ( cos > 0.9 || cos < - 0.9);
            int i = 1; // jet #

            //if(dot>=0){i=0;}
            //if(dot<0){i=1;}
            if (!isDarkMatter){

                vec[i][0][0]+= part4mom[0]; 
                vec[i][0][1]+= part4mom[1];
                vec[i][0][2]+= part4mom[2];
                vec[i][0][3]+= part4mom[3];
                if (!isNeutrino){
                    vec[i][1][0]+= part4mom[0];
                    vec[i][1][1]+= part4mom[1];
                    vec[i][1][2]+= part4mom[2];
                    vec[i][1][3]+= part4mom[3];
                    if(!isForward){
                        vec[i][2][0]+= part4mom[0];
                        vec[i][2][1]+= part4mom[1];
                        vec[i][2][2]+= part4mom[2];
                        vec[i][2][3]+= part4mom[3];
                    }
                }
            }       
        }
    }

    //int d = _thrustDetectability;
    double beta = (vec[0][d][0]-vec[1][d][0])/(vec[0][d][3]-vec[1][d][3]); // beta using right particles 
    //double beta2 = pow(beta,2);
    beta2 = pow(beta,2);
    double gamma = 1/(sqrt(1-beta2));
    for (int n=0;n<_inParVec->getNumberOfElements() ;n++){

        MCParticle* aPart = dynamic_cast<MCParticle*>( _inParVec->getElementAt(n) );
        const double* parp = aPart->getMomentum();

        try{
            id = aPart->getPDG();
            stat = aPart->getGeneratorStatus();
        }
        catch(const std::exception& e){
            cout << "exception caught with message " << e.what() << "\n";
        }
        double part4Vec[4] = {aPart->getEnergy(), parp[0], parp[1], parp[2] };
        double R4Vec[4] = {gamma*part4Vec[0]-gamma*beta*part4Vec[3], part4Vec[1], part4Vec[2], 
            -gamma*beta*part4Vec[0]+gamma*part4Vec[3] }; 
        bool isDarkMatter = (id == 1000022);

        bool isNeutrino = (
                id == 12 || id == -12 ||
                id == 14 || id == -14 ||
                id == 16 || id == -16 ||
                id == 18 || id == -18);
        double cos = parp[2]/(sqrt(parp[0]*parp[0]+parp[1]*parp[1]+parp[2]*parp[2]));
        bool isForward = (cos > 0.9 || cos < - 0.9);
        bool isDetectable = (!isDarkMatter && !isNeutrino);
        bool isDetected = (isDetectable && !isForward); 
        if(stat ==1){
            //cout << "id: "<<id<<endl;
            //cout << parp<< endl;
            if(_jetDetectability == 0){
                if(!isDarkMatter){
                    Rvec[d][0]+=R4Vec[0];
                    Rvec[d][1]+=R4Vec[1];
                    Rvec[d][2]+=R4Vec[2];
                    Rvec[d][3]+=R4Vec[3];
                }
            }
            if(_jetDetectability ==1){
                if(isDetectable){
                    Rvec[d][0]+=R4Vec[0];
                    Rvec[d][1]+=R4Vec[1];
                    Rvec[d][2]+=R4Vec[2];
                    Rvec[d][3]+=R4Vec[3];
                }
            }
            if(_jetDetectability == 2){
                if(isDetected){
                    Rvec[d][0]+=R4Vec[0];
                    Rvec[d][1]+=R4Vec[1];
                    Rvec[d][2]+=R4Vec[2];
                    Rvec[d][3]+=R4Vec[3];
                }
            }   
        }
    }
    double ETM[2] = {-Rvec[d][1], - Rvec[d][2]};
    double ETMmag = sqrt(ETM[0]*ETM[0]+ETM[1]*ETM[1]);

    double ptj1mag = sqrt(vec[0][d][1]*vec[0][d][1]+vec[0][d][2]*vec[0][d][2]);
    double ptj2mag = sqrt(vec[1][d][1]*vec[1][d][1]+vec[1][d][2]*vec[1][d][2]);
    double ptmagsum = ptj1mag + ptj2mag; 

    double ptvecsum[2] = {vec[0][d][1]+vec[1][d][1], vec[0][d][2]+vec[1][d][2]};
    double ETMdotptvecsum = ETM[0]*ptvecsum[0]+ETM[1]*ptvecsum[1];
    double MTR = sqrt((ETMmag*ptmagsum-ETMdotptvecsum)/2);

    double pj1 = sqrt(vec[0][d][1]*vec[0][d][1]+vec[0][d][2]*vec[0][d][2]+vec[0][d][3]*vec[0][d][3]);
    R = MTR/(2*pj1);

    if(beta2<=1){
        if(d==0){_R_T->Fill(R);}
        if(d==1){_R_DAB->Fill(R);}
        if(d==2){_R_DED->Fill(R);}
        if(parp1){
            // cout << "there are valid beta events that have only 1 particle" << endl; 
        }
        if(!parp1){
            //   cout << "there are valid beta events that have not only 1 particle" <<endl; 
        }
    }
    else{
        if(d==0){_R_T->Fill(-2);}
        if(d==1){_R_DAB->Fill(-2);}
        if(d==2){_R_DED->Fill(-2);}

        if(parp1){
            //cout << "there are inval beta events that have only 1 particle"<< endl;
        }
        if(!parp1){
            //cout << "there are inval beta events that have not only 1 particle"<< endl;
            if(parp0){
                //cout <<" zero particles" << endl;
            }
            if(!parp0){
                //cout << "not zero particles"<< endl;
            }
        }

        betaCheck += " ";
        betaCheck += std::to_string(_nEvt);
        betaCheck += " ";  
    }  
    cout << "End EVENT "<< _nEvt<< endl;

    _nEvt ++ ; // different from original-moved out of for loop - summer 
}




void JetRazor::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void JetRazor::end(){ 
    _rootfile->Write();
    cout << parpCheck << endl;
    cout << betaCheck << endl; 
}

