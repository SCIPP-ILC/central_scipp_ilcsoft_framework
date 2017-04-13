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
 *
 * author Summer Zuber 
 * March 4, 2017 
 *
 * In order to analyze the kinematic event oversables 
 * 'S', 'V', and 'M' of the hadronic system of two-photon events 
 * (The hadronic system being all final state particles less the initial 
 * scattered electron and positron)
 * S = sum of the magnitudes of the transerve momentum
 * of each particle 
 * V = magnitude of the vector sum of teh transverse mometum 
 * of each particle 
 * M = sqrt((sum(energies))^(2)-((sum(px))^(2)+(sum(py))^(2)+(sum(pz))^(2)) )
 * (mass of event) 
 *
 * We analyze these at three different levels of "detectability"
 * True: entire hadronic system
 * Detectable: True less neutrinos 
 * Detected: Detectable less forward particle with |cos(theta)| > 0.9 
 *
 */

#include "TwoPhotonSVM.h"
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

// #include <CLHEP/Vector/ThreeVector.h>
// #include <CLHEP/Random/RanluxEngine.h>

#include <EVENT/LCCollection.h>
#include <EVENT/LCIO.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/LCTOOLS.h>


using namespace lcio;
using namespace marlin;
using namespace std;


static TFile* _rootfile;

static TH1F* _S_TRU;
static TH1F* _S_DAB;
static TH1F* _S_DED;
static TH1F* _V_TRU;
static TH1F* _V_DAB;
static TH1F* _V_DED;
static TH1F* _M_TRU;
static TH1F* _M_DAB;
static TH1F* _M_DED;

TwoPhotonSVM TwoPhotonSVM;

TwoPhotonSVM::TwoPhotonSVM() : Processor("TwoPhotonSVM") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}

void TwoPhotonSVM::init() { 
    streamlog_out(DEBUG)  << "   init called  " << std::endl ;
    
    _rootfile = new TFile("TwoPhotonSVM_eW.pW.I39212.root", "RECREATE");

    _S_TRU = new TH1F("S_TRU", "True Scalar",80,0,20);
    _S_DAB = new TH1F("S_DAB", "Detectable Scalar",80,0,20); 
    _S_DED = new TH1F("S_DED", "Detected Scalar",80,0,20); 
   
    _V_TRU = new TH1F("V_TRU", "True Vector",80,0,20);
    _V_DAB = new TH1F("V_DAB", "Detectable Vector",80,0,20); 
    _V_DED = new TH1F("V_DED", "Detected Vector",80,0,20); 
    
    _M_TRU = new TH1F("M_TRU", "True Mass",80,0,20);
    _M_DAB = new TH1F("M_DAB", "Detectable Mass",80,0,20); 
    _M_DED = new TH1F("M_DED", "Detected Mass",80,0,20); 
    
    //irameters() ;
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
    }
    // if file not existusually a good idea to
    //printParameters() ;
    _nEvt = 0 ;
}



void TwoPhotonSVM::processRunHeader( LCRunHeader* run) { 
    //run->parameters().setValue("thrust",12300321);
    //    _nRun++ ; 

} 

void TwoPhotonSVM::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    _inParVec = evt->getCollection( _colName) ; 

    int id, stat;
    
    double SVM[3][5]; // TRU, DAB, DED : E, px, py, pz, scalar
    int electronI = -1;
    int positronI = -1;
    double electronEnergy = 0;
    double positronEnergy = 0;
   
    for (int n=0;n<_inParVec->getNumberOfElements(); n++){
        MCParticle* aPart = dynamic_cast<MCParticle*>(_inParVec->getElementAt(n) ); 
        try {
            id = aPart->getPDG();
            stat = aPart->getGeneratorStatus();
        }
        catch(const std::exception& e){
            cout << "exception caught with message "<< e.what() <<"\n";
        }
        if (stat != 1) continue;
        double energy = aPart->getEnergy();
        if(id == 11 && energy > electronEnergy){
            electronI = n;
            electronEnergy = energy;
        }
        if(id== -11 && energy > positronEnergy){
            positronI = n;
            positronEnergy = energy;
        }
    }
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
        if(stat==1){
            if(n==electronI) continue ; // not including the 
            if(n==positronI) continue ; // initial electron and positron 
            cout << "id: "<< id<< endl;
            const double* P = aPart->getMomentum();
            double PMag = sqrt(P[0]*P[0]+P[1]*P[1]+P[2]*P[2]);
            double E = aPart->getEnergy();
            double scalar = sqrt(P[0]*P[0]+P[1]*P[1]);
            double cos = P[2]/PMag;
            bool isForward = (cos > 0.9 || cos < -0.9); 
            bool isNeutrino = (
                    id == 12 || id == -12 ||
                    id == 14 || id == -14 ||
                    id == 16 || id == -16 ||
                    id == 18 || id == -18);
            bool isDetectable = (!isNeutrino);
            bool isDetected = (isDetectable && !isForward); 
            SVM[0][0]+=E;
            SVM[0][1]+=P[0];
            SVM[0][2]+=P[1];
            SVM[0][3]+=P[2];
            SVM[0][4]+=scalar;
            if(isDetectable){
                SVM[1][0]+=E;
                SVM[1][1]+=P[0];
                SVM[1][2]+=P[1];
                SVM[1][3]+=P[2];
                SVM[1][4]+=scalar;
            }
            if(isDetected){
                SVM[2][0]+=E;
                SVM[2][1]+=P[0];
                SVM[2][2]+=P[1];
                SVM[2][3]+=P[2];
                SVM[2][4]+=scalar;
            }
        } //stat==1
    } // for particle
    
    // for all particles in event 
    double total_TRU_S = SVM[0][4];
    double total_DAB_S = SVM[1][4];
    double total_DED_S = SVM[2][4];

    double total_TRU_V_sq = SVM[0][1]*SVM[0][1]+SVM[0][2]*SVM[0][2]; 
    double total_DAB_V_sq = SVM[1][1]*SVM[1][1]+SVM[1][2]*SVM[1][2];
    double total_DED_V_sq = SVM[2][1]*SVM[2][1]+SVM[2][2]*SVM[2][2];

    double total_TRU_V = sqrt(total_TRU_V_sq);
    double total_DAB_V = sqrt(total_DAB_V_sq);
    double total_DED_V = sqrt(total_DED_V_sq);

    double total_TRU_M_sq = SVM[0][0]*SVM[0][0] - (SVM[0][1]*SVM[0][1]+SVM[0][2]*SVM[0][2]+SVM[0][3]*SVM[0][3]);
    double total_DAB_M_sq = SVM[1][0]*SVM[1][0] - (SVM[1][1]*SVM[1][1]+SVM[1][2]*SVM[1][2]+SVM[1][3]*SVM[1][3]);
    double total_DED_M_sq = SVM[2][0]*SVM[2][0] - (SVM[2][1]*SVM[2][1]+SVM[2][2]*SVM[2][2]+SVM[2][3]*SVM[2][3]);

    double total_TRU_M = sqrt(total_TRU_M_sq);
    double total_DAB_M = sqrt(total_DAB_M_sq);
    double total_DED_M = sqrt(total_DED_M_sq);

    _S_TRU->Fill(total_TRU_S);
    _S_DAB->Fill(total_DAB_S);
    _S_DED->Fill(total_DED_S);
    
    _V_TRU->Fill(total_TRU_V);
    _V_DAB->Fill(total_DAB_V);
    _V_DED->Fill(total_DED_V);
    
    _M_TRU->Fill(total_TRU_M);
    _M_DAB->Fill(total_DAB_M);
    _M_DED->Fill(total_DED_M);
    
    cout << "end Event "<<_nEvt <<endl; 
    _nEvt ++ ; // different from original-moved out of for loop - summer                     
}  // end process Event 




void TwoPhotonSVM::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void TwoPhotonSVM::end(){ 
    _rootfile->Write();
}

