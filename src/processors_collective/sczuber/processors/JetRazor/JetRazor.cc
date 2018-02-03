#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0
/* 
 * author Summer Zuber 
 * January 18, 2017 
 *
 * Using FastJet
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

// #include <CLHEP/Vector/ThreeVector.h>
// #include <CLHEP/Random/RanluxEngine.h>

#include <EVENT/LCCollection.h>
#include <EVENT/LCIO.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/LCTOOLS.h>

#include "fastjet/ClusteringSequence.hh"

using namespace lcio;
using namespace marlin;
using namespace std;
using namespace fastjet;


static TFile* _rootfile;

static TH1F* _R_T;
static TH1F* _R_DAB;
static TH1F* _R_DED;
static TH1F* _MR_T;
static TH1F* _MR_DAB;
static TH1F* _MR_DED; 

JetRazor JetRazor;

JetRazor::JetRazor() : Processor("JetRazor") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
    registerProcessorParameter( "typeOfJetRazorFinder" ,
            "Type of thrust reconstruction algorithm to be used:\n#\t1 : Tasso algorithm\n#\t2 : JetSet algorithm"  ,
            _typeOfJetRazorFinder , 2 ) ;
    registerProcessorParameter( "thrustDetectability",
            "Detectability of the Thrust Axis/Value to be used:\n#\t0 : True \n#t1 : Detectable \n#t2 : Detected" ,
            _thrustDetectability, 2 );
}



void JetRazor::init() { 
    streamlog_out(DEBUG)  << "   init called  " << std::endl ;
    cout << "initialized" << endl;
    if(_thrustDetectability==0){_rootfile = new TFile("JetRazor_.39133._T.root","RECREATE");
        _R_T = new TH1F("R_T", "R =MTR/MR",130,-3,10);
        _MR_T = new TH1F("MR_T","MR", 100, 0 ,10); 
    }
    if(_thrustDetectability==1){_rootfile = new TFile("JetRazor_.39133._DAB.root","RECREATE");
        _MR_DAB = new TH1F("MR_DAB","MR", 100, 0 ,10); 
        _R_DAB = new TH1F("R_DAB", "R =MTR/MR",130,-3,10);
    }
    if(_thrustDetectability==2){_rootfile = new TFile("JetRazor_.39133._DED.root","RECREATE");
        _MR_DED = new TH1F("MR_DED","MR", 100, 0 ,10); 
        _R_DED = new TH1F("R_DED", "R =MTR/MR",130,-3,10);
    }
    
    freopen( "JetRazor.log", "w", stdout );
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
}

void JetRazor::processRunHeader( LCRunHeader* run) { 
    //run->parameters().setValue("thrust",12300321);
    //    _nRun++ ;
} 

void JetRazor::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    partMom1 = false;
    partMom0 = false;  
    // usually the working horse ...
    cout << "EVENT: " << _nEvt << endl; 
    _inParVec = evt->getCollection( _colName) ;
    cout << "num of elements " << _inParVec->getNumberOfElements() << endl;
    if (!_partMom.empty()) _partMom.clear();

    int id, stat;
    cout << "loop #1"<< endl; 
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

        const double* partMom = aPart->getMomentum();
        double partMomMag = sqrt(partMom[0]*partMom[0]+partMom[1]*partMom[1]+partMom[2]*partMom[2]);

        if(stat==1){
            cout << "id: " << id<< endl;
            cout << "mom: "<< partMom[0]<<" "<< partMom[1]<<" "<<partMom[2]<<endl;
            bool isDarkMatter = (id == 1000022);
            bool isNeutrino = (
                    id == 12 || id == -12 ||
                    id == 14 || id == -14 ||
                    id == 16 || id == -16 ||
                    id == 18 || id == -18);

            double cos = partMom[2]/partMomMag;
            bool isForward = ( cos > 0.9 || cos < - 0.9);
            bool isDetectable = (!isDarkMatter && !isNeutrino);
            bool isDetected = (isDetectable &&  !isForward  );
            if(_thrustDetectability == 0){
                if(!isDarkMatter){
                    _partMom.push_back( Hep3Vector(partMom[0], partMom[1], partMom[2]) );
                }
            }
            if(_thrustDetectability == 1){
                if(isDetectable){
                    _partMom.push_back( Hep3Vector(partMom[0], partMom[1], partMom[2]) );
                }
            }

            if(_thrustDetectability == 2){ 
                if(isDetected){ 
                    _partMom.push_back( Hep3Vector(partMom[0], partMom[1], partMom[2]) ); 
                }
            }
        } // stat = 1
    } // for particle 
    cout << "end loop #1"<<endl;

    //reset variables for output   
    _principleJetRazorValue = -1;
    _majorJetRazorValue     = -1;
    _minorJetRazorValue     = -1;
    _principleJetRazorAxis.set(0,0,0);
    _majorJetRazorAxis.set(0,0,0);
    _minorJetRazorAxis.set(0,0,0);

    // Switch to the desired type of thrust finder
    if (_typeOfJetRazorFinder == 1)
    {
        cout << "type of Thrust Razor Finder = 1 : Tasso Thrust Razor " << endl; 
        TassoJetRazor();
        cout << "type of Thrust Razor Finder = 1 : Tasso Thrust Razor" << endl;
    }
    else if (_partMom.size()<=1)
    {
        partMomCheck += " ";
        partMomCheck += std::to_string(_nEvt);
        partMomCheck += " ";  
        TassoJetRazor();
    }
    else if (_typeOfJetRazorFinder == 2)
    {
        cout << "type of Thrust Razor Finder = 2 : Jetset Thrust Razor" << endl; 
        JetsetJetRazor();
    }
    // ###write
    //    evt->parameters().setValue("thrust",_principleJetRazorValue);

    FloatVec thrax;
    thrax.clear();
    thrax.push_back(_principleJetRazorAxis.x());
    thrax.push_back(_principleJetRazorAxis.y());
    thrax.push_back(_principleJetRazorAxis.z());

    _inParVec->parameters().setValue("principleJetRazorValue",_principleJetRazorValue);
    _inParVec->parameters().setValues("principleJetRazorAxis",thrax);

    if (_typeOfJetRazorFinder == 2)
    {
        thrax.clear();
        thrax.push_back(_majorJetRazorAxis.x());
        thrax.push_back(_majorJetRazorAxis.y());
        thrax.push_back(_majorJetRazorAxis.z());

        _inParVec->parameters().setValue("majorJetRazorValue",_majorJetRazorValue);
        _inParVec->parameters().setValues("majorJetRazorAxis",thrax);

        thrax.clear();
        thrax.push_back(_minorJetRazorAxis.x());
        thrax.push_back(_minorJetRazorAxis.y());
        thrax.push_back(_minorJetRazorAxis.z());

        _inParVec->parameters().setValue("minorJetRazorValue",_minorJetRazorValue);
        _inParVec->parameters().setValues("minorJetRazorAxis",thrax);

        float Oblateness;
        Oblateness = _majorJetRazorValue - _minorJetRazorValue;
        _inParVec->parameters().setValue("Oblateness",Oblateness);
        if ( (_majorJetRazorValue < 0) || (_minorJetRazorValue < 0) )
        {
            _inParVec->parameters().setValue("Oblateness",-1);
        }
    }


    // these are the final values:
    streamlog_out( DEBUG4 ) << " thrust: " << _principleJetRazorValue << " TV: " << _principleJetRazorAxis << endl;
    streamlog_out( DEBUG4 ) << "  major: " << _majorJetRazorValue << " TV: " << _majorJetRazorAxis << endl;
    streamlog_out( DEBUG4 ) << "  minor: " << _minorJetRazorValue << " TV: " << _minorJetRazorAxis << endl;
    //cout << "EVENT: " << _nEvt << endl;
    cout << " thrust: " << _principleJetRazorValue << " TV: " << _principleJetRazorAxis << endl;
    cout <<"                       "<< _principleJetRazorAxis.x()<<","<< _principleJetRazorAxis.y()<< ","<<_principleJetRazorAxis.z()<<endl;
    cout << "  major: " << _majorJetRazorValue << " TV: " << _majorJetRazorAxis << endl;
    cout << "  minor: " << _minorJetRazorValue << " TV: " << _minorJetRazorAxis << endl;

    double ptaX = _principleJetRazorAxis.x();
    double ptaY = _principleJetRazorAxis.y();
    double ptaZ = _principleJetRazorAxis.z();

    double vec[2][3][4]; // jet 1, jet 2 : true detectable, detected : energy, momx, momy, momz
    double Rvec[3][4]; // true, detectable, detected : energy, px, py, pz 
    int d = _thrustDetectability;
    double R;
    double beta2;
    if(_principleJetRazorValue==-1){
        R = -1; 
    }
    else{
        //int id, stat;
        cout << "start loop 2" << endl;  
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
                const double* partMom = aPart->getMomentum();
                double part4mom[4];

                part4mom[0] = aPart->getEnergy(); 
                part4mom[1] = partMom[0];
                part4mom[2] = partMom[1]; 
                part4mom[3] = partMom[2];   
                double pta[3] = {ptaX, ptaY, ptaZ};

                cout << "id      : " << id << endl;  
                cout << "Momentum: " << partMom[0] <<" "<< partMom[1] <<" "<< partMom[2]<< endl;
                cout << "Thrust A: " << ptaX << " "<< ptaY << " " << ptaZ << endl;
                double dot = ptaX*partMom[0]+ptaY*partMom[1]+ptaZ*partMom[2];
                cout << "dot " << dot << endl;

                // need momentum and energy of entire jet 

                bool isDarkMatter = (id == 1000022);
                bool isNeutrino = (
                        id == 12 || id == -12 ||
                        id == 14 || id == -14 ||
                        id == 16 || id == -16 ||
                        id == 18 || id == -18);
                double cos = partMom[2]/(sqrt(partMom[0]*partMom[0]+partMom[1]*partMom[1]+partMom[2]*partMom[2]));
                bool isForward = ( cos > 0.9 || cos < - 0.9);
                int i; // jet #
                cout << "dot: " <<dot << endl; 
                if(dot>=0){i=0;}
                if(dot<0){i=1;}
                if (!isDarkMatter){
                    cout << "i: "<< i << endl;  
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
        cout << "end loop 3 "<< endl; 
        //int d = _thrustDetectability;
        double beta = (vec[0][d][0]-vec[1][d][0])/(vec[0][d][3]-vec[1][d][3]); // beta using right particles 
        //double beta2 = pow(beta,2);
        beta2 = pow(beta,2);
        double gamma = 1/(sqrt(1-beta2));
        for (int n=0;n<_inParVec->getNumberOfElements() ;n++){

            MCParticle* aPart = dynamic_cast<MCParticle*>( _inParVec->getElementAt(n) );
            const double* partMom = aPart->getMomentum();

            try{
                id = aPart->getPDG();
                stat = aPart->getGeneratorStatus();
            }
            catch(const std::exception& e){
                cout << "exception caught with message " << e.what() << "\n";
            }
            double part4Vec[4] = {aPart->getEnergy(), partMom[0], partMom[1], partMom[2] };
            double R4Vec[4] = {gamma*part4Vec[0]-gamma*beta*part4Vec[3], part4Vec[1], part4Vec[2], 
                -gamma*beta*part4Vec[0]+gamma*part4Vec[3] }; 
            bool isDarkMatter = (id == 1000022);

            bool isNeutrino = (
                    id == 12 || id == -12 ||
                    id == 14 || id == -14 ||
                    id == 16 || id == -16 ||
                    id == 18 || id == -18);
            double cos = partMom[2]/(sqrt(partMom[0]*partMom[0]+partMom[1]*partMom[1]+partMom[2]*partMom[2]));
            bool isForward = (cos > 0.9 || cos < - 0.9);
            bool isDetectable = (!isDarkMatter && !isNeutrino);
            bool isDetected = (isDetectable && !isForward); 
            if(stat ==1){
            cout << "id: "<<id<<endl;
            cout << partMom<< endl;
                if(_thrustDetectability == 0){
                    if(!isDarkMatter){
                        Rvec[d][0]+=R4Vec[0];
                        Rvec[d][1]+=R4Vec[1];
                        Rvec[d][2]+=R4Vec[2];
                        Rvec[d][3]+=R4Vec[3];
                    }
                }
                if(_thrustDetectability ==1){
                    if(isDetectable){
                        Rvec[d][0]+=R4Vec[0];
                        Rvec[d][1]+=R4Vec[1];
                        Rvec[d][2]+=R4Vec[2];
                        Rvec[d][3]+=R4Vec[3];
                    }
                }
                if(_thrustDetectability == 2){
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
    }
    if(beta2<=1){
        if(d==0){_R_T->Fill(R);}
        if(d==1){_R_DAB->Fill(R);}
        if(d==2){_R_DED->Fill(R);}
        if(partMom1){
            cout << "there are valid beta events that have only 1 particle" << endl; 
        }
        if(!partMom1){
            cout << "there are valid beta events that have not only 1 particle" <<endl; 
        }
    }
    else{
        if(d==0){_R_T->Fill(-2);}
        if(d==1){_R_DAB->Fill(-2);}
        if(d==2){_R_DED->Fill(-2);}
        
        if(partMom1){
            cout << "there are inval beta events that have only 1 particle"<< endl;
        }
        if(!partMom1){
            cout << "there are inval beta events that have not only 1 particle"<< endl;
            if(partMom0){
                cout <<" zero particles" << endl;
            }
            if(!partMom0){
                cout << "not zero particles"<< endl;
            }
        }
        
        betaCheck += " ";
        betaCheck += std::to_string(_nEvt);
        betaCheck += " ";  
    } 

    if (_principleJetRazorValue >= _max) _max = _principleJetRazorValue;
    if (_principleJetRazorValue <= _min) _min = _principleJetRazorValue;
    cout << "End EVENT "<< _nEvt<< endl;

    _nEvt ++ ; // different from original-moved out of for loop - summer 
}




void JetRazor::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void JetRazor::end(){ 
    _rootfile->Write();
    cout << partMomCheck << endl;
    cout << betaCheck << endl; 
}

