#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0
/* 
 * Razor.cc
 * author Summer Zuber 
 * January 24, 2017 
 *
 * Copied from ThurstRazor.cc. Will use the JetFinder utility (translated to C++/ROOT by M. Iwasaki) to identify jets for use in the Razor variables. 
 */

#include <stdlib.h>
#include "Razor.h"

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

//#include <TLorentzVector.h>
#include <TObjArray.h>
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

static TH1F* _R_T;
static TH1F* _R_DAB;
static TH1F* _R_DED;
static TH1F* _MR_T;
static TH1F* _MR_DAB;
static TH1F* _MR_DED; 

Razor Razor;
//JetFinder JetFinder; 

Razor::Razor() : Processor("Razor") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
    registerProcessorParameter( "typeOfRazorFinder" ,
            "Type of thrust reconstruction algorithm to be used:\n#\t1 : Tasso algorithm\n#\t2 : JetSet algorithm"  ,
            _typeOfRazorFinder , 2 ) ;
    registerProcessorParameter( "thrustDetectability",
            "Detectability of the Thrust Axis/Value to be used:\n#\t0 : True \n#t1 : Detectable \n#t2 : Detected" ,
            _thrustDetectability, 2 );
}



void Razor::init() { 
    streamlog_out(DEBUG)  << "   init called  " << std::endl ;
    cout << "initialized" << endl;
    if(_thrustDetectability==0){_rootfile = new TFile("Razor_.39133._T.root","RECREATE");
        _R_T = new TH1F("R_T", "R =MTR/MR",130,-3,10);
        _MR_T = new TH1F("MR_T","MR", 100, 0 ,10); 
    }
    if(_thrustDetectability==1){_rootfile = new TFile("Razor_.39133._DAB.root","RECREATE");
        _MR_DAB = new TH1F("MR_DAB","MR", 100, 0 ,10); 
        _R_DAB = new TH1F("R_DAB", "R =MTR/MR",130,-3,10);
    }
    if(_thrustDetectability==2){_rootfile = new TFile("Razor_.39133._DED.root","RECREATE");
        _MR_DED = new TH1F("MR_DED","MR", 100, 0 ,10); 
        _R_DED = new TH1F("R_DED", "R =MTR/MR",130,-3,10);
    }
    
    freopen( "Razor.log", "w", stdout );
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

void Razor::processRunHeader( LCRunHeader* run) { 
    //run->parameters().setValue("thrust",12300321);
    //    _nRun++ ;
} 

void Razor::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    parp1 = false;
    parp0 = false;  
    // usually the working horse ...
    cout << "EVENT: " << _nEvt << endl; 
    _inParVec = evt->getCollection( _colName) ;
    cout << "num of elements " << _inParVec->getNumberOfElements() << endl;
    //if (!_parp.empty()) _parp.clear(); NEED TO REPLACE! 

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
        double parpMag = sqrt(parp[0]*parp[0]+parp[1]*parp[1]+parp[2]*parp[2]);
        double par4p[4];
        par4p[0] += aPart->getEnergy();
        _parp->SetOwner(kTRUE); 
        if(stat==1){ 
           
            bool isDarkMatter = (id == 1000022);
            bool isNeutrino = (
                    id == 12 || id == -12 ||
                    id == 14 || id == -14 ||
                    id == 16 || id == -16 ||
                    id == 18 || id == -18);

            double cos = parp[2]/parpMag;
            bool isForward = ( cos > 0.9 || cos < - 0.9);
            bool isDetectable = (!isDarkMatter && !isNeutrino);
            bool isDetected = (isDetectable &&  !isForward  );
            if(_thrustDetectability == 0){
                if(!isDarkMatter){
                    TLorentzVector* v = new TLorentzVector(parp[0], parp[1], parp[2], par4p[0]);
                    _parp->Add(v);
                    //cout << _parp->getName();
                    //_parp.Add(TLorentzVector(parp[0], parp[1], parp[2], par4p[0])) ;
                }
            }
            /*if(_thrustDetectability == 1){
                if(isDetectable){
                    _parp.Add( TLorentzVector(parp[0], parp[1], parp[2], par4p[0]) );
                }
            }

            if(_thrustDetectability == 2){ 
                if(isDetected){ 
                    _parp.Add( TLorentzVector(parp[0], parp[1], parp[2], par4p[0]) ); 
                }
            }*/
        } // stat = 1
    } // for particle 
  

    // JetFinder stuff -----------------
    //jetFinder.setPartList(_parp);
    
    //reset variables for output   
    _principleRazorValue = -1;
    _majorRazorValue     = -1;
    _minorRazorValue     = -1;
    _principleRazorAxis.set(0,0,0);
    _majorRazorAxis.set(0,0,0);
    _minorRazorAxis.set(0,0,0);

    // Switch to the desired type of thrust finder
    if (_typeOfRazorFinder == 1)
    {
        cout << "type of Thrust Razor Finder = 1 : Tasso Thrust Razor " << endl; 
        TassoRazor();
        cout << "type of Thrust Razor Finder = 1 : Tasso Thrust Razor" << endl;
    }
    else if (_parp->GetEntries()<=1) // edited for TObjArray function GetEntries()
    {
        parpCheck += " ";
        parpCheck += std::to_string(_nEvt);
        parpCheck += " ";  
        TassoRazor();
    }
    else if (_typeOfRazorFinder == 2)
    {
        cout << "type of Thrust Razor Finder = 2 : Jetset Thrust Razor" << endl; 
        JetsetRazor();
    }
    // ###write
    //    evt->parameters().setValue("thrust",_principleRazorValue);

    FloatVec thrax;
    thrax.clear();
    thrax.push_back(_principleRazorAxis.x());
    thrax.push_back(_principleRazorAxis.y());
    thrax.push_back(_principleRazorAxis.z());

    _inParVec->parameters().setValue("principleRazorValue",_principleRazorValue);
    _inParVec->parameters().setValues("principleRazorAxis",thrax);

    if (_typeOfRazorFinder == 2)
    {
        thrax.clear();
        thrax.push_back(_majorRazorAxis.x());
        thrax.push_back(_majorRazorAxis.y());
        thrax.push_back(_majorRazorAxis.z());

        _inParVec->parameters().setValue("majorRazorValue",_majorRazorValue);
        _inParVec->parameters().setValues("majorRazorAxis",thrax);

        thrax.clear();
        thrax.push_back(_minorRazorAxis.x());
        thrax.push_back(_minorRazorAxis.y());
        thrax.push_back(_minorRazorAxis.z());

        _inParVec->parameters().setValue("minorRazorValue",_minorRazorValue);
        _inParVec->parameters().setValues("minorRazorAxis",thrax);

        float Oblateness;
        Oblateness = _majorRazorValue - _minorRazorValue;
        _inParVec->parameters().setValue("Oblateness",Oblateness);
        if ( (_majorRazorValue < 0) || (_minorRazorValue < 0) )
        {
            _inParVec->parameters().setValue("Oblateness",-1);
        }
    }


    // these are the final values:
    streamlog_out( DEBUG4 ) << " thrust: " << _principleRazorValue << " TV: " << _principleRazorAxis << endl;
    streamlog_out( DEBUG4 ) << "  major: " << _majorRazorValue << " TV: " << _majorRazorAxis << endl;
    streamlog_out( DEBUG4 ) << "  minor: " << _minorRazorValue << " TV: " << _minorRazorAxis << endl;
    //cout << "EVENT: " << _nEvt << endl;
    cout << " thrust: " << _principleRazorValue << " TV: " << _principleRazorAxis << endl;
    cout <<"                       "<< _principleRazorAxis.x()<<","<< _principleRazorAxis.y()<< ","<<_principleRazorAxis.z()<<endl;
    cout << "  major: " << _majorRazorValue << " TV: " << _majorRazorAxis << endl;
    cout << "  minor: " << _minorRazorValue << " TV: " << _minorRazorAxis << endl;

    double ptaX = _principleRazorAxis.x();
    double ptaY = _principleRazorAxis.y();
    double ptaZ = _principleRazorAxis.z();

    double vec[2][3][4]; // jet 1, jet 2 : true detectable, detected : energy, momx, momy, momz
    double Rvec[3][4]; // true, detectable, detected : energy, px, py, pz 
    int d = _thrustDetectability;
    double R;
    double beta2;
    if(_principleRazorValue==-1){
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
                const double* parp = aPart->getMomentum();
                double par4p[4];

                par4p[0] = aPart->getEnergy(); 
                par4p[1] = parp[0];
                par4p[2] = parp[1]; 
                par4p[3] = parp[2];   
                double pta[3] = {ptaX, ptaY, ptaZ};

                cout << "id      : " << id << endl;  
                cout << "Momentum: " << parp[0] <<" "<< parp[1] <<" "<< parp[2]<< endl;
                cout << "Thrust A: " << ptaX << " "<< ptaY << " " << ptaZ << endl;
                double dot = ptaX*parp[0]+ptaY*parp[1]+ptaZ*parp[2];
                cout << "dot " << dot << endl;

                // need momentum and energy of entire jet 

                bool isDarkMatter = (id == 1000022);
                bool isNeutrino = (
                        id == 12 || id == -12 ||
                        id == 14 || id == -14 ||
                        id == 16 || id == -16 ||
                        id == 18 || id == -18);
                double cos = parp[2]/(sqrt(parp[0]*parp[0]+parp[1]*parp[1]+parp[2]*parp[2]));
                bool isForward = ( cos > 0.9 || cos < - 0.9);
                int i; // jet #
                cout << "dot: " <<dot << endl; 
                if(dot>=0){i=0;}
                if(dot<0){i=1;}
                if (!isDarkMatter){
                    cout << "i: "<< i << endl;  
                    vec[i][0][0]+= par4p[0]; 
                    vec[i][0][1]+= par4p[1];
                    vec[i][0][2]+= par4p[2];
                    vec[i][0][3]+= par4p[3];
                    if (!isNeutrino){
                        vec[i][1][0]+= par4p[0];
                        vec[i][1][1]+= par4p[1];
                        vec[i][1][2]+= par4p[2];
                        vec[i][1][3]+= par4p[3];
                        if(!isForward){
                            vec[i][2][0]+= par4p[0];
                            vec[i][2][1]+= par4p[1];
                            vec[i][2][2]+= par4p[2];
                            vec[i][2][3]+= par4p[3];
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
            cout << "id: "<<id<<endl;
            cout << parp<< endl;
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
        if(parp1){
            cout << "there are valid beta events that have only 1 particle" << endl; 
        }
        if(!parp1){
            cout << "there are valid beta events that have not only 1 particle" <<endl; 
        }
    }
    else{
        if(d==0){_R_T->Fill(-2);}
        if(d==1){_R_DAB->Fill(-2);}
        if(d==2){_R_DED->Fill(-2);}
        
        if(parp1){
            cout << "there are inval beta events that have only 1 particle"<< endl;
        }
        if(!parp1){
            cout << "there are inval beta events that have not only 1 particle"<< endl;
            if(parp0){
                cout <<" zero particles" << endl;
            }
            if(!parp0){
                cout << "not zero particles"<< endl;
            }
        }
        
        betaCheck += " ";
        betaCheck += std::to_string(_nEvt);
        betaCheck += " ";  
    } 

    if (_principleRazorValue >= _max) _max = _principleRazorValue;
    if (_principleRazorValue <= _min) _min = _principleRazorValue;
    cout << "End EVENT "<< _nEvt<< endl;

    _nEvt ++ ; // different from original-moved out of for loop - summer 
}




void Razor::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void Razor::end(){ 
    _rootfile->Write();
    cout << parpCheck << endl;
    cout << betaCheck << endl; 
}

