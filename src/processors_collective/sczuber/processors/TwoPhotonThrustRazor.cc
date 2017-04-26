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
 * January 18, 2017 
 *
 * Using thrust reconstruction algorithm by: ____
 * Used to develope razor variables which here are applied to 
 * the hadronic system of two-photon events. 
 */

#include "TwoPhotonThrustRazor.h"
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

static TH1F* _R_T;
static TH1F* _R_DAB;
static TH1F* _R_DED;

TwoPhotonThrustRazor TwoPhotonThrustRazor;

TwoPhotonThrustRazor::TwoPhotonThrustRazor() : Processor("TwoPhotonThrustRazor") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
    registerProcessorParameter( "typeOfTwoPhotonThrustRazorFinder" ,
            "Type of thrust reconstruction algorithm to be used:\n#\t1 : Tasso algorithm\n#\t2 : JetSet algorithm"  ,
            _typeOfTwoPhotonThrustRazorFinder , 2 ) ;
    registerProcessorParameter( "ThrustDetectability" ,
            "Detectability Level of the Thrust and Thrust Axis:\n#\t0 : True\n#\t1 : Detectable\n#\t2 : Detected" ,
            _ThrustDetectability, 2  );
}



void TwoPhotonThrustRazor::init() { 
    streamlog_out(DEBUG)  << "   init called  " << std::endl ;


    if(_ThrustDetectability==0){_rootfile = new TFile("TwoPhotonThrustRazor_.eW.pW.I39212._T.root","RECREATE");
        _R_T = new TH1F("R_T", "R=MTR/MR",130,-3,10);}
    if(_ThrustDetectability==1){_rootfile = new TFile("TwoPhotonThrustRazor_.eW.pW.I39212._DAB.root","RECREATE");
        _R_DAB = new TH1F("R_DAB", "R=MTR/MR",130,-3,10);}
    if(_ThrustDetectability==2){_rootfile = new TFile("TwoPhotonThrustRazor_.eW.pW.I39212._DED.root","RECREATE");
        _R_DED = new TH1F("R_DED", "R=MTR/MR",130,-3,10);} 

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
    } // if file not existusually a good idea to
    //printParameters() ;
    _nEvt = 0 ;
}



void TwoPhotonThrustRazor::processRunHeader( LCRunHeader* run) { 
    //run->parameters().setValue("thrust",12300321);
    //    _nRun++ ; 

} 



void TwoPhotonThrustRazor::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    _inParVec = evt->getCollection( _colName) ;
    //cout << " get collection " << endl;
    //cout << _inParVec->getNumberOfElements() << endl;
    if (!_partMom.empty()) _partMom.clear();
    partMom_0 = false;
    partMom_1 = false; 

    MCParticle* high_e;
    MCParticle* high_p; 

    int id, stat; 

    // setting high_e and high_p to a final state electron and positron
    cout << "Loop #1: set high_e and high_p " << endl; 
    for(int n=0; n<_inParVec->getNumberOfElements(); n++)
    {
        MCParticle* hit = dynamic_cast<MCParticle*>(_inParVec->getElementAt(n));
        id = hit->getPDG();
        stat = hit->getGeneratorStatus();
        if(stat==1){
            if(id==11){
                high_e = hit;
            }
            if(id==-11){
                high_p = hit;
            }
        }// end final state 
    } // end for loop
    cout<< "Loop #2: set high_e and high_p to highest energy ones" << endl; 
    for (int n = 0; n<_inParVec->getNumberOfElements(); n++)
    {
        MCParticle* hit = dynamic_cast<MCParticle*>(_inParVec->getElementAt(n));
        id = hit->getPDG();
        stat = hit->getGeneratorStatus();
        if(stat==1){
            if(id==11){
                if(hit->getEnergy()>high_e->getEnergy()){
                    high_e = hit;
                }
            }
            if(id == -11){
                if(hit->getEnergy()>high_p->getEnergy()){
                    high_p = hit;
                }
            }
        }// end final state 

    }//end for loop 
    cout << " Start Loop #3" << endl; 
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
            const double* partMom = aPart->getMomentum();
            double momMag = sqrt(partMom[0]*partMom[0]+partMom[1]*partMom[1]+partMom[2]*partMom[2]);
            double cos = partMom[2]/momMag;
            bool isForward = (cos > 0.9 || cos < -0.9);
            bool isDarkMatter = (id == 1000022);
            bool isNeutrino = (
                    id == 12 || id == -12 ||
                    id == 14 || id == -14 ||
                    id == 16 || id == -16 ||
                    id == 18 || id == -18);
            bool isDetectable = (!isDarkMatter && !isNeutrino);
            bool isDetected = (isDetectable && !isForward);
            if(aPart != high_e && aPart != high_p){ 

                if(_ThrustDetectability == 0){ 
                    _partMom.push_back( Hep3Vector(partMom[0], partMom[1], partMom[2]) ); 
                }
                if(_ThrustDetectability == 1){
                    if(isDetectable){
                        _partMom.push_back(Hep3Vector(partMom[0], partMom[1], partMom[2]));
                    }
                }
                if(_ThrustDetectability == 2){
                    if(isDetected){
                        _partMom.push_back(Hep3Vector(partMom[0], partMom[1], partMom[2]));
                    }
                }
            } // end particle not original e+e-
        } //stat==1
    } // for particle 
    _nEvt ++ ; // different from original-moved out of for loop - summer 
    //reset variables for output   
    _principleTwoPhotonThrustRazorValue = -1;
    _majorTwoPhotonThrustRazorValue     = -1;
    _minorTwoPhotonThrustRazorValue     = -1;
    _principleTwoPhotonThrustRazorAxis.set(0,0,0);
    _majorTwoPhotonThrustRazorAxis.set(0,0,0);
    _minorTwoPhotonThrustRazorAxis.set(0,0,0);

    // Switch to the desired type of thrust finder
    if (_typeOfTwoPhotonThrustRazorFinder == 1)
    { 
        TassoThrustRazor();
    }
    else if (_partMom.size()<=1)
    {
        cout << " part Mom size less or equal 1 " << endl; 
        TassoThrustRazor();
    }
    else if (_typeOfTwoPhotonThrustRazorFinder == 2)
    {
        JetsetTwoPhotonThrustRazor();
    }
    // ###write
    //    evt->parameters().setValue("thrust",_principleTwoPhotonThrustRazorValue);

    FloatVec thrax;
    thrax.clear();
    thrax.push_back(_principleTwoPhotonThrustRazorAxis.x());
    thrax.push_back(_principleTwoPhotonThrustRazorAxis.y());
    thrax.push_back(_principleTwoPhotonThrustRazorAxis.z());

    _inParVec->parameters().setValue("principleTwoPhotonThrustRazorValue",_principleTwoPhotonThrustRazorValue);
    _inParVec->parameters().setValues("principleTwoPhotonThrustRazorAxis",thrax);

    if (_typeOfTwoPhotonThrustRazorFinder == 2)
    {
        thrax.clear();
        thrax.push_back(_majorTwoPhotonThrustRazorAxis.x());
        thrax.push_back(_majorTwoPhotonThrustRazorAxis.y());
        thrax.push_back(_majorTwoPhotonThrustRazorAxis.z());

        _inParVec->parameters().setValue("majorTwoPhotonThrustRazorValue",_majorTwoPhotonThrustRazorValue);
        _inParVec->parameters().setValues("majorTwoPhotonThrustRazorAxis",thrax);

        thrax.clear();
        thrax.push_back(_minorTwoPhotonThrustRazorAxis.x());
        thrax.push_back(_minorTwoPhotonThrustRazorAxis.y());
        thrax.push_back(_minorTwoPhotonThrustRazorAxis.z());

        _inParVec->parameters().setValue("minorTwoPhotonThrustRazorValue",_minorTwoPhotonThrustRazorValue);
        _inParVec->parameters().setValues("minorTwoPhotonThrustRazorAxis",thrax);

        float Oblateness;
        Oblateness = _majorTwoPhotonThrustRazorValue - _minorTwoPhotonThrustRazorValue;
        _inParVec->parameters().setValue("Oblateness",Oblateness);
        if ( (_majorTwoPhotonThrustRazorValue < 0) || (_minorTwoPhotonThrustRazorValue < 0) )
        {
            _inParVec->parameters().setValue("Oblateness",-1);
        }
    }

    streamlog_out( DEBUG4 ) << " thrust: " << _principleTwoPhotonThrustRazorValue << " TV: " << _principleTwoPhotonThrustRazorAxis << endl;
    streamlog_out( DEBUG4 ) << "  major: " << _majorTwoPhotonThrustRazorValue << " TV: " << _majorTwoPhotonThrustRazorAxis << endl;
    streamlog_out( DEBUG4 ) << "  minor: " << _minorTwoPhotonThrustRazorValue << " TV: " << _minorTwoPhotonThrustRazorAxis << endl;
    cout << "EVENT: " << _nEvt << endl;
    cout << " thrust: " << _principleTwoPhotonThrustRazorValue << " TV: " << _principleTwoPhotonThrustRazorAxis << endl;
    cout <<"                       "<< _principleTwoPhotonThrustRazorAxis.x()<<","<< _principleTwoPhotonThrustRazorAxis.y()<< ","<<_principleTwoPhotonThrustRazorAxis.z()<<endl;
    cout << "  major: " << _majorTwoPhotonThrustRazorValue << " TV: " << _majorTwoPhotonThrustRazorAxis << endl;
    cout << "  minor: " << _minorTwoPhotonThrustRazorValue << " TV: " << _minorTwoPhotonThrustRazorAxis << endl;

    double ptaX = _principleTwoPhotonThrustRazorAxis.x();
    double ptaY = _principleTwoPhotonThrustRazorAxis.y();
    double ptaZ = _principleTwoPhotonThrustRazorAxis.z();

    double pta[3] = {ptaX, ptaY, ptaZ};

    double vec[2][3][4]; // jet 1, jet 2 : true detectable, detected : energy, momx, momy, momz 
    double Rvec[3][4]; // true, detectable, detected : energy, px, py, pz 
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
            if(aPart != high_e && aPart != high_p){
                const double* partMom = aPart->getMomentum();
                double part4mom[4] = {aPart->getEnergy(), partMom[0], partMom[1], partMom[2]} ;

                cout << "id      : " << id << endl;  
                cout << "Momentum: " << partMom[0] <<" "<< partMom[1] <<" "<< partMom[2]<< endl;

                double dot = ptaX*partMom[0]+ptaY*partMom[1]+ptaZ*partMom[2]; 

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
                if(dot>0){i=0;}
                if(dot<0){i=1;}
                if(dot == 0){i=0;} 

                vec[i][0][0]+= part4mom[0]; 
                vec[i][0][1]+= part4mom[1];
                vec[i][0][2]+= part4mom[2];
                vec[i][0][3]+= part4mom[3];
                if (!isDarkMatter && !isNeutrino){
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
            }// end particle not original 
        } //stat == 1 
    } // for particle
    int d = _ThrustDetectability;
    //cout << "d "<<d<<endl; 
    double beta = (vec[0][d][0]-vec[1][d][0])/(vec[0][d][3]-vec[1][d][3]); // beta using detectable particles
    //cout << "BETA " << beta << endl;  
    double beta2 = pow(beta,2);
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
        double cos = part4Vec[3]/sqrt(part4Vec[1]*part4Vec[1]+part4Vec[2]*part4Vec[2]+part4Vec[3]*part4Vec[3]);
        bool isForward = (cos > 0.9 || cos < -0.9);
        if(stat ==1){
            if(aPart != high_e && aPart != high_p){
                Rvec[0][0]+=R4Vec[0];
                Rvec[0][1]+=R4Vec[1];
                Rvec[0][2]+=R4Vec[2];
                Rvec[0][3]+=R4Vec[3];
                if(!isDarkMatter && !isNeutrino){
                    Rvec[1][0]+=R4Vec[0];
                    Rvec[1][1]+=R4Vec[1];
                    Rvec[1][2]+=R4Vec[2];
                    Rvec[1][3]+=R4Vec[3];
                    if(!isForward){
                        Rvec[2][0]+=R4Vec[0];
                        Rvec[2][1]+=R4Vec[1];
                        Rvec[2][2]+=R4Vec[2];
                        Rvec[2][3]+=R4Vec[3];

                    }
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
    double R = MTR/(2*pj1);
    if(beta2<=1){
        if(d==0){_R_T->Fill(R);}
        if(d==1){_R_DAB->Fill(R);}
        if(d==2){_R_DED->Fill(R);}
    }
    if(beta2>1){
        if(d==0){_R_T  ->Fill(-2);}
        if(d==1){_R_DAB->Fill(-2);}
        if(d==2){_R_DED->Fill(-2);}
    }
    if(partMom_0==true){
        if(d==0){_R_T  ->Fill(-.5);}
        if(d==1){_R_DAB->Fill(-.5);}
        if(d==2){_R_DED->Fill(-.5);}
    }
    if(partMom_1==true){
        if(d==0){_R_T  ->Fill(-1.5);}
        if(d==1){_R_DAB->Fill(-1.5);}
        if(d==2){_R_DED->Fill(-1.5);}
    }

    if (_principleTwoPhotonThrustRazorValue >= _max) _max = _principleTwoPhotonThrustRazorValue;
    if (_principleTwoPhotonThrustRazorValue <= _min) _min = _principleTwoPhotonThrustRazorValue;
}




void TwoPhotonThrustRazor::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void TwoPhotonThrustRazor::end(){ 
    _rootfile->Write();
}

int TwoPhotonThrustRazor::TassoThrustRazor(){
    int ThrustError = 0;
    if(_partMom.size()<=0)
    {
        ThrustError = -1;
        cout<< "NO PARTICLES IN EVENT AAAAAH" << endl;
        _principleTwoPhotonThrustRazorValue = -1;
        _principleTwoPhotonThrustRazorAxis.set(0,0,0);
        partMom_0 = true;
    }
    if(_partMom.size()==1)
    {
        cout << "partMom.size"<< _partMom.size() << endl;
        _principleTwoPhotonThrustRazorValue = 1;
        _principleTwoPhotonThrustRazorAxis = _partMom[0];
        partMom_1 = true; 
    }
    return ThrustError;
}

int TwoPhotonThrustRazor::JetsetTwoPhotonThrustRazor(){
    const int nwork=11,iFastMax = 4,iGood=2;
    const float dConv=0.0001; // 0.0001
    int sgn;
    double theta=0,phi=0;
    double thp,thps,tds,tmax,dOblateness;
    vector<Hep3Vector> TAxes(3),Fast(iFastMax+1),Workv(nwork);
    vector<double> Workf(nwork),dTwoPhotonThrustRazor(3);
    Hep3Vector tdi,tpr,mytest;

    tmax = 0;
    for ( unsigned int i=0; i < _partMom.size(); i++)
        tmax += _partMom[i].mag();

    // pass = 0: find thrust axis
    // pass = 1: find major axis
    for ( int pass=0; pass <= 1; pass++ )
    {
        if ( pass == 1 )
        {
            phi   = TAxes[0].phi();
            theta = TAxes[0].theta();
            for ( unsigned  int i = 0;i < _partMom.size(); i++)
            {
                _partMom[i].rotateZ(-phi);
                _partMom[i].rotateY(-theta);
            }
            TAxes[0].set(0,0,1);
        } // if pass == 1

        // Find the ifast highest momentum particles and
        // put the highest in Fast[0], next in Fast[1],....Fast[iFast-1].
        // Fast[iFast] is just a workspace.
        for ( unsigned  int i = 0; i < Fast.size(); i++ )
            Fast[i].set(0,0,0);

        for ( unsigned int i = 0; i < _partMom.size(); i++ )
        {
            for ( int ifast = iFastMax -1; ifast >= 0 ; ifast-- )
            {
                if (_partMom[i].mag2() > Fast[ifast].mag2() )
                {
                    Fast[ifast + 1] = Fast[ifast];
                    if (ifast == 0) Fast[ifast] = _partMom[i];
                }
                else
                {
                    Fast[ifast + 1] = _partMom[i];
                    break;
                } // if p>p_fast
            } // for ifast 
        } // for i 

        // Find axis with highest thrust (case 0)/ highest major (case 1).
        for ( unsigned int iw = 0; iw < Workv.size(); iw++ )
        {
            Workf[iw] = 0.;
        }
        int p = (int) min( iFastMax,_partMom.size() ) - 1 ;
        int nc = 1 << p;
        for ( int n = 0; n < nc; n++ )
        {
            tdi.set(0,0,0);
            for (int i = 0; i < min(iFastMax,nc) ; i++)
            {
                if ( (1 << (i+1)) * ( (n + (1<<i)) / (1<<(i+1)) ) >= n+1) //i+1 
                { sgn = -1;} else {sgn = 1;}
                tdi += sgn*Fast[i];
                if (pass==1) tdi.setZ(0);
            } // for i 
            tds = tdi.mag2();
            for ( int iw = (int) min(n,9); iw >= 0; iw-- )
            {
                if (tds > Workf[iw])
                {
                    Workf[iw+1] = Workf[iw];
                    Workv[iw+1] = Workv[iw];
                    if (iw == 0)
                    { Workv[iw] = tdi; Workf[iw] = tds;}
                }
                else // if tds 
                {
                    Workv[iw+1] = tdi;
                    Workf[iw+1] = tds;
                } // if tds 
            } // for iw
        } // for n 

        // Iterate direction of axis until stable maximum.
        dTwoPhotonThrustRazor[pass] = 0;
        int nagree = 0;
        for ( int iw = 0; iw < min(nc,10) && nagree < iGood; iw++ )
        {
            thp = 0;
            thps = -99999.;
            while ( thp > thps + dConv )
            {
                thps = thp;
                if ( thp <= 1E-10 )
                { tdi = Workv[iw]; } else { tdi=tpr; }
                tpr.set(0,0,0);
                for ( unsigned int i = 0; i < _partMom.size(); i++ )
                {
                    sgn = (int) sign(1,tdi.dot(_partMom[i]));
                    tpr += sgn*_partMom[i];
                    if (pass == 1) { tpr.setZ(0); } // ###
                } // for i 
                thp = tpr.mag()/tmax;
            } // while 
            // Save good axis. Try new initial axis until enough
            // tries agree.
            if ( thp < dTwoPhotonThrustRazor[pass] - dConv ) continue;
            if ( thp > dTwoPhotonThrustRazor[pass] + dConv )
            {
                nagree = 0;
                //          if (myrnd.flat() > 0.49999)
                //        {sgn = 1;} else {sgn=-1;}
                sgn = 1;
                TAxes[pass] = sgn*tpr/(tmax*thp);
                dTwoPhotonThrustRazor[pass] = thp;
            } // if thp
            nagree++;
        } // for iw (2)
    } // for pass ...
    // Find minor axis and value by orthogonality.
    if (myrnd.flat() > 0.49999)
    {sgn = 1;} else {sgn=-1;}
    TAxes[2].set( -sgn*TAxes[1].y(), sgn*TAxes[1].x(), 0);
    thp = 0.;
    for ( unsigned int i = 0; i < _partMom.size(); i++ )
    {
        thp += fabs(TAxes[2].dot(_partMom[i]) );
    } // for i 
    dTwoPhotonThrustRazor[2] = thp/tmax;

    // Rotate back to original coordinate system.
    for ( unsigned int i = 0;i < TAxes.size(); i++)
    {
        TAxes[i].rotateY(theta);
        TAxes[i].rotateZ(phi);
    }
    dOblateness = dTwoPhotonThrustRazor[1] - dTwoPhotonThrustRazor[2];

    _principleTwoPhotonThrustRazorValue = dTwoPhotonThrustRazor[0];
    _majorTwoPhotonThrustRazorValue     = dTwoPhotonThrustRazor[1];
    _minorTwoPhotonThrustRazorValue     = dTwoPhotonThrustRazor[2];
    _principleTwoPhotonThrustRazorAxis  =   TAxes[0];
    _majorTwoPhotonThrustRazorAxis      =   TAxes[1];
    _minorTwoPhotonThrustRazorAxis      =   TAxes[2];


    return  0;
}

//______________________________________________________________
// helper function to get sign of b
double TwoPhotonThrustRazor::sign(double a, double b)
{
    if ( b < 0 )
    { return -fabs(a); } else { return fabs(a); }
}
//______________________________________________________________
double TwoPhotonThrustRazor::min(double a, double b)
{
    if ( a < b )
    { return a; } else { return b; }
}
//______________________________________________________________


