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
 * Using thrust reconstruction code by: _____
 */

#include "ThrustRazor.h"
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

static TH1F* _RPlot;

ThrustRazor ThrustRazor;

ThrustRazor::ThrustRazor() : Processor("ThrustRazor") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
    registerProcessorParameter( "typeOfThrustRazorFinder" ,
            "Type of thrust reconstruction algorithm to be used:\n#\t1 : Tasso algorithm\n#\t2 : JetSet algorithm"  ,
            _typeOfThrustRazorFinder , 2 ) ;
}



void ThrustRazor::init() { 
    streamlog_out(DEBUG)  << "   init called  " << std::endl ;

    _rootfile = new TFile("ThrustRazor_.39133.root","RECREATE");
    _RPlot = new TH1F("RPlot", "R =MTR/MR",100,0,10);
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



void ThrustRazor::processRunHeader( LCRunHeader* run) { 
    //run->parameters().setValue("thrust",12300321);
    //    _nRun++ ;
} 



void ThrustRazor::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    _inParVec = evt->getCollection( _colName) ;
    cout << _inParVec->getNumberOfElements() << endl;
    if (!_partMom.empty()) _partMom.clear();

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
        if(stat==1){
            bool isDarkMatter = (id == 1000022);
            bool isNeutrino = (
                    id == 12 || id == -12 ||
                    id == 14 || id == -14 ||
                    id == 16 || id == -16 ||
                    id == 18 || id == -18);

            if(!isDarkMatter && !isNeutrino){  
                const double* partMom = aPart->getMomentum();
                _partMom.push_back( Hep3Vector(partMom[0], partMom[1], partMom[2]) ); 
            }
        }
    }
    _nEvt ++ ; // different from original-moved out of for loop - summer 
    //reset variables for output   
    _principleThrustRazorValue = -1;
    _majorThrustRazorValue     = -1;
    _minorThrustRazorValue     = -1;
    _principleThrustRazorAxis.set(0,0,0);
    _majorThrustRazorAxis.set(0,0,0);
    _minorThrustRazorAxis.set(0,0,0);

    // Switch to the desired type of thrust finder
    if (_typeOfThrustRazorFinder == 1)
    { 
        TassoThrustRazor();
    }
    else if (_partMom.size()<=1)
    { 
        TassoThrustRazor();
    }
    else if (_typeOfThrustRazorFinder == 2)
    {
        JetsetThrustRazor();
    }
    // ###write
    //    evt->parameters().setValue("thrust",_principleThrustRazorValue);

    FloatVec thrax;
    thrax.clear();
    thrax.push_back(_principleThrustRazorAxis.x());
    thrax.push_back(_principleThrustRazorAxis.y());
    thrax.push_back(_principleThrustRazorAxis.z());

    _inParVec->parameters().setValue("principleThrustRazorValue",_principleThrustRazorValue);
    _inParVec->parameters().setValues("principleThrustRazorAxis",thrax);

    if (_typeOfThrustRazorFinder == 2)
    {
        thrax.clear();
        thrax.push_back(_majorThrustRazorAxis.x());
        thrax.push_back(_majorThrustRazorAxis.y());
        thrax.push_back(_majorThrustRazorAxis.z());

        _inParVec->parameters().setValue("majorThrustRazorValue",_majorThrustRazorValue);
        _inParVec->parameters().setValues("majorThrustRazorAxis",thrax);

        thrax.clear();
        thrax.push_back(_minorThrustRazorAxis.x());
        thrax.push_back(_minorThrustRazorAxis.y());
        thrax.push_back(_minorThrustRazorAxis.z());

        _inParVec->parameters().setValue("minorThrustRazorValue",_minorThrustRazorValue);
        _inParVec->parameters().setValues("minorThrustRazorAxis",thrax);

        float Oblateness;
        Oblateness = _majorThrustRazorValue - _minorThrustRazorValue;
        _inParVec->parameters().setValue("Oblateness",Oblateness);
        if ( (_majorThrustRazorValue < 0) || (_minorThrustRazorValue < 0) )
        {
            _inParVec->parameters().setValue("Oblateness",-1);
        }
    }

    streamlog_out( DEBUG4 ) << " thrust: " << _principleThrustRazorValue << " TV: " << _principleThrustRazorAxis << endl;
    streamlog_out( DEBUG4 ) << "  major: " << _majorThrustRazorValue << " TV: " << _majorThrustRazorAxis << endl;
    streamlog_out( DEBUG4 ) << "  minor: " << _minorThrustRazorValue << " TV: " << _minorThrustRazorAxis << endl;
    cout << "EVENT: " << _nEvt << endl;
    cout << " thrust: " << _principleThrustRazorValue << " TV: " << _principleThrustRazorAxis << endl;
    cout <<"                       "<< _principleThrustRazorAxis.x()<<","<< _principleThrustRazorAxis.y()<< ","<<_principleThrustRazorAxis.z()<<endl;
    cout << "  major: " << _majorThrustRazorValue << " TV: " << _majorThrustRazorAxis << endl;
    cout << "  minor: " << _minorThrustRazorValue << " TV: " << _minorThrustRazorAxis << endl;

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
            const double* partMom = aPart->getMomentum();
            double part4mom[4] = {aPart->getEnergy(), partMom[0], partMom[1], partMom[2]};

            double ptaX = _principleThrustRazorAxis.x();
            double ptaY = _principleThrustRazorAxis.y();
            double ptaZ = _principleThrustRazorAxis.z();

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
            if(dot>0){
                vec[0][0][0]+= part4mom[0]; 
                vec[0][0][1]+= part4mom[1];
                vec[0][0][2]+= part4mom[2];
                vec[0][0][3]+= part4mom[3];
                if (!isDarkMatter && !isNeutrino){
                    vec[0][1][0]+= part4mom[0];
                    vec[0][1][1]+= part4mom[1];
                    vec[0][1][2]+= part4mom[2];
                    vec[0][1][3]+= part4mom[3];
                    if(!isForward){
                        vec[0][2][0]+= part4mom[0];
                        vec[0][2][1]+= part4mom[1];
                        vec[0][2][2]+= part4mom[2];
                        vec[0][2][3]+= part4mom[3];
                    }
                }
            }
            if(dot<0){ 
                vec[1][0][0]+= part4mom[0]; 
                vec[1][0][1]+= part4mom[1];
                vec[1][0][2]+= part4mom[2];
                vec[1][0][3]+= part4mom[3];
                if (!isDarkMatter && !isNeutrino){
                    vec[1][1][0]+=part4mom[0];
                    vec[1][1][1]+=part4mom[1];
                    vec[1][1][2]+=part4mom[2];
                    vec[1][1][3]+=part4mom[3];
                    if(!isForward){
                        vec[1][2][0]+=part4mom[0];
                        vec[1][2][1]+=part4mom[1];
                        vec[1][2][2]+=part4mom[2];
                        vec[1][2][3]+=part4mom[3];
                    }
                }
            }
            if(dot==0){
                int jet = rand() % 2;
                cout << " RANDOM " << jet << endl; 
                vec[1][0][0]+= part4mom[0]; 
                vec[1][0][1]+= part4mom[1];
                vec[1][0][2]+= part4mom[2];
                vec[1][0][3]+= part4mom[3];
                if (!isDarkMatter && !isNeutrino){
                    vec[1][1][0]+=part4mom[0];
                    vec[1][1][1]+=part4mom[1];
                    vec[1][1][2]+=part4mom[2];
                    vec[1][1][3]+=part4mom[3];
                    if(!isForward){
                        vec[1][2][0]+=part4mom[0];
                        vec[1][2][1]+=part4mom[1];
                        vec[1][2][2]+=part4mom[2];
                        vec[1][2][3]+=part4mom[3];
                    }
                }
            }
        }
    }
    double beta = (vec[0][1][0]-vec[1][1][0])/(vec[0][1][3]-vec[1][1][3]); // beta using detectable particles
    cout << "BETA " << beta << endl;  
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
        if(stat ==1){
            if(!isDarkMatter && !isNeutrino){
                Rvec[1][0]+=R4Vec[0];
                Rvec[1][1]+=R4Vec[1];
                Rvec[1][2]+=R4Vec[2];
                Rvec[1][3]+=R4Vec[3];
            }
        }


    }
    double ETM[2] = {-Rvec[1][1], - Rvec[1][2]};
    double ETMmag = sqrt(ETM[0]*ETM[0]+ETM[1]*ETM[1]);

    double ptj1mag = sqrt(vec[0][1][1]*vec[0][1][1]+vec[0][1][2]*vec[0][1][2]);
    double ptj2mag = sqrt(vec[1][1][1]*vec[1][1][1]+vec[1][1][2]*vec[1][1][2]);
    double ptmagsum = ptj1mag + ptj2mag; 

    double ptvecsum[2] = {vec[0][1][1]+vec[1][1][1], vec[0][1][2]+vec[1][1][2]};
    double ETMdotptvecsum = ETM[0]*ptvecsum[0]+ETM[1]*ptvecsum[1];
    double MTR = sqrt((ETMmag*ptmagsum-ETMdotptvecsum)/2);

    double pj1 = sqrt(vec[0][1][1]*vec[0][1][1]+vec[0][1][2]*vec[0][1][2]+vec[0][1][3]*vec[0][1][3]);
    double R = MTR/(2*pj1);
    if(beta2<=1){
        _RPlot->Fill(R);
    }

    if (_principleThrustRazorValue >= _max) _max = _principleThrustRazorValue;
    if (_principleThrustRazorValue <= _min) _min = _principleThrustRazorValue;
}




void ThrustRazor::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void ThrustRazor::end(){ 
    _rootfile->Write();
}

int ThrustRazor::JetsetThrustRazor(){
    const int nwork=11,iFastMax = 4,iGood=2;
    const float dConv=0.0001; // 0.0001
    int sgn;
    double theta=0,phi=0;
    double thp,thps,tds,tmax,dOblateness;
    vector<Hep3Vector> TAxes(3),Fast(iFastMax+1),Workv(nwork);
    vector<double> Workf(nwork),dThrustRazor(3);
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
        dThrustRazor[pass] = 0;
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
            if ( thp < dThrustRazor[pass] - dConv ) continue;
            if ( thp > dThrustRazor[pass] + dConv )
            {
                nagree = 0;
                //          if (myrnd.flat() > 0.49999)
                //        {sgn = 1;} else {sgn=-1;}
                sgn = 1;
                TAxes[pass] = sgn*tpr/(tmax*thp);
                dThrustRazor[pass] = thp;
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
    dThrustRazor[2] = thp/tmax;

    // Rotate back to original coordinate system.
    for ( unsigned int i = 0;i < TAxes.size(); i++)
    {
        TAxes[i].rotateY(theta);
        TAxes[i].rotateZ(phi);
    }
    dOblateness = dThrustRazor[1] - dThrustRazor[2];

    _principleThrustRazorValue = dThrustRazor[0];
    _majorThrustRazorValue     = dThrustRazor[1];
    _minorThrustRazorValue     = dThrustRazor[2];
    _principleThrustRazorAxis  =   TAxes[0];
    _majorThrustRazorAxis      =   TAxes[1];
    _minorThrustRazorAxis      =   TAxes[2];


    return  0;
}

//______________________________________________________________
// helper function to get sign of b
double ThrustRazor::sign(double a, double b)
{
    if ( b < 0 )
    { return -fabs(a); } else { return fabs(a); }
}
//______________________________________________________________
double ThrustRazor::min(double a, double b)
{
    if ( a < b )
    { return a; } else { return b; }
}
//______________________________________________________________


