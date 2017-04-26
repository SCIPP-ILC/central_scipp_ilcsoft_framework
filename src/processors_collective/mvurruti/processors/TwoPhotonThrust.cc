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
 * Using thrust reconstruction algorithm by: _____
 */

#include "TwoPhotonThrust.h"
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

static TH1F* _TV_T;
static TH1F* _TV_DAB;
static TH1F* _TV_DED;

static TH1F* _TA_T;
static TH1F* _TA_DAB;
static TH1F* _TA_DED;

TwoPhotonThrust TwoPhotonThrust;

TwoPhotonThrust::TwoPhotonThrust() : Processor("TwoPhotonThrust") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
    registerProcessorParameter( "typeOfTwoPhotonThrustFinder" ,
            "Type of thrust reconstruction algorithm to be used:\n#\t1 : Tasso algorithm\n#\t2 : JetSet algorithm"  ,
            _typeOfTwoPhotonThrustFinder , 2 ) ;
    registerProcessorParameter( "ThrustDetectability" ,
            "Detectability Level of the Thrust and Thrust Axis:\n#\t0 : True\n#\t1 : Detectable\n#\t2 : Detected" ,
            _ThrustDetectability, 0  );
}



void TwoPhotonThrust::init() { 
    streamlog_out(DEBUG)  << "   init called  " << std::endl ;


    if(_ThrustDetectability==0){
        _rootfile = new TFile("TwoPhotonThrust_eW.pW.I39212._T.root","RECREATE");
        _TV_T = new TH1F("TV_T", "Thrust Value",100,0,1); 
        _TA_T = new TH1F("TA_T", "Cosine(Thrust Angle)",100,0,1);
    }
    if(_ThrustDetectability==1){
        _rootfile = new TFile("TwoPhotonThrust_eW.pW.I39212._DAB.root","RECREATE");
        _TV_DAB = new TH1F("TV_DAB", "Thrust Value",100,0,1);
        _TA_DAB = new TH1F("TA_DAB", "Cosine(Thrust Angle)",100,0,1);
    }
    if(_ThrustDetectability==2){
        _rootfile = new TFile("TwoPhotonThrust_eW.pW.I39212._DED.root","RECREATE");
        _TV_DED = new TH1F("TV_DED", "Thrust Value",100,0,1);
        _TA_DED = new TH1F("TA_DED", "Cosine(Thrust Angle)",100,0,1);
    } 

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



void TwoPhotonThrust::processRunHeader( LCRunHeader* run) { 
    //run->parameters().setValue("thrust",12300321);
    //    _nRun++ ; 

} 



void TwoPhotonThrust::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    _inParVec = evt->getCollection( _colName) ;
    //cout << " get collection " << endl;
    //cout << _inParVec->getNumberOfElements() << endl;
    if (!_partMom.empty()) _partMom.clear();

    MCParticle* high_e;
    MCParticle* high_p;

    int id, stat;
    for (int n=0; n<_inParVec->getNumberOfElements();n++)
     {
         MCParticle* part = dynamic_cast<MCParticle*>(_inParVec->getElementAt(n));
         id = part->getPDG();
         stat = part->getGeneratorStatus();
         if(stat == 1){
             if(id==11){
                 high_e = part;
             }
             if(id==-11){
                 high_p = part;
             }
         } //end final state 
     } //end for loop
    for (int n = 0; n<_inParVec->getNumberOfElements() ; n++){
        MCParticle* part = dynamic_cast<MCParticle*>(_inParVec->getElementAt(n));
        id = part->getPDG();
        stat = part->getGeneratorStatus();
        if(stat==1){
            if(id==11){
                if(part->getEnergy()>high_e->getEnergy()){
                    high_e = part;
                }
            }
            if(id==-11){
                if(part->getEnergy()>high_p->getEnergy()){
                    high_p = part;
                }
            }
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
            }
        } //stat==1
    } // for particle 
    _nEvt ++ ; // different from original-moved out of for loop - summer 
    //reset variables for output   
    _principleTwoPhotonThrustValue = -1;
    _majorTwoPhotonThrustValue     = -1;
    _minorTwoPhotonThrustValue     = -1;
    _principleTwoPhotonThrustAxis.set(0,0,0);
    _majorTwoPhotonThrustAxis.set(0,0,0);
    _minorTwoPhotonThrustAxis.set(0,0,0);

    // Switch to the desired type of thrust finder
    if (_typeOfTwoPhotonThrustFinder == 1)
    { 
        //TassoTwoPhotonThrust();
    }
    else if (_partMom.size()<=1)
    {
        cout << " part Mom size less or equal 1 " << endl; 
        //TassoTwoPhotonThrust();
    }
    else if (_typeOfTwoPhotonThrustFinder == 2)
    {
        JetsetTwoPhotonThrust();
    }
    // ###write
    //    evt->parameters().setValue("thrust",_principleTwoPhotonThrustValue);

    FloatVec thrax;
    thrax.clear();
    thrax.push_back(_principleTwoPhotonThrustAxis.x());
    thrax.push_back(_principleTwoPhotonThrustAxis.y());
    thrax.push_back(_principleTwoPhotonThrustAxis.z());

    _inParVec->parameters().setValue("principleTwoPhotonThrustValue",_principleTwoPhotonThrustValue);
    _inParVec->parameters().setValues("principleTwoPhotonThrustAxis",thrax);

    if (_typeOfTwoPhotonThrustFinder == 2)
    {
        thrax.clear();
        thrax.push_back(_majorTwoPhotonThrustAxis.x());
        thrax.push_back(_majorTwoPhotonThrustAxis.y());
        thrax.push_back(_majorTwoPhotonThrustAxis.z());

        _inParVec->parameters().setValue("majorTwoPhotonThrustValue",_majorTwoPhotonThrustValue);
        _inParVec->parameters().setValues("majorTwoPhotonThrustAxis",thrax);

        thrax.clear();
        thrax.push_back(_minorTwoPhotonThrustAxis.x());
        thrax.push_back(_minorTwoPhotonThrustAxis.y());
        thrax.push_back(_minorTwoPhotonThrustAxis.z());

        _inParVec->parameters().setValue("minorTwoPhotonThrustValue",_minorTwoPhotonThrustValue);
        _inParVec->parameters().setValues("minorTwoPhotonThrustAxis",thrax);

        float Oblateness;
        Oblateness = _majorTwoPhotonThrustValue - _minorTwoPhotonThrustValue;
        _inParVec->parameters().setValue("Oblateness",Oblateness);
        if ( (_majorTwoPhotonThrustValue < 0) || (_minorTwoPhotonThrustValue < 0) )
        {
            _inParVec->parameters().setValue("Oblateness",-1);
        }
    }

    streamlog_out( DEBUG4 ) << " thrust: " << _principleTwoPhotonThrustValue << " TV: " << _principleTwoPhotonThrustAxis << endl;
    streamlog_out( DEBUG4 ) << "  major: " << _majorTwoPhotonThrustValue << " TV: " << _majorTwoPhotonThrustAxis << endl;
    streamlog_out( DEBUG4 ) << "  minor: " << _minorTwoPhotonThrustValue << " TV: " << _minorTwoPhotonThrustAxis << endl;
    cout << "EVENT: " << _nEvt << endl;
    cout << " thrust: " << _principleTwoPhotonThrustValue << " TV: " << _principleTwoPhotonThrustAxis << endl;
    cout <<"                       "<< _principleTwoPhotonThrustAxis.x()<<","<< _principleTwoPhotonThrustAxis.y()<< ","<<_principleTwoPhotonThrustAxis.z()<<endl;
    cout << "  major: " << _majorTwoPhotonThrustValue << " TV: " << _majorTwoPhotonThrustAxis << endl;
    cout << "  minor: " << _minorTwoPhotonThrustValue << " TV: " << _minorTwoPhotonThrustAxis << endl;

    double ptaX = _principleTwoPhotonThrustAxis.x();
    double ptaY = _principleTwoPhotonThrustAxis.y();
    double ptaZ = _principleTwoPhotonThrustAxis.z();

    double pta[3] = {ptaX, ptaY, ptaZ};

    double cosTT = abs(ptaZ);

    if(_ThrustDetectability == 0){
        _TV_T->Fill(_principleTwoPhotonThrustValue);
        _TA_T->Fill(cosTT);
    }
    if(_ThrustDetectability == 1){
        _TV_DAB->Fill(_principleTwoPhotonThrustValue);
        _TA_DAB->Fill(cosTT);
    }
    if(_ThrustDetectability == 2){
        _TV_DED->Fill(_principleTwoPhotonThrustValue);
        _TA_DED->Fill(cosTT);
    }



    if (_principleTwoPhotonThrustValue >= _max) _max = _principleTwoPhotonThrustValue;
    if (_principleTwoPhotonThrustValue <= _min) _min = _principleTwoPhotonThrustValue;
    cout << "End EVENT " <<_nEvt << endl; 
}




void TwoPhotonThrust::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void TwoPhotonThrust::end(){ 
    _rootfile->Write();
}
int TwoPhotonThrust::TassoThrust(){
    int ThrustError = 0;
    if(_partMom.size()<=0)
    {
        //partMom0 = true;
        ThrustError = -1;
        cout<< "NO PARTICLES IN EVENT AAAAAH" << endl;
        _principleTwoPhotonThrustValue = -1;
        _principleTwoPhotonThrustAxis.set(0,0,0);
    }
    if(_partMom.size()==1)
    {
        cout << "partMom.size"<< _partMom.size() << endl;
        //partMom1 = true;
        _principleTwoPhotonThrustValue = 1;
        _principleTwoPhotonThrustAxis = _partMom[0];
    }
    return ThrustError;
}

int TwoPhotonThrust::JetsetTwoPhotonThrust(){
    const int nwork=11,iFastMax = 4,iGood=2;
    const float dConv=0.0001; // 0.0001
    int sgn;
    double theta=0,phi=0;
    double thp,thps,tds,tmax,dOblateness;
    vector<Hep3Vector> TAxes(3),Fast(iFastMax+1),Workv(nwork);
    vector<double> Workf(nwork),dTwoPhotonThrust(3);
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
        dTwoPhotonThrust[pass] = 0;
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
            if ( thp < dTwoPhotonThrust[pass] - dConv ) continue;
            if ( thp > dTwoPhotonThrust[pass] + dConv )
            {
                nagree = 0;
                //          if (myrnd.flat() > 0.49999)
                //        {sgn = 1;} else {sgn=-1;}
                sgn = 1;
                TAxes[pass] = sgn*tpr/(tmax*thp);
                dTwoPhotonThrust[pass] = thp;
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
    dTwoPhotonThrust[2] = thp/tmax;

    // Rotate back to original coordinate system.
    for ( unsigned int i = 0;i < TAxes.size(); i++)
    {
        TAxes[i].rotateY(theta);
        TAxes[i].rotateZ(phi);
    }
    dOblateness = dTwoPhotonThrust[1] - dTwoPhotonThrust[2];

    _principleTwoPhotonThrustValue = dTwoPhotonThrust[0];
    _majorTwoPhotonThrustValue     = dTwoPhotonThrust[1];
    _minorTwoPhotonThrustValue     = dTwoPhotonThrust[2];
    _principleTwoPhotonThrustAxis  =   TAxes[0];
    _majorTwoPhotonThrustAxis      =   TAxes[1];
    _minorTwoPhotonThrustAxis      =   TAxes[2];


    return  0;
}

//______________________________________________________________
// helper function to get sign of b
double TwoPhotonThrust::sign(double a, double b)
{
    if ( b < 0 )
    { return -fabs(a); } else { return fabs(a); }
}
//______________________________________________________________
double TwoPhotonThrust::min(double a, double b)
{
    if ( a < b )
    { return a; } else { return b; }
}
//______________________________________________________________


