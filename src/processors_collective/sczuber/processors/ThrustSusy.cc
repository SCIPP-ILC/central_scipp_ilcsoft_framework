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
 * March 3, 2017 
 *
 * Using thrust reconstruction code by: _____
 */

#include "ThrustSusy.h"
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
static TH2F* _TVTA_T;
static TH2F* _TVTA_DAB;
static TH2F* _TVTA_DED;

ThrustSusy ThrustSusy;

ThrustSusy::ThrustSusy() : Processor("ThrustSusy") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
    registerProcessorParameter( "typeOfThrustSusyFinder" ,
            "Type of thrust reconstruction algorithm to be used:\n#\t1 : Tasso algorithm\n#\t2 : JetSet algorithm"  ,
            _typeOfThrustSusyFinder , 2 ) ;
    registerProcessorParameter( "thrustDetectability",
            "Detectability of the Thrust Axis/Value to be used:\n#\t0 : True \n#t1 : Detectable \n#t2 : Detected" ,
            _thrustDetectability, 2 );
    
}



void ThrustSusy::init() {
    cout << "entered init" << endl;  
    streamlog_out(DEBUG)  << "   init called  " << std::endl ;

    if(_thrustDetectability==0){

        _rootfile = new TFile("ThrustSusy_.39133._T.root","RECREATE");

        _TV_T = new TH1F("TV_T", "Thrust Value", 100, -1.5,1);
        _TA_T = new TH1F("TA_T", "Cosine of Thrust Angle", 100, 0,1);
        _TVTA_T = new TH2F("TVTA_T", "Thrust Value v Thrust Angle",100, 0,1,100,0,1);  
    }
    if(_thrustDetectability==1){

        _rootfile = new TFile("ThrustSusy_.39133._DAB.root","RECREATE");

        _TV_DAB = new TH1F("TV_DAB", "Thrust Value", 100, -1.5,1);
        _TA_DAB = new TH1F("TA_DAB", "Cosine of Thrust Angle", 100, 0,1);
        _TVTA_DAB = new TH2F("TVTA_DAB", "Thrust Value v Thrust Angle",100, 0,1,100,0,1);  

    }
    if(_thrustDetectability==2){
        _rootfile = new TFile("ThrustSusy_.39133._DED.root","RECREATE");

        _TV_DED = new TH1F("TV_DED", "Thrust Value", 100, -1.5,1);
        _TA_DED = new TH1F("TA_DED", "Cosine of Thrust Angle", 100, 0,1);
        _TVTA_DED = new TH2F("TVTA_DED", "Thrust Value v Thrust Angle",100, 0,1,100,0,1);  

    }

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



void ThrustSusy::processRunHeader( LCRunHeader* run) { 
    //run->parameters().setValue("thrust",12300321);
    //    _nRun++ ;
    cout << "Entered processRunHeader" << endl; 
} 



void ThrustSusy::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...
    cout << "EVENT: "<< _nEvt << endl; 
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
    _nEvt ++ ; // different from original-moved out of for loop - summer 
    //reset variables for output   
    _principleThrustSusyValue = -1;
    _majorThrustSusyValue     = -1;
    _minorThrustSusyValue     = -1;
    _principleThrustSusyAxis.set(0,0,0);
    _majorThrustSusyAxis.set(0,0,0);
    _minorThrustSusyAxis.set(0,0,0);

    // Switch to the desired type of thrust finder
    if (_typeOfThrustSusyFinder == 1)
    {
        cout << "type of Thrust Razor Finder = 1 : Tasso Thrust Razor " << endl; 
        TassoThrustSusy();
        cout << "type of Thrust Razor Finder = 1 : Tasso Thrust Razor" << endl;
    }
    else if (_partMom.size()<=1)
    {
        partMomCheck += " ";
        partMomCheck += std::to_string(_nEvt);
        partMomCheck += " ";  
        TassoThrustSusy();
    }
    else if (_typeOfThrustSusyFinder == 2)
    {
        cout << "type of Thrust Razor Finder = 2 : Jetset Thrust Razor" << endl; 
        JetsetThrustSusy();
    }
    // ###write
    //    evt->parameters().setValue("thrust",_principleThrustSusyValue);

    FloatVec thrax;
    thrax.clear();
    thrax.push_back(_principleThrustSusyAxis.x());
    thrax.push_back(_principleThrustSusyAxis.y());
    thrax.push_back(_principleThrustSusyAxis.z());

    _inParVec->parameters().setValue("principleThrustSusyValue",_principleThrustSusyValue);
    _inParVec->parameters().setValues("principleThrustSusyAxis",thrax);

    if (_typeOfThrustSusyFinder == 2)
    {
        thrax.clear();
        thrax.push_back(_majorThrustSusyAxis.x());
        thrax.push_back(_majorThrustSusyAxis.y());
        thrax.push_back(_majorThrustSusyAxis.z());

        _inParVec->parameters().setValue("majorThrustSusyValue",_majorThrustSusyValue);
        _inParVec->parameters().setValues("majorThrustSusyAxis",thrax);

        thrax.clear();
        thrax.push_back(_minorThrustSusyAxis.x());
        thrax.push_back(_minorThrustSusyAxis.y());
        thrax.push_back(_minorThrustSusyAxis.z());

        _inParVec->parameters().setValue("minorThrustSusyValue",_minorThrustSusyValue);
        _inParVec->parameters().setValues("minorThrustSusyAxis",thrax);

        float Oblateness;
        Oblateness = _majorThrustSusyValue - _minorThrustSusyValue;
        _inParVec->parameters().setValue("Oblateness",Oblateness);
        if ( (_majorThrustSusyValue < 0) || (_minorThrustSusyValue < 0) )
        {
            _inParVec->parameters().setValue("Oblateness",-1);
        }
    }


    // these are the final values:
    streamlog_out( DEBUG4 ) << " thrust: " << _principleThrustSusyValue << " TV: " << _principleThrustSusyAxis << endl;
    streamlog_out( DEBUG4 ) << "  major: " << _majorThrustSusyValue << " TV: " << _majorThrustSusyAxis << endl;
    streamlog_out( DEBUG4 ) << "  minor: " << _minorThrustSusyValue << " TV: " << _minorThrustSusyAxis << endl;
    cout << " thrust: " << _principleThrustSusyValue << " TV: " << _principleThrustSusyAxis << endl;
    cout <<"                       "<< _principleThrustSusyAxis.x()<<","<< _principleThrustSusyAxis.y()<< ","<<_principleThrustSusyAxis.z()<<endl;
    cout << "  major: " << _majorThrustSusyValue << " TV: " << _majorThrustSusyAxis << endl;
    cout << "  minor: " << _minorThrustSusyValue << " TV: " << _minorThrustSusyAxis << endl;

    double ptaX = _principleThrustSusyAxis.x();
    double ptaY = _principleThrustSusyAxis.y();
    double ptaZ = _principleThrustSusyAxis.z();

    double cosTT = abs(ptaZ); // cosine of the theta angle of thrust axis
    cout << "cos   "<< ptaZ<< endl;
    cout << "|cos| "<<abs(ptaZ) << endl; 

    if(_thrustDetectability == 0){
        _TV_T->Fill(_principleThrustSusyValue);
        _TA_T->Fill(cosTT);
        _TVTA_T->Fill(cosTT,_principleThrustSusyValue);
    }
    if(_thrustDetectability == 1){

        _TV_DAB->Fill(_principleThrustSusyValue);
        _TA_DAB->Fill(cosTT);
        _TVTA_DAB->Fill(cosTT,_principleThrustSusyValue);
    }
    if(_thrustDetectability == 2){

        _TV_DED->Fill(_principleThrustSusyValue);
        _TA_DED->Fill(cosTT);
        _TVTA_DED->Fill(cosTT,_principleThrustSusyValue);
    }    

    if (_principleThrustSusyValue >= _max) _max = _principleThrustSusyValue;
    if (_principleThrustSusyValue <= _min) _min = _principleThrustSusyValue;
    cout << "End EVENT "<< _nEvt<< endl;
}




void ThrustSusy::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void ThrustSusy::end(){ 
    _rootfile->Write();
    cout << partMomCheck << endl; 
}
int ThrustSusy::TassoThrustSusy(){
    int ThrustError = 0;
    if(_partMom.size()<=0)
    {
        //partMom0 = true;
        ThrustError = -1;
        //cout<< "NO PARTICLES IN EVENT AAAAAH" << endl;
        _principleThrustSusyValue = -1;
        _principleThrustSusyAxis.set(0,0,0);
    }
    if(_partMom.size()==1)
    {
        cout << "partMom.size"<< _partMom.size() << endl;
        //partMom1 = true;
        _principleThrustSusyValue = 1;
        _principleThrustSusyAxis = _partMom[0];
    }
    return ThrustError;
}

/*
   int ThrustSusy::TassoThrustSusy(){
   int ThrustError = 0;
   Hep3Vector tvec;

// No particle in Event: Error
if (_inParVec->getNumberOfElements()<=0)
{
ThrustError = -1;
_principleThrustSusyValue = 0;
_principleThrustSusyAxis.set(0,0,0);
}
// only one Particle in Event: Thrust direction = direction of particle
else if (_inParVec->getNumberOfElements()==1)
{
_principleThrustSusyValue = 1;
_principleThrustSusyAxis = _partMom[0];
}
else
{
Hep3Vector ptm, ptot, pt;
std::vector<Hep3Vector> pc;
float sp,u, pp, tmax, t;

sp = 0;
for (int i=0;i < _inParVec->getNumberOfElements();i++)
{
pp = _partMom[i].mag();
sp += pp;
ptot += _partMom[i];
} // for i
// ###
for (int m = 0; m <= 2; m++ )
ptot[m] *= 0.5;
tmax = 0;
for (int k = 1; k < _inParVec->getNumberOfElements(); k++)
{
for (int j = 0; j <= k-1;j++)
{
// cross product
tvec = _partMom[j].cross(_partMom[k]);
pt = -1 * ptot;
for (int l = 0; l < _inParVec->getNumberOfElements(); l++)
{
if (l==k) continue;
if (l==j) continue;
u = _partMom[l] * tvec;
if (u<0) continue;
pt += _partMom[l];
} // for l

while(!pc.empty())
{
pc.pop_back();
}
// note: the order is important!!!
pc.push_back(pt);
pc.push_back(pt + _partMom[k]);
pc.push_back(pt + _partMom[j]);
pc.push_back(pc[2] + _partMom[k]);
for (int m = 0; m <= 3; m++ )
{
t = pc[m].mag2();
if (t <= tmax) continue;
tmax = t;
ptm = pc[m];
} // for m
} // for j
} // for k
_principleThrustSusyValue = 2 * sqrt(tmax) / sp;
tvec = ptm;
} // end else

// Normalization of thrust vector
double ax = 0;
ax = tvec.mag();
if (ax != 0)
{
    ax = 1/ax;
    _principleThrustSusyAxis = ax * tvec;
}
else
{
    ThrustError = -1;
    _principleThrustSusyValue = -1;
    _principleThrustSusyAxis.set(0,0,0);
}
return ThrustError;
}
*/
int ThrustSusy::JetsetThrustSusy(){
    const int nwork=11,iFastMax = 4,iGood=2;
    const float dConv=0.0001; // 0.0001
    int sgn;
    double theta=0,phi=0;
    double thp,thps,tds,tmax,dOblateness;
    vector<Hep3Vector> TAxes(3),Fast(iFastMax+1),Workv(nwork);
    vector<double> Workf(nwork),dThrustSusy(3);
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
        dThrustSusy[pass] = 0;
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
            if ( thp < dThrustSusy[pass] - dConv ) continue;
            if ( thp > dThrustSusy[pass] + dConv )
            {
                nagree = 0;
                //          if (myrnd.flat() > 0.49999)
                //        {sgn = 1;} else {sgn=-1;}
                sgn = 1;
                TAxes[pass] = sgn*tpr/(tmax*thp);
                dThrustSusy[pass] = thp;
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
    dThrustSusy[2] = thp/tmax;

    // Rotate back to original coordinate system.
    for ( unsigned int i = 0;i < TAxes.size(); i++)
    {
        TAxes[i].rotateY(theta);
        TAxes[i].rotateZ(phi);
    }
    dOblateness = dThrustSusy[1] - dThrustSusy[2];

    _principleThrustSusyValue = dThrustSusy[0];
    _majorThrustSusyValue     = dThrustSusy[1];
    _minorThrustSusyValue     = dThrustSusy[2];
    _principleThrustSusyAxis  =   TAxes[0];
    _majorThrustSusyAxis      =   TAxes[1];
    _minorThrustSusyAxis      =   TAxes[2];


    return  0;
}

//______________________________________________________________
// helper function to get sign of b
double ThrustSusy::sign(double a, double b)
{
    if ( b < 0 )
    { return -fabs(a); } else { return fabs(a); }
}
//______________________________________________________________
double ThrustSusy::min(double a, double b)
{
    if ( a < b )
    { return a; } else { return b; }
}
//______________________________________________________________


