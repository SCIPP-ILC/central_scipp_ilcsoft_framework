#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0
/* 
 * author Summer Zuber 
 * January 17, 2018 
 *
 * This is a copy of ThrustRazor (for susy files) which I am going to edit
 * so that instead of using the thrust to define the jets it uses a jet algorithm. This
 * way the jets are properly defined for use in the Razor Variable calculation.
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

//ClassImp(JetFinder)
using namespace lcio;
using namespace marlin;
using namespace std;

const Int_t JetRazor::DURHAM = 1; // jetFinder
const Int_t JetRazor::JADE   = 2; // jetFinder
const Int_t JetRazor::JADEE  = 3; // jetFinder  

const Int_t JetRazor::UNASSOC = -999; //jetFinder 

/*JetRazor::JetRazor(Double_t ycut): //changed from jetFinder
  m_evis(0.),m_dycut(ycut),
  m_algorithm(JetFinder::DURHAM){ // Default constructor  uses Durham
  m_4vec = new TObjArray();
  }*/
static TFile* _rootfile;

static TH1F* _R_T;
static TH1F* _R_DAB;
static TH1F* _R_DED;
static TH1F* _MR_T;
static TH1F* _MR_DAB;
static TH1F* _MR_DED; 

JetRazor JetRazor;
//JetFinder JetFinder;

JetRazor::JetRazor() : Processor("JetRazor") {
    cout << "ummm what is happening " << endl; 
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
    registerProcessorParameter( "typeofJetFinder", 
            "Type of Jet Finding Algorithm to be used:\n#\t1 : DURHAM \n#t2 : JADE \n#t3 : JADEE", 
            _typeofJetFinder, 1 ); // added this, needs finish, need? 
    registerProcessorParameter( "jetDetectability",
            "Detectability level used for jet algorithm:\n#\t0 : True \n#t1 : Detectable \n#t2 : Detected" ,
            _jetDetectability, 2 );
}

void JetRazor::init() { 
    streamlog_out(DEBUG)  << "   init called  " <<  std::endl ;
    cout << "initialized" << endl;
    if(_jetDetectability==0){_rootfile = new TFile("JetRazor_.39133._T.root","RECREATE");
        _R_T = new TH1F("R_T", "R =MTR/MR",130,-3,10);
        _MR_T = new TH1F("MR_T","MR", 100, 0 ,10); 
    }
    if(_jetDetectability==1){_rootfile = new TFile("JetRazor_.39133._DAB.root","RECREATE");
        _MR_DAB = new TH1F("MR_DAB","MR", 100, 0 ,10); 
        _R_DAB = new TH1F("R_DAB", "R =MTR/MR",130,-3,10);
    }
    if(_jetDetectability==2){_rootfile = new TFile("JetRazor_.39133._DED.root","RECREATE");
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
    parp1 = false;
    parp0 = false; 
    cout << "-------------------------------------------------------------------------" << endl;  
    cout << "EVENT: " << _nEvt << endl; 
    _inParVec = evt->getCollection( _colName) ;
    cout << "num of elements " << _inParVec->getNumberOfElements() << endl;
    if (!_parp.empty()) _parp.clear(); // _parp global Hep3Vector

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

        const double* parp = aPart->getMomentum(); // particle momentum vector

        if(stat==1){
            cout << "id: " << id<< endl;
            

            bool isDarkMatter = (id == 1000022);
            bool isNeutrino = (
                id == 12 || id == -12 ||
                id == 14 || id == -14 ||
                id == 16 || id == -16 ||
                id == 18 || id == -18);

            double parpMag = sqrt(parp[0]*parp[0]+parp[1]*parp[1]+parp[2]*parp[2]);
            double cos = parp[2]/parpMag;
            bool isForward = ( cos > 0.9 || cos < - 0.9);
            bool isDetectable = (!isDarkMatter && !isNeutrino);
            bool isDetected = (isDetectable &&  !isForward  );
            
            if(_jetDetectability == 0){
                if(!isDarkMatter){
                    _parp.push_back( Hep3Vector(parp[0], parp[1], parp[2]) );
                }
            }
            if(_jetDetectability == 1){
                if(isDetectable){
                    _parp.push_back( Hep3Vector(parp[0], parp[1], parp[2]) );
                }
            }

            if(_jetDetectability == 2){ 
                if(isDetected){ 
                    _parp.push_back( Hep3Vector(parp[0], parp[1], parp[2]) ); 
                }
            }
        } // stat = 1
    } // for particle 
    cout << "end loop #1"<<endl; 
    //JetRazor.setPartList(_parp);
    
    double vec[2][3][4]; // jet 1, jet 2 : true detectable, detected : energy, momx, momy, momz
    double Rvec[3][4]; // true, detectable, detected : energy, px, py, pz 
    int d = _jetDetectability;
    double R;
    double beta2;

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

    double beta = (vec[0][d][0]-vec[1][d][0])/(vec[0][d][3]-vec[1][d][3]); // beta using right particles  
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

    cout << "End EVENT "<< _nEvt<< endl;
    _nEvt ++ ;  
}


void JetRazor::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void JetRazor::end(){ 
    _rootfile->Write();
    cout << parpCheck << endl;
    cout << betaCheck << endl; 
}
/*struct parDetectability
{
    bool isDarkMatter;
    bool isTrue;
    bool isNeutrino;
    bool isDetectable;
    bool isForward;
    bool isDetected;
}
parDetectability JetRazor::findParDetectability(int id, const double* parp){

    bool isDarkMatter = (id == 1000022);
    bool isNeutrino = (
            id == 12 || id == -12 ||
            id == 14 || id == -14 ||
            id == 16 || id == -16 ||
            id == 18 || id == -18);

    double parpMag = sqrt(parp[0]*parp[0]+parp[1]*parp[1]+parp[2]*parp[2]);
    double cos = parp[2]/parpMag;
    bool isForward = ( cos > 0.9 || cos < - 0.9);
    bool isDetectable = (!isDarkMatter && !isNeutrino);
    bool isDetected = (isDetectable &&  !isForward  );
    struct parDetectability parD;
    parD.isDarkMatter=isDarkMatter;
    parD.isTrue=!isDarkMatter;
    parD.isNeutrino = isNeutrino;
    parD.isDetectable = isDetectable;
    parD.isForward = isForward;
    parD.isDetected = isDetected; 
    return parD; 
}*/ //finish this later 
// helper function to get sign of b
double JetRazor::sign(double a, double b)
{
    if ( b < 0 )
    { return -fabs(a); } else { return fabs(a); }
}
//______________________________________________________________
double JetRazor::min(double a, double b)
{
    if ( a < b )
    { return a; } else { return b; }
}
//______________________________________________________________


void JetRazor::doFindJets(){ // from JetFinder

    Int_t np = m_4vec->GetEntries();  
    if (np<2) return;

    TObjArray* part = new TObjArray();
    for (Int_t Ipart=0; Ipart <np; Ipart++) {
        TVector3 vec3(((TLorentzVector*)m_4vec->At(Ipart))->X(),
                ((TLorentzVector*)m_4vec->At(Ipart))->Y(),
                ((TLorentzVector*)m_4vec->At(Ipart))->Z());
        TLorentzVector* vec4 = 
            new TLorentzVector(vec3,
                    ((TLorentzVector *)(m_4vec->At(Ipart)))->T());
        part->Add(vec4);
    }

    m_ipart_jet_assoc.Reset();
    m_ipart_jet_assoc.Set(np);
    for (Int_t m=0; m<np; m++) m_ipart_jet_assoc[m] = UNASSOC;

    m_njets = 0;


    // create invariant mass pair array.
    //
    TMatrix ymass = TMatrix(np,np);

    for (Int_t i1 = 0; i1 < np - 1; i1++ ) {
        for (Int_t i2 = i1 + 1 ; i2 < np ; i2++ ) {
            TLorentzVector &jeti1 = *(TLorentzVector*)part->At(i1);
            TLorentzVector &jeti2 = *(TLorentzVector*)part->At(i2);

            ymass(i1,i2) = calcinvmass(jeti1, jeti2);

        }
    }

    Double_t masscut = m_dycut * m_evis * m_evis ;

    for (;;)
    {
        Int_t im = -1;
        Int_t jm = -1;
        Double_t minmass = 100000000.;
        if ( minmass <= masscut ) minmass = masscut * 10.;
        //
        // find least invariant mass pair.
        //
        for(Int_t i = 0 ; i < np - 1 ; i++ ) {
            for(Int_t j = i + 1 ; j < np ; j++ ) {
                if (m_ipart_jet_assoc[i] != JetRazor::UNASSOC) continue; // jetFinder
                if (m_ipart_jet_assoc[j] != JetRazor::UNASSOC) continue; // jetFinder
                if (ymass(i,j) > minmass) continue;

                minmass = ymass(i,j);
                im = i;  jm = j;
            }
        }

        if (minmass > masscut) break;

        // combine particles im and jm.
        //
        *(TLorentzVector*)part->At(im) += *(TLorentzVector*)part->At(jm);

        for(Int_t ipart = 0; ipart < np; ipart++ ){
            if(m_ipart_jet_assoc[ipart] == jm) m_ipart_jet_assoc[ipart] = im;
        }
        //
        // Recalculate invariant masses for newly combined particle
        //
        m_ipart_jet_assoc[jm] = im;
        for (Int_t ipar = 0; ipar < np ; ipar++) {
            if (ipar == im) continue;
            if (m_ipart_jet_assoc[ipar] != UNASSOC ) continue;

            Int_t imin = TMath::Min(ipar,im);
            Int_t imax = TMath::Max(ipar,im);

            TLorentzVector &jetimin = *(TLorentzVector*)part->At(imin);
            TLorentzVector &jetimax = *(TLorentzVector*)part->At(imax);

            ymass(imin,imax) = calcinvmass(jetimin, jetimax);

        }

    }

    // finish up by filling jet array.


    for(Int_t ip = 0 ; ip < np ; ip++) {
        if (m_ipart_jet_assoc[ip] == UNASSOC) m_njets++;            
    }

    m_jet.Delete();
    m_inparts_per_jet.Reset();
    m_inparts_per_jet.Set(m_njets);  

    Int_t nj = 0;
    Int_t npart;
    m_ifewest_parts = 100000; // Starting min value
    for(Int_t i = 0 ; i < np ; i++ ){
        if (m_ipart_jet_assoc[i] != UNASSOC) continue;

        TVector3 vec3(((TLorentzVector*)part->At(i))->X(),
                ((TLorentzVector*)part->At(i))->Y(),
                ((TLorentzVector*)part->At(i))->Z());     
        TLorentzVector* jet = 
            new TLorentzVector(vec3,((TLorentzVector*)part->At(i))->T());
        m_jet.Add(jet);

        npart = 1;
        for (Int_t j = 0 ; j < np ; j++) {
            if(m_ipart_jet_assoc[j] == i) {
                m_ipart_jet_assoc[j] = nj;
                npart++;
            }
        }
        m_ipart_jet_assoc[i] = nj;
        m_inparts_per_jet[nj] = npart;
        if( npart < m_ifewest_parts) m_ifewest_parts = npart;
        nj++;
    }  
    part->Delete(); delete part;  
}
//
JetRazor::~JetRazor() { //jetFinder
    m_jet.Delete();
    m_4vec->Delete(); delete m_4vec;
}
//______________________________________________________

TLorentzVector* JetRazor::jet4vec(Int_t index) { // jetFinder
    return (TLorentzVector*)m_jet.At(index);
}
//_______________________________________________________

Int_t JetRazor::nParticlesPerJet(Int_t index) { //jetFunder
    return m_inparts_per_jet[index];
}
//_______________________________________________________

void JetRazor::setYCut(Double_t ycut) { //jetFinder
    m_dycut = ycut;
}
//_______________________________________________________

// Input the particle 4(3)-vector list
// e: 4-vector  TLorentzVector ..(px,py,pz,E) or
//    3-vector  TVector3       ..(px,py,pz) 
// If you input TVector3, the energy of particle
// will be E = sqrt(px**2 + py**2 + pz**2)
void JetRazor::setPartList(TObjArray* e) { //jetFinder 

    m_evis = 0;
    m_4vec->Delete();

    Int_t ne = e->GetEntries();
    for(Int_t i=0;i<ne;i++) {
        TObject* o = e->At(i);
        TString nam(o->IsA()->GetName());
        if (nam.Contains("TLorentzVector")) {
            TVector3 vec3(((TLorentzVector*) o)->X(),
                    ((TLorentzVector*) o)->Y(),
                    ((TLorentzVector*) o)->Z());    
            TLorentzVector* in = 
                new TLorentzVector(vec3,((TLorentzVector *) o)->T());
            m_evis += in->T();
            m_4vec->Add(in);
        }    
        else if (nam.Contains("TVector3")) {
            TVector3 vec3(((TVector3 *) o)->X(),
                    ((TVector3 *) o)->Y(),
                    ((TVector3 *) o)->Z());
            TLorentzVector* in = 
                new TLorentzVector(vec3,((TVector3 *) o)->Mag());
            m_evis += in->T();
            m_4vec->Add(in);
        }
        else {
            printf("JetFinder::setEvent input is not a TVector3 or a TLorentzVector\n");
        }
    }
}

void JetRazor::setDURHAM(){ // jetFinder 
    m_algorithm = JetRazor::DURHAM; //jetFinder
}
void JetRazor::setJADE(){ // jetFinder
    m_algorithm = JetRazor::JADE; //jetFinder 
}
void JetRazor::setJADEE(){ // jetFinder 
    m_algorithm = JetRazor::JADEE; // jetFinder 
}

Double_t JetRazor::calcinvmass(const TLorentzVector& jet1, // jetFinder 
        const TLorentzVector& jet2){
    TVector3 P_jet1 = jet1.Vect();
    TVector3 P_jet2 = jet2.Vect();
    Double_t costh = (P_jet1 * P_jet2)/P_jet1.Mag()/P_jet2.Mag();

    if     (m_algorithm == JetRazor::DURHAM) {  // DURHAM jetFinder
        Double_t minT = TMath::Min(jet1.E(),jet2.E());
        return 2. * minT*minT * (1.-costh);
    }  
    else if (m_algorithm == JetRazor::JADE)   {  // JADE jetFinder
        return 2. * (jet1.E())*(jet2.E()) * (1.-costh) ; 
    }  
    else if (m_algorithm == JetRazor::JADEE)  {  // JADE E jetFinder
        return (jet1 + jet2).M2();
    } 
    printf(" Strange Algorithm!!!! \n");
    return 0.;     
}
