#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0

/*
 * author Summer Zuber
 * November 23, 2016
 * This is an initial attempt to calculate and use Razor Variables with our degenerate Susy Events 
 */

#include "SusyRazorVariables.h"
#include "scipp_ilc_utilities.h"
#include <iostream>
#include <cmath>

#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/MCParticle.h>

#include <TFile.h>
#include <TH2D.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include <list>


using namespace lcio;
using namespace marlin;
using namespace std;

SusyRazorVariables SusyRazorVariables;

static TFile* _rootfile;

static TH1F* _FirstRazorPlot;
static TH1F* _RPlot;

SusyRazorVariables::SusyRazorVariables() : Processor("SusyRazorVariables") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
}



void SusyRazorVariables::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("SusyRazorVariables.root","RECREATE");
    _FirstRazorPlot = new TH1F("FirstRazorPlot","My First Razor Plot", 100,0,20); 
    _RPlot = new TH1F("RPlot", " R = MTR/MR",100,0,10);
    // usually a good idea to
    //printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;
}



void SusyRazorVariables::processRunHeader( LCRunHeader* run) { 
    //    _nRun++ ;
} 



void SusyRazorVariables::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    LCCollection* col = evt->getCollection( _colName ) ;
    cout << endl;
    cout << endl;
    cout << endl;
    cout << "event = " << _nEvt << endl;

    //these will be the two tau's for the event:
    MCParticle* tau1;
    MCParticle* tau2;
    
    // these will be my momentum and energy sums for the event: 
    double vec[4][3];  // momentum 3 vectors of: tau 1, tau 2, lsp 1, lsp 2
    double scalarPT[4]; // magnitude of the PT of the same categories 
    double energy[4];  // energy of same categories 
    double vecFinals[1][3]; // the momentum 3 vectors of final state particles: for not just all final state particles

    //particle identifiers
    int id, stat; 

    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;

        int i = 0; // this is to count the 2 tau leptons 
        // For each particle in Event ...
        for(int particleIndex = 0; particleIndex < nElements ; particleIndex++){
            MCParticle* particle = dynamic_cast<MCParticle*>( col->getElementAt(particleIndex) ); 
            try{ 
                id = particle->getPDG(); 
                stat = particle->getGeneratorStatus();
            }
            catch(const std::exception& e){
                cout << "exception caught with message " << e.what() << "\n";
            } 
            if (id==15 || id == -15){
                //cout << particle << " " << particle->getPDG() << endl;
                for(MCParticle* parent : particle->getParents()){
                    //cout << "parent: parent, id" << parent << " " << parent->getPDG() << endl;

                    if(parent->getPDG() == 1000015 || parent->getPDG() == -1000015){
                        
                        if(i==0){tau1 = particle;}
                        if(i==1){tau2 = particle;}
                        i++;
                    }
                }

            } 
        }//end for
        cout << "tau 1 " << tau1 << " " << tau1->getPDG() << endl;
        cout << "tau 2 " << tau2 << " " << tau2->getPDG() << endl;
        
        vec[0][0]+=tau1->getMomentum()[0];
        vec[0][1]+=tau1->getMomentum()[1];
        vec[0][2]+=tau1->getMomentum()[2];
        vec[1][0]+=tau2->getMomentum()[0];
        vec[1][1]+=tau2->getMomentum()[1];
        vec[1][2]+=tau2->getMomentum()[2];

        energy[0]+= tau1->getEnergy(); 
        energy[1]+= tau2->getEnergy();

        scalarPT[0]+=sqrt(vec[0][0]*vec[0][0]+vec[0][1]*vec[0][1]); //tau 1
        scalarPT[1]+=sqrt(vec[1][0]+vec[1][0]+vec[1][1]*vec[1][1]); //tau 2

        //transform to R frame 
        double beta = (tau1->getEnergy() - tau2->getEnergy())/(vec[0][2] - vec[1][2]);
        cout << "BETA :"<< beta<<endl;
        
        double beta2 = pow(beta,2);
        double gamma = 1/(sqrt(1-beta2));
        // want to check whether this transformation equilizes 3-momentum of tau 1 and tau 2 
        // I don't need the energy, here are transformed three vectors
        double vecT1[3] = {vec[0][0], vec[0][1], - gamma*beta*energy[0]+gamma*vec[0][2]}; 
        double vecT2[3] = {vec[1][0], vec[1][1], - gamma*beta*energy[1]+gamma*vec[1][2]};
        double magT1 = sqrt(vecT1[0]*vecT1[0]+vecT1[1]*vecT1[1]+vecT1[2]*vecT1[2]);
        double magT2 = sqrt(vecT2[0]*vecT2[0]+vecT2[1]*vecT2[1]+vecT2[2]*vecT2[2]);
        double dif3mag = magT1 - magT2;
        cout << dif3mag << endl;
        
        for(int particleIndex = 0; particleIndex < nElements ; particleIndex++){
            MCParticle* particle = dynamic_cast<MCParticle*>( col->getElementAt(particleIndex) );
                                
            try{
                id = particle->getPDG();
                stat = particle->getGeneratorStatus();
            }
            catch(const std::exception& e){
                cout << "exception caught with message "<< e.what() << "\n";
            }
            
            double part4Vector[4] = {particle->getEnergy(), particle->getMomentum()[0], 
                                         particle->getMomentum()[1], particle->getMomentum()[2]};
            double R4Vector[4] = {gamma*part4Vector[0]-gamma*beta*part4Vector[3],part4Vector[1], part4Vector[2], 
                                    -gamma*beta*part4Vector[0]+gamma*part4Vector[3]};
            //double *R4Vector[4] = {Transform2RFrame( particle4Vector, beta )};
            //cout << "Four Vec    " << endl;
            //cout << particle4Vector[0] <<" "<<particle4Vector[1]<<" "<<particle4Vector[2]<<" "<<particle4Vector[3]<< endl;
            //cout << "Transformed " << R4Vector[0]        <<" "<<R4Vector[1]       <<" "<<R4Vector[2]       <<" "<<R4Vector[3]<<" "<<R4Vector[4]<< endl;
            if(stat==1){
                vecFinals[0][0]+=R4Vector[1];
                vecFinals[0][1]+=R4Vector[2];
                vecFinals[0][2]+=R4Vector[3];
            }
        } //end for

        // create MTR variable:
        // missing transverse energy (i.e. missing transverse momentum)
        double ETM[2] = { -vecFinals[0][0], -vecFinals[0][1]}; 
        double ETMmag = sqrt(ETM[0]*ETM[0]+ETM[1]*ETM[1]);
        
        // pt j1 + pt j2 vector :
        double pt_taus[2] = {vec[0][0]+vec[1][0], vec[0][1]+vec[1][1]};
        // dot product of ETM and pt_taus:
        double dotProd = ETM[0]*pt_taus[0]+ETM[1]*pt_taus[1]; 
        double MTR = sqrt((ETMmag*(scalarPT[0]+scalarPT[1])-dotProd)/2);
        double R = MTR/(2*magT1);
        if(beta2<=1){
            cout << "filling  " << MTR << endl;
            _FirstRazorPlot->Fill(MTR);
            _RPlot->Fill(R);
        }
    }//ind if col

    _nEvt ++;

    cout << "event "<< _nEvt <<" finished " << endl;
}//end process


// function to transform into R frame (unnecessary) 
/* 
double *SusyRazorVariables::Transform2RFrame(double in[4], double beta){
    cout << "----------------------------"<<endl;
    cout << "Running Transform Function!"<<endl;
    double beta2 = pow(beta,2);
    cout << "BETA SQUARED: "<<beta2<<endl; 
    double gamma = 1/(sqrt(1-pow(beta,2)));
    cout << "GAMMA: "<<gamma<<endl;
    cout << "In 4 Vec: " <<endl;
    cout << in[0] <<" "<<in[1]<<" "<<in[2]<<" "<<in[3]<<endl;
    double out[4] = {gamma*in[0]-gamma*beta*in[3], in[1], in[2], -gamma*beta*in[0]+gamma*in[3]}; 
    cout << "Out:"<<endl;
    cout << out[0]<<" "<<out[1]<<" "<<out[2]<<" "<<out[3]<<endl;
    cout << "----------------------------"<<endl; 
    return out;     
    
}
*/

void SusyRazorVariables::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void SusyRazorVariables::end(){ 

    _rootfile->Write();
}
