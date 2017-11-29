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
 */

#include "SV_Overlay.h"
#include "scipp_ilc_utilities.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/MCParticle.h>

#include <TFile.h>
#include <TH2D.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"



using namespace lcio;
using namespace marlin;
using namespace std;
using namespace scipp_ilc;

SV_Overlay SV_Overlay;

static TFile* _rootfile;
static TH1F* _V_Dtd;
static TH1F* _V_Dbl;
static TH1F* _V_Tru;
static TH1F* _S_Dtd;
static TH1F* _S_Dbl;
static TH1F* _S_Tru;
static TH1F* _M_Dtd;
static TH1F* _M_Dbl;
static TH1F* _M_Tru;
static bool v = false;


// below is a prototype. This tells the program there is a print statement that we want to implement through the code even though the if statement is all the way to the end
void print (string input="") ;

SV_Overlay::SV_Overlay() : Processor("SV_Overlay") {
  // modify processor description
  _description = "Protype Processor" ;

  // register steering parameters: name, description, class-variable, default value
  registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
  registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}



void SV_Overlay::init() { 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  _rootfile = new TFile("SV_Overlay.root","RECREATE");
    _V_Dtd = new TH1F("V_Dtd","Detected Vector",100,0,150);
    _V_Dbl = new TH1F("V_Dbl","Detectable Vector",100,0,200);
    _V_Tru = new TH1F("V_Tru","True Vector",100,0,150);

    _S_Dtd = new TH1F("S_Dtd","Detected Scalar",100,0,200);
    _S_Dbl = new TH1F("S_Dbl","Detectable Scalar",100,0,200);
    _S_Tru = new TH1F("S_Tru","True Scalar",100,0,200);
 
  
    _M_Dtd = new TH1F("M_Dtd","Detected Mass",100,0,150);
    _M_Dbl = new TH1F("M_Dbl","Detectable Mass",100,0,60);
    _M_Tru = new TH1F("M_Tru","True Mass",100,0,4.0);
  
    //_S = new TH1D("S", "Scalar", 200, 0.0, 23.0);
    //_V = new TH1D("V", "Vector", 1000, 0.0, 1.5); 
    //_M = new TH1D("M", "Mass", 1000, 0.0, 150.0); 


  // usually a good idea to
  //printParameters() ;
  _nEvt = 0 ;

  cout << "Anything to let me know " << endl;

}



void SV_Overlay::processRunHeader( LCRunHeader* run) { 
  //    _nRun++ ;
} 



void SV_Overlay::processEvent( LCEvent * evt ) { 
  // this gets called for every event 
  // usually the working horse ...
  LCCollection* col = evt->getCollection( _colName ) ;
  //cout << "# of Elements "<<col->getNumberOfElements()<<endl;
  
  double vec[4][3];
  double scalars[4];
  double energy[4];


  double tot_mom[]={0, 0, 0, 0};
  double tot_mag[]={0, 0, 0, 0};
  double S=0,V=0,M=0; 
  int stat;
  // this will only be entered if the collection is available
  if( col != NULL ){
    int nElements = col->getNumberOfElements()  ;
    print();
    print();
    // print("=====EVENT: " + to_string( _nEvt ) + " ===== " );
    vector<MCParticle*> system;
    //Pre analysis - I need to identify particle, status, and energy- If stat =1 (final state) highest energy of electron and positron, remove. Also energy
    
    double MaxEnerge=0.0, MaxEnergp=0.0;

    for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
      MCParticle* particle = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
      stat = particle->getGeneratorStatus();
      int id;

       try{ 
            id = particle->getPDG();
          //  cout << id << endl;  
       }
       catch(const std::exception& e){
	 //         cout << "exception caught with message " << e.what() << "\n";
       }  
      //yo dis is kewl. like it says "if not final state, skip" so ya. 
      //TWO EQUALS ARE EQUAL> ONE EQUAL IS DEFINING
       //      if(stat!=1) continue;
       //      double ID = particle->getPDG();
       //      double energy=particle->getEnergy();
     
       //      if( ID==11 && energy>MaxEnerge){
       //	  MaxEnerge=energy;
       //      }else if( ID==-11 && energy>MaxEnergp){
       //	  MaxEnergp=energy;
       // }
        if(stat==1){

                bool isDarkMatter = (id == 1000022);
                if(isDarkMatter) continue ;
                double E = particle->getEnergy();
                const double* P = particle->getMomentum();
                double Pmag = sqrt(P[0]*P[0]+P[1]*P[1]+P[2]*P[2]);
                double cos = P[2]/Pmag;
                double scalar = sqrt(P[0]*P[0]+P[1]*P[1]); 
                bool isNeutrino = (
                    id == 12 || id == -12 ||
                    id == 14 || id == -14 ||
                    id == 16 || id == -16 ||
                    id == 18 || id == -18);
                bool isForward = ( cos > 0.9 || cos < -0.9);               
		
		//cout << "anything u want "<< scalar << endl;

                scalars[0]+=scalar; //true
                vec[0][0]+=P[0];
                vec[0][1]+=P[1];
                vec[0][2]+=P[2];
                energy[0]+=E;                        
                if(!isDarkMatter && !isNeutrino){ //detectable
                    scalars[2]+=scalar;
                    vec[2][0]+=P[0];
                    vec[2][1]+=P[1];
                    vec[2][2]+=P[2];
                    energy[2]+=E;
                    if(!isForward){
                        scalars[1]+=scalar; //detected
                        vec[1][0]+=P[0];
                        vec[1][1]+=P[1];
                        vec[1][2]+=P[2];
                        energy[1]+=E;      
                    }
                }
                 
	}//end final state


 
    }//end for



    // cout << "event "<< _nEvt <<" finished " << endl;
    //all
    double total_true_scalar = scalars[0];
    double total_detected_scalar = scalars[1];
    double total_detectable_scalar = scalars[2];

    double total_true_mass_squared = energy[0]*energy[0]-
      (vec[0][0]*vec[0][0]+vec[0][1]*vec[0][1]+
       vec[0][2]*vec[0][2]);
    double total_true_mass = sqrt(total_true_mass_squared);
    double total_detected_mass_squared = energy[1]*energy[1]-
      (vec[1][0]*vec[1][0]+vec[1][1]*vec[1][1]+
       vec[1][2]*vec[1][2]);
    double total_detected_mass = sqrt(total_detected_mass_squared);
    double total_detectable_mass_squared = energy[2]*energy[2]-
      (vec[2][0]*vec[2][0]+vec[2][1]*vec[2][1]+
       vec[2][2]*vec[2][2]);
    double total_detectable_mass = sqrt(total_detectable_mass_squared);
            
    double total_true_vector = sqrt(vec[0][0]*vec[0][0]+vec[0][1]*vec[0][1]+vec[0][2]*vec[0][2]); 
    double total_detected_vector = sqrt(vec[1][0]*vec[1][0]+vec[1][1]*vec[1][1]+vec[1][2]*vec[1][2]); 
    double total_detectable_vector = sqrt(vec[2][0]*vec[2][0]+vec[2][1]*vec[2][1]+vec[2][2]*vec[2][2]);
         
    _V_Dtd->Fill(total_detected_vector);
    _V_Dbl->Fill(total_detectable_vector);
    _V_Tru->Fill(total_true_vector);
        
    _S_Dtd->Fill(total_detected_scalar);
    _S_Dbl->Fill(total_detectable_scalar);
    _S_Tru->Fill(total_true_scalar);

    _M_Dtd->Fill(total_detected_mass);
    _M_Dbl->Fill(total_detectable_mass);
    _M_Tru->Fill(total_true_mass);
  }



  //    for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
  //  MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
  //  double energy=hit->getEnergy();

  //  stat = hit->getGeneratorStatus();
  //  if(stat==1){ 
  //	double ID = hit->getPDG();
  //	if( ID==11 && energy==MaxEnerge){
	  
  //	  continue;
  //	}
  //	if( ID==-11 && energy==MaxEnergp) {
  //	  continue;
  //	}


  //	double mom[3] = {0};
  //	mom[0] = hit->getMomentum()[0];
  //	mom[1] = hit->getMomentum()[1];
  //	mom[2] = hit->getMomentum()[2];

				
  //	double mag = sqrt(pow(mom[0], 2)+pow(mom[1], 2));
  //	S+=mag;


  //	tot_mom[0]+=mom[0];
  //	tot_mom[1]+=mom[1];
  //	tot_mom[2]+=mom[2];
  //	tot_mom[3]+=hit->getEnergy();
  
  //	tot_mag[0]+=abs(mom[0]);
  //	tot_mag[1]+=abs(mom[1]);
  //	tot_mag[2]+=abs(mom[2]);
  //	tot_mag[3]+=hit->getEnergy();


	// }//end final state
	//}//end for
    //for(MCParticle* particle : system)

  //    V=sqrt(pow(tot_mom[0], 2)+pow(tot_mom[1], 2));
    //M=sqrt(pow(tot_mom[3], 2)-pow(tot_mom[0], 2) -pow(tot_mom[1], 2) -pow(tot_mom[2], 2));
  //    M=sqrt(pow(tot_mag[3], 2)-pow(tot_mag[0], 2) -pow(tot_mag[1], 2) -pow(tot_mag[2], 2));

  //  _S->Fill(S); 

  //  _V->Fill(V);
  //  _M->Fill(M);
  //   }

  // _nEvt ++ ;
  // }

}

void SV_Overlay::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void SV_Overlay::end(){ 
  _rootfile->Write();
}

void print (string input) {
  if (v) cout << input << endl;
}
