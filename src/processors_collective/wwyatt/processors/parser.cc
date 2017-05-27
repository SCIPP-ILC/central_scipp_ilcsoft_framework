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

#include "parser.h"
#include "scipp_ilc_utilities.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/MCParticle.h>

#include <TFile.h>
#include <TH2D.h>

#include <cmath>
#include <vector>
#include <map>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"


using namespace lcio;
using namespace marlin;
using namespace std;


parser parser;

static TFile* _rootfile;
static int nBhabha=0;
static int nBase=0;
static int nTwoPhoton=0;


parser::parser() : Processor("parser") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
    
    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}


void parser::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;
    _rootfile = new TFile("parser.root","RECREATE");
    _nEvt = 0 ;
}

void parser::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 

void parser::processEvent( LCEvent * evt ) { 
    LCCollection* col = evt->getCollection( _colName );
    int stat, id =0;
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
	int trees = numberOfTrees(evt, false);
	cout << "There are " << trees << " trees." << endl;
	if(nElements==4){
	  ++nBase;
	  return;
	}
	if(true){
	  ++nBhabha;
	}else{
	  ++nTwoPhoton;
	}

	
	/*	int eps = 0;
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );        
	   
            id = hit->getPDG();
            stat = hit->getGeneratorStatus();
            if(stat==1){
            }//end final state   
	    }//end for*/
    }
    _nEvt ++ ;
}


void parser::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void parser::end(){ 
  cout << "Number of empty events: " << nBase << endl;
  cout << "Number of Bhabha events: " << nBhabha << endl;
  cout << "Number of Two Photon events: " << nTwoPhoton << endl;
  _rootfile->Write();
}

int parser::numberOfTrees(LCEvent * evt, bool v){
  LCCollection* col = evt->getCollection( _colName ) ;
  int id, stat;

  //v is verbosity, if v is true, then this function will print a lot a data to the console.
  //I do not know how to send data only to the log or error console.
  if(v)cout <<endl<< "***** Event " << _nEvt << ". *****" << endl;

  if( col != NULL ){
    int nElements = col->getNumberOfElements()  ;
    vector<pair<int,double>> children_index; //used to index the children map.
    vector<pair<int,double>> parent_index; //used to index the children map.
    map<pair<int,double>, MCParticle *> parents; //Set of all particles with PDG 0
    map<pair<int,double>, MCParticle *> children; // All other particles
    //There maps have key of a pair of int as their PDG and a double as their energy.
    //The idea is that each even has a unique particle with a unique enegy.
    //If not an error will throw and I will find a way to make them unique.

    //Looping through all MCParitles
    for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
      MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex));
      id = hit->getPDG();
      stat = hit->getGeneratorStatus();      

      if(v)cout << "-- Particle index " << hitIndex << " --" << endl;

      //If the gen status is zero then it is a initial particle and add it to the parent map.
      if(stat == 0){
	pair<int,double> key(id, hit->getEnergy());
	parent_index.push_back(key);
	pair<pair<int,double>,MCParticle *> ret(key, hit);
	bool duplicate=parents.insert(ret).second;
	//Check to see if the particle is already in the map
	if(duplicate == false){
	  //Paritle is already there, if this prints then pair<id, enery> is not unique.
	  cout << "Error: there is multiple of particle PDG " << id << "." << endl;
	}
      }else if(stat == 2){
	pair<int,double> key(id, hit->getEnergy());
	children_index.push_back(key);
	pair<pair<int,double>, MCParticle *> ret(key, hit);
	bool duplicate = children.insert(ret).second;
	//Check to see if the particle is already in the map
	if(duplicate == false){
	  //Paritle is already there, if this prints then pair<id, enery> is not unique.
	  cout << "Error: there is multiple of particle PDG " << id << "." << endl;
	}
      }
      //Verbose Print statements
      if(v && stat != 3){
	cout << "This particle, " << id << ", with gen status of " << stat<< " has " << hit->getParents().size() << " parents and has the energy " << hit->getEnergy() << "."<< endl;
	for(MCParticle* parent : hit->getParents()){
	  cout << "For the particle with PGD: " << id<< " this the parent id: " << parent->getPDG() << " and energy " << parent->getEnergy() << endl;
	}	     
      }
    }

    cout << "number of parents: " << parents.size() << endl;
    cout << "number of children: " << children.size() << endl;

    int numTrees = 0;
  }
  return 1;
}
