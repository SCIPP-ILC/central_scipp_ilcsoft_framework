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
#include <algorithm>

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
static int nCombo=0;


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

//Returns true if it finds final state electrons or positrons.
bool parser::isBhabha(LCCollection* col, vector<vector<MCParticle*>>* trees){
  //It says there are 30 items in trees but only loops through one...
  //Probably something wrong with this pointer syntax.
  for(const auto tree: *trees){
    cout << "size: " << tree.size() << endl;    
    for(const auto hit: tree){
      cout << "id: " << hit->getPDG() << endl;
      if(hit->getGeneratorStatus() == 1){
	int id = hit->getPDG();
	if (id != 11 && id != -11 && id != 22){
	  return false;
	}
      }
    }
  }
  return true;
}

void parser::processEvent( LCEvent * evt ) { 
    LCCollection* col = evt->getCollection( _colName );
    int stat, id =0;
    bool v=true;
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
	if(v) cout << "##### New Event (  "<< nElements<< " elements )#####" << endl;
	vector<vector<MCParticle*>>* trees = nTrees(evt);
	if(v) cout << "There are " << trees->size() << " trees." << endl;
	if(trees->size() == 1){
	  if (isBhabha(col, trees)) ++nBhabha;
	  else ++nBase;
	}else if(isBhabha(col, trees)) ++nCombo;
	else ++nTwoPhoton;
    }
}


void parser::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void parser::end(){ 
  cout << "Number of empty events: " << nBase << endl;
  cout << "Number of Bhabha events: " << nBhabha << endl;
  cout << "Number of Two Photon events: " << nTwoPhoton << endl;
  cout << "Number of Two Photon & Bhabha events: " << nCombo << endl;
  _rootfile->Write();
}

void parser::addToTree(MCParticle* obj, MCParticle* associate, vector<vector<MCParticle*>>* trees, vector<MCParticle*> &all){
    if(find(all.begin(), all.end(),obj) == all.end() && find(all.begin(), all.end(),associate) == all.end()){
    //Not in tree
    all.push_back(obj);
    all.push_back(associate);
    vector<MCParticle*>* arr = new vector<MCParticle*>;
    arr->push_back(obj);
    arr->push_back(associate);
    trees->push_back(*arr);
  }else{
    for(auto tree: *trees){
      auto i = find(tree.begin(), tree.end(), obj);
      auto j = find(tree.begin(), tree.end(), associate);
      if( i != tree.end() && j != tree.end()){
	return;
      }else if( i == tree.end() ){
	all.push_back(obj);
	tree.push_back(obj);
      }else{
	all.push_back(associate);
	tree.push_back(associate);
      }
    }
  }
}


vector<vector<MCParticle*>>* parser::nTrees(LCEvent *evt, bool v){
  vector<vector<MCParticle *>>* trees = new vector<vector<MCParticle *>>;
  int numTrees = 0;
  LCCollection* col = evt->getCollection( _colName ) ;
  int id, stat;
  if( col != NULL ){
    int nElements = col->getNumberOfElements();
    vector<MCParticle *> trees_all;
    for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
      MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex));
      vector<MCParticle*> parents = hit->getParents();
      vector<MCParticle*> daughters = hit->getDaughters();
      if (hit->getGeneratorStatus()==3) continue;
      for(auto  parent: parents){	
	addToTree(hit, parent, trees, trees_all);
      }
      for(auto  child: daughters){
	addToTree(hit, child, trees, trees_all);
      }
    }
    return trees;
  }
  return trees;
}
