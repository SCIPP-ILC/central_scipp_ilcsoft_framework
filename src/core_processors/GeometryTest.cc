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

#include <iostream>

#include <vector>
#include <unordered_map>

#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/MCParticle.h>

#include "scipp_ilc_utilities.h"
#include "simple_list_geometry.h"

#include "GeometryTest.h"


// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"



using namespace lcio;
using namespace marlin;
using namespace std;

using namespace scipp_ilc::beamcal_recon;


GeometryTest GeometryTest;

static TFile* _rootfile;
static TH2F* _geomtest;

GeometryTest::GeometryTest() : Processor("GeometryTest") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
    registerProcessorParameter( "BeamcalGeometryFile" , "input file"  , _beamcal_geometry_file_name , std::string("input.xml") ) ;
    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}



void GeometryTest::init() { 
    initialize_geometry(_beamcal_geometry_file_name);

    float bound = 160;
    float res = 1000;
    float increment = 2*bound / res;

    _rootfile = new TFile(_root_file_name.c_str(),"RECREATE");
    _geomtest = new TH2F("geomtest","Geom Test",res,-bound,bound,res,-bound,bound);

    int prevID = 0;
    for (float x = -bound; x <= bound; x += increment) {
        for (float y = -bound; y <= bound; y += increment) {
            int ID = getID(x,y);
            
            if (ID != prevID) {
                prevID = ID;
                _geomtest->Fill(x,y);
            }
        }
    }
    cout << "finished x ...";

    prevID = 0;
    for (float y = -bound; y <= bound; y += increment) {
        for (float x = -bound; x <= bound; x += increment) {
            int ID = getID(x,y);
            
            if (ID != prevID) {
                prevID = ID;
                _geomtest->Fill(x,y);
            }
        }
    }
    cout << "finished y\n";
    _rootfile->Write();

}



void GeometryTest::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void GeometryTest::processEvent( LCEvent * evt ) { 
}



void GeometryTest::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void GeometryTest::end(){ 
}
