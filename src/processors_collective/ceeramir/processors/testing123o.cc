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

//#include <ctime>    //************************************************************
#include <iostream>
#include <chrono>
#include <string>
#include <sstream>
#include <unordered_map>
#include <cmath>

#include "testing123o.h"
#include "scipp_ilc_utilities.h"
#include "polar_coords.h"
//#include "beamcal_reconstructor.h"
#include "include/beamcal_reconstructor_xy.h"
//#include "beamcal_scanner.h"
#include "include/beamcal_scanner_xy.h"

#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/MCParticle.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

// ----- all for ploting -----
#include <TFile.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPaveStats.h>

#include <TStyle.h>
#include <TColor.h>
#include <TLegend.h>
#include <TH1I.h>

typedef std::chrono::high_resolution_clock Clock;

using namespace lcio;
using namespace marlin;
using namespace std;


testing123o testing123o;

static TFile* _rootfile;
static TProfile* _radeff;
static int _detected_num = 0;
static int _test_num = 0;

static TProfile2D* _hitmap_bgd;
static TProfile2D* _hitmap_zeros;
static TProfile2D* _test_slice;

static TFile* _f1;
static TFile* _f2;

static TH2F* _h1;
static TH2F* _h2;
static TH2F* _h1_test;
static TH2F* _h2_test;

static TCanvas* _c1;
static TCanvas* _c2;




testing123o::testing123o() : Processor("testing123o") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
    registerProcessorParameter( "BeamcalGeometryFile" , "input file"  , _beamcal_geometry_file_name , std::string("input.xml") ) ;
    registerProcessorParameter( "BackgroundEventList" , "input file"  , _background_event_list , std::string("input.xml") ) ;
    registerProcessorParameter( "BackgroundEventsToRead" , "number"  , _num_bgd_events_to_read , 10 ) ;
    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}


/*
void testing123o::RootPlotTitle(stringstream s1, string ){
    s1.str("");
}
 */


void testing123o::RootPlot(TH2F* graph){                              // This function edits a root plot passed from init

  graph->GetXaxis()->SetTitle("X axis (mm)");
  graph->GetYaxis()->SetTitle("Y axis (mm)");
  graph->GetZaxis()->SetTitle("Efficiency");

  graph->GetXaxis()->CenterTitle();
  graph->GetYaxis()->CenterTitle();
  graph->GetZaxis()->CenterTitle();

  graph->GetXaxis()->SetTitleOffset(1.4);
  graph->GetYaxis()->SetTitleOffset(1.6);

}



void testing123o::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile(_root_file_name.c_str(),"RECREATE");
    _radeff = new TProfile("radeff","Radial Efficiency",14*2,0.0,140.0,0.0,1.0);

    _hitmap_bgd = new TProfile2D("hitmap_bgd","Hit Distribution",300.0,-150.0,150.0,300.0,-150.0,150.0);
    _hitmap_zeros = new TProfile2D("hitmap_zeros","Hit Distribution",300.0,-150.0,150.0,300.0,-150.0,150.0);
    _test_slice = new TProfile2D("hitmap_slice","Hit Distribution",300.0,-150.0,150.0,300.0,-150.0,150.0);

    _c1 = new TCanvas("c1","c1",600,400);


    int LEGObins = 60;
    string bgd_events = "bgd,";

    std::stringstream s1;

    // ------ LEGO GRAPHS 1-2 ------
    s1 << "LEGO 1s,"<< _num_bgd_events_to_read << "events," << LEGObins << "bin";
    const char* LEGOtitle = s1.str().c_str();
    _h1 = new TH2F("h1", LEGOtitle ,LEGObins ,-150,150,LEGObins,-150,150);
    RootPlot(_h1);

    s1.str("");
    s1 << "LEGO 0s,"<< _num_bgd_events_to_read << "events," << LEGObins << "bin";
    const char* LEGOtitle2 = s1.str().c_str();
    _h2 = new TH2F("h2", LEGOtitle2, LEGObins,-150,150,LEGObins,-150,150);
    RootPlot(_h2);
    // ------ end LEGO GRAPH 1-2 ------



    // ------ LEGO test GRAPHS 1 & 2 ------
    s1.str("");
    s1 << "LEGO 0s,"<< _num_bgd_events_to_read << "events," << LEGObins << "bin";
    const char* LEGOTitleTest1 = s1.str().c_str();
    _h1_test = new TH2F("h1_test", LEGOTitleTest1, LEGObins,-150,150,LEGObins,-150,150);
    RootPlot(_h1_test);

    s1.str("");
    s1 << "LEGO 0s,"<< _num_bgd_events_to_read << "events," << LEGObins << "bin";
    const char* LEGOTitleTest2 = s1.str().c_str();
    _h2_test = new TH2F("h2_test", LEGOTitleTest2, LEGObins,-150,150,LEGObins,-150,150);
    RootPlot(_h2_test);
    // ------ end test GRAPHS 1 & 2 ------  


    _f1 = TFile::Open("~/ilc_main/output/BeamCalRecon_xy_50_events.root");
    _c1 = (TCanvas*)_f1->Get("c1");
    _h1 = (TH2F*)_c1->FindObject("hlego");

    _f2 = TFile::Open("~/ilc_main/output/BeamCalRecon_xy_50bgd_varpolar.root");
    _c2 = (TCanvas*)_f2->Get("c1");
    _h2 = (TH2F*)_c2->FindObject("hlego");


    _nRun = 0 ;
    _nEvt = 0 ;
}



void testing123o::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void testing123o::processEvent( LCEvent* signal_event ) {
    //Make sure we are using an electron that actually hits the Positive BeamCal
  MCParticle* electron = NULL;
  bool detectable_electron = scipp_ilc::get_detectable_signal_event(signal_event,electron);
  if ( not detectable_electron ) return;


  scipp_ilc::testing123o::beamcal_cluster* signal_cluster;
  signal_cluster = scipp_ilc::testing123o::reconstruct_beamcal_event(signal_event);
  bool detected = signal_cluster->exceeds_sigma_cut;

  _detected_num += detected;
  if(_detected_num%100==0){
    cout << _detected_num << endl;
  }
}



void testing123o::check( LCEvent * evt ){
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void testing123o::end(){ 
    _h1_test->Add(_h1);
    _h1_test->Add(_h2,-1);

    _h2_test->Add(_h2);
    _h2_test->Add(_h1,-1);

    /*  _hlego_test->Add(_hlego_zeros);
  _hlego_test->Add(_hlego);
  _hlego_inefficiency->Divide(_hlego_test);
    */
  cout << "\n ********* end *********" << endl;    
    _rootfile->Write();
}
