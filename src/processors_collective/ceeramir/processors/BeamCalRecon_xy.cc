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

#include "BeamCalRecon_xy.h"
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


BeamCalRecon_xy BeamCalRecon_xy;

static TFile* _rootfile;
static TProfile* _radeff;
static int _detected_num = 0;
static int _test_num = 0;


//static TH2F* _hitmap_bgd;
//static TH2D* _hitmap_bgd;
static TProfile2D* _hitmap_bgd;
static TProfile2D* _hitmap_zeros;
static TProfile2D* _test_slice;
static TProfile2D* _hitmap_signal_electrons;

static TH1F* _1DRadHitsSigE_wCut;
static TH1F* _1DRadHitsSigE_wBGD;
static TH1F* _1DRadHitsSigE_wDiv;

static TH2F* _hlego;
static TH2F* _hlego_zeros;
static TH2F* _hlego_inefficiency;
static TH2F* _hlego_test;

//static TH2I* _h2;
static TH2F* _hlego_pol1;
static TH2F* _hlego_pol2;
static TH2F* _hlego_pol3;
static TH2F* _hlego_pol4;


static TH2D* _hcol1;

//static TGraphPolar* _hlegoo;

//static TH2F* _hlego_var;
//static TH2F* _hlego_zeros_var;
//static TH2F* _hlego_inefficiency_var;
//static TH2F* _hlego_test_var;

static TCanvas* _c2;
static TCanvas* _c1;

static unordered_map<pair<float,float>,double>* _all_map;
static unordered_map<pair<float,float>,double>* _zeros_map;
//vector<pair<float,float>> bgd_plot_xy;


//auto _t1 = Clock::now();
//auto _t2 = Clock::now();


//const bool _polar_coord_ID = true;

//const bool scipp_ilc::simple_list_geometry_xy::_polar_coord_ID = true;
//scipp_ilc::beamcal_recon_xy::

bool _test_bool;
bool _polar_coord_ID;

BeamCalRecon_xy::BeamCalRecon_xy() : Processor("BeamCalRecon_xy") {
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
void BeamCalRecon_xy::RootPlotTitle(stringstream s1, string ){
    s1.str("");
}
 */


void BeamCalRecon_xy::PlotTH2F(TH2F* graph){                              // This function edits a root plot passed from init

  graph->GetXaxis()->SetTitle("X axis (mm)");
  graph->GetYaxis()->SetTitle("Y axis (mm)");
  graph->GetZaxis()->SetTitle("Efficiency");

  graph->GetXaxis()->CenterTitle();
  graph->GetYaxis()->CenterTitle();
  graph->GetZaxis()->CenterTitle();

  graph->GetXaxis()->SetTitleOffset(1.4);
  graph->GetYaxis()->SetTitleOffset(1.6);

  //  _c1->SetStatX(0.15);

  //  graph->SetTheta(90);
  /* 
     gStyle->SetOptStat(1111111);            // Set stat options
     gStyle->SetStatY(0.9);                  // Set y-position (fraction of pad size)
     gStyle->SetStatX(0.9);                  // Set x-position (fraction of pad size)
     gStyle->SetStatW(0.4);                  // Set width of stat-box (fraction of pad size)
     gStyle->SetStatH(0.2);                  // Set height of stat-box (fraction of pad size)
  */

  //  graph->SetPhi(60);
  //  TPaveStats *st = (TPaveStats*)graph->FindObject("stats");
  //  st->SetX1NDC(0.0);
  //  st->SetX2NDC(0.0);
  //    gStyle->SetOptStat(0);
}

void BeamCalRecon_xy::PlotTH1F(TH1F* graph){                              // This function edits a root plot passed from init 
  graph->GetXaxis()->SetTitle("Radius (mm)");
  graph->GetYaxis()->SetTitle("e Count");
  graph->GetXaxis()->CenterTitle();
  graph->GetYaxis()->CenterTitle();
  graph->GetXaxis()->SetTitleOffset(1.0);
  graph->GetYaxis()->SetTitleOffset(1.0);

}

void BeamCalRecon_xy::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile(_root_file_name.c_str(),"RECREATE");
    _radeff = new TProfile("radeff","Radial Efficiency",14*2,0.0,140.0,0.0,1.0);
    _hitmap_bgd = new TProfile2D("hitmap_bgd","Hit Distribution",300.0,-150.0,150.0,300.0,-150.0,150.0);
    _hitmap_zeros = new TProfile2D("hitmap_zeros","Hit Distribution",300.0,-150.0,150.0,300.0,-150.0,150.0);
    _test_slice = new TProfile2D("hitmap_slice","Hit Distribution",300.0,-150.0,150.0,300.0,-150.0,150.0);
    _hitmap_signal_electrons = new TProfile2D("hitmap_es","Hit Distribution",300.0,-150.0,150.0,300.0,-150.0,150.0);

    //    _c2 = new TCanvas("c2","c2",300,300);
    _c1 = new TCanvas("c1","c1",600,400);

    _1DRadHitsSigE_wCut = new TH1F("radHitsWCut","Radial eHits",150,0,300);
    _1DRadHitsSigE_wBGD = new TH1F("radHitsWBgd","Radial eHits",150,0,300);
    _1DRadHitsSigE_wDiv = new TH1F("radHitsWDiv","Radial eHits",150,0,300);

    PlotTH1F(_1DRadHitsSigE_wCut);
    PlotTH1F(_1DRadHitsSigE_wBGD);
    PlotTH1F(_1DRadHitsSigE_wDiv);

    Double_t theta[8];
    Double_t radius[8];
    Double_t e_theta[8];
    Double_t e_radius[8];

    /*
    _hlegoo = new TgraphPolar(8, theta, radius, e_theta, e_radius);

    _hlegoo->SetMarkerStyle(20);
    _hlegoo->SetMarkerSize(20);
    _hlegoo->SetMarkerColor(20);
    _hlegoo->SetLineColor(20);
    _hlegoo->SetLineWidth(20);

    _h2 = new TH1I("h2","Gaus",100,-5,5);
    _h2->GetXaxis()->SetTitle("Standard deviation #sigma");
    _h2->GetYaxis()->SetTitle("dN/d#sigma");

    */

    int LEGObins = 60;
    int POLARbins = 40;
    const Int_t NBINS = 68;
    //    static float edges[] = { -199.50, -189.57, -181.02, -172.60,
    Double_t _var_edges[NBINS + 1] = { -199.50, -189.57, -181.02, -172.60,
			   -164.32,    -156.17,    -148.16,    -140.30,    -132.57,
			   -125.00,    -117.58,    -110.30,    -103.19,    -96.23,
			   -89.44,    -82.82,    -76.37,    -70.09,    -64.00,
			   -58.09,    -52.38,    -46.87,    -41.57,    -36.48,
			   -31.62,    -27.00,    -22.63,    -18.52,    -14.70,
			   -11.18,    -8.00,    -5.20,    -2.83,    -1.00,
			   0.00,    1.00,    2.83,    5.20,    8.00,
			   11.18,    14.70,    18.52,    22.63,    27.00,
			   31.62,    36.48,    41.57,    46.87,    52.38,
			   58.09,    64.00,    70.09,    76.37,    82.82,
			   89.44,    96.23,    103.19,    110.30,    117.58,
			   125.00,    132.57,    140.30,    148.16,    156.17,
			   164.32, 172.60, 181.02, 189.57, 199.50};

    std::stringstream s1;
    string bgd_events = "bgd,";    

    // ------ LEGO hits GRAPH ------
    s1 << "LEGO 1s,"<< _num_bgd_events_to_read << "events," << LEGObins << "bin";
    const char* LEGOtitle = s1.str().c_str();
    _hlego = new TH2F("hlego", LEGOtitle ,LEGObins ,-150,150,LEGObins,-150,150);
    PlotTH2F(_hlego);
    // ------ end LEGO GRAPH ------

    // ------ POLAR GRAPH ------
    s1.str("");
    s1 << "LEGO 1s,"<< _num_bgd_events_to_read << bgd_events << LEGObins << "bin";
    const char* LEGOtitle2 = s1.str().c_str();
    _hlego_pol1 = new TH2F("hlego_pol1", LEGOtitle2 ,POLARbins ,-150,150,POLARbins,-150,150);
    PlotTH2F(_hlego_pol1);
    // ------ end POLAR GRAPH ------

    // ------ LEGO non-hits GRAPH ------
    s1.str("");
    s1 << "LEGO 0s,"<< _num_bgd_events_to_read << "events," << LEGObins << "bin";
    const char* LEGOtitleZeros = s1.str().c_str();
    _hlego_zeros = new TH2F("hlego_0s", LEGOtitleZeros, LEGObins,-150,150,LEGObins,-150,150);
    PlotTH2F(_hlego_zeros);
    // ------ end LEGO non-hits GRAPH ------


    // ------ LEGO inefficiency GRAPH ------
    s1.str("");
    s1 << "LEGO,"<< _num_bgd_events_to_read << "events," << LEGObins << "bin";
    const char* LEGOtitleInefficiency = s1.str().c_str();
    _hlego_inefficiency = new TH2F("hlego_inefficiency", LEGOtitleInefficiency, LEGObins,-150,150,LEGObins,-150,150);
    //    _hlego_inefficiency_var = new TH2F("hlego_inefficiency_var", LEGOTitleInefficiency, 67, _var_edges, 67, _var_edges);
    PlotTH2F(_hlego_inefficiency);

    s1.str("");
    s1 << "LEGO , test"<< _num_bgd_events_to_read << "events," << LEGObins << "bin";
    const char* LEGOTestTitle = s1.str().c_str();
    _hlego_test = new TH2F("hlego_test", LEGOTestTitle ,LEGObins,-150,150,LEGObins,-150,150);
    //    _hlego_test_var = new TH2F("hlego_test_var", LEGOTestTitle, 67, _var_edges, 67, _var_edges);
    PlotTH2F(_hlego_test);
    // ------ end LEGO inefficiency GRAPH ------


    /*
      // ------ Variable binning ------
    s1.str("");
    s1 << "LEGO 1s, var"<< _num_bgd_events_to_read << "events," << LEGObins << "bin";
    const char* LEGOtitleVar = s1.str().c_str();
    _hlego_var = new TH2F("hlego_var", LEGOtitleVar, 67, _var_edges, 67, _var_edges);

    s1.str("");
    s1 << "LEGO 0s, var"<< _num_bgd_events_to_read << "events," << LEGObins << "bin";
    const char* LEGOtitlezVar = s1.str().c_str();
    _hlego_zeros_var = new TH2F("hlego_0s_var", LEGOtitlezVar, 67, _var_edges, 67, _var_edges);
     */


    //*********************************************************************************
    //    TCanvas *c1 = new TCanvas("c1","c1",600,400);
    //    TH2F *_hcol1 = new TH2F("hcol1","Option COLor combined with POL",40,-4,4,40,-4,4);
    _hcol1 = new TH2D("hcol1","Option COLor combined with POL",40,-150,150,40,-150,150);

    Float_t px, py;
    //    for (Int_t i = 0; i < 25000; i++) {
    //      gRandom->Rannor(px,py);

    //    }
    //    gStyle->SetPalette(kBird);

    //*********************************************************************************

    //    hcol1->Fill(px,py);
    //    hcol1->Draw("COLZPOL");


    TLegend *legend = new TLegend(0.05,0.05,0.06,0.1);


    // ------ clock start ------
    //    _t1 = Clock::now();

    //Load up all the bgd events, and initialize the reconstruction algorithm.

    scipp_ilc::beamcal_recon_xy::initialize_beamcal_reconstructor(_beamcal_geometry_file_name, _background_event_list, _num_bgd_events_to_read);

    _max_radius = 0.0;
    _nRun = 0 ;
    _nEvt = 0 ;
}



void BeamCalRecon_xy::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 


void BeamCalRecon_xy::processEvent( LCEvent* signal_event ) {
    //Make sure we are using an electron that actually hits the Positive BeamCal
    //  _hitmap_bgd->Fill
  _polar_coord_ID = true;
  _test_bool = true;
  //  cout << " BeamCal Recon test_bool: " << _test_bool << endl;

    MCParticle* electron = NULL;
    bool detectable_electron = scipp_ilc::get_detectable_signal_event(signal_event,electron);
    if ( not detectable_electron ) return;

    double electron_energy = electron->getEnergy();
    //Get the radius at which the signal electron hit
    const double* endpoint = electron->getEndpoint();
    double end_x = (endpoint[0] - 0.007*endpoint[2]);
    double end_y = endpoint[1];
    double endx = end_x;
    double endy = end_y;
    double radius,phi;
    scipp_ilc::cartesian_to_polar(end_x,end_y,radius,phi);

    //2D hitmap with ring at about 60mm
    _hitmap_signal_electrons->Fill(endx,endy,1);
    _1DRadHitsSigE_wCut->Fill(radius,detectable_electron);
    if(radius > _max_radius){
      _max_radius = radius;
    }

    //Perform the reconstrunction algorithm, determine if the algorithm
    //detected the electron.
    scipp_ilc::beamcal_recon_xy::beamcal_cluster* signal_cluster;
    signal_cluster = scipp_ilc::beamcal_recon_xy::reconstruct_beamcal_event(signal_event);
    bool detected = signal_cluster->exceeds_sigma_cut;
    _1DRadHitsSigE_wBGD->Fill(radius,detected);
    //      cout << "******************************************************************" << endl;
    //      cout << "detected:  " << detected<< "        x-y: "<< endx << "\t" <<endy << endl;
    //      cout << "******************************************************************" << endl;


    // ------ set up map ------
    pair<float,float> pos;
    pos.first = (float) endx;
    pos.second = (float) endy;
    
    string endx_s = std::to_string(pos.first);
    string endy_s = std::to_string(pos.second);
    //    cout << "endx string"<< endx_s << endl;

    string ID = endx_s + "," + endy_s;
    //    cout << "ID string"<< ID << endl;

    Float_t px,py;
    px = endx;
    py =endy;
    //	int ID = scipp_ilc::beamcal_recon_xy::getID(end_x,end_y);
    //	int ID = scipp_ilc::simple_list_geometry_xy::getID(end_x,end_y);
    // ------  map  ------


    //Plot our results with respect to the radius of the signal electron.
    _radeff->Fill(radius,detected);             //bools and ints are basically interchangeable...
    _detected_num += detected;

    if(detected && endx > 0 && endy < 0){       //Graph of slice of beamcal
      _test_num += detected;
      _test_slice->Fill(endx,endy,detected);
    }
    if(detected){                               //Graph of detected
      _hitmap_bgd->Fill(endx,endy,detected);    //      _hitmap_bgd->Fill(endx,endy);
      _hlego->Fill(endx,endy,detected);
      _hlego_pol1->Fill(px,py);
      _hcol1->Fill((float_t(px)),(float_t(py)));
      //      _h2->Fill(electron_energy);
      //      _hlego_var->Fill(endx,endy,detected);
    }else{                                      //Graph of not detected
      _hitmap_zeros->Fill(endx,endy,true);
      _hlego_zeros->Fill(endx,endy,true);
      //      _hlego_zeros_var->Fill(endx,endy,true);
      
      //      if((_zeros_map)[pos]>=1.0){
      //	(_zeros_map)[pos]+= 1.0;
      //      }else{
      //	(_zeros_map)[pos] = 1.0;
      //      }
    }
    //    if(!detected){                              //Print out xy-coords of not detected
    //      cout << "detected:  " << detected<< "        x-y: "<< endx << "\t" <<endy << endl;
    //    }
    //    (*_all_map)[ID]+= (detected || !detected);

    //    gStyle->SetPalette(kBird);
    _hlego->SetFillColor(kYellow);
    _nEvt++;
    if(_nEvt%100==0){
      cout << _nEvt << endl;
    }
}


void BeamCalRecon_xy::check( LCEvent * evt ){
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void BeamCalRecon_xy::end(){ 
  // ------ ------
  //  for(auto bit:*_zeros_map){
  //    pair<float,float> ID = bit.first;
  //    double bit_hit = bit.second;
  //    double all_hit = *_all_map[ID];
  //  _hlego_inefficiency->Fill(ID.first,ID.second,bit_hit/all_hit)
  //  }

  _hlego_inefficiency->Add(_hlego_zeros);
  /*  _hlego_test->Add(_hlego_zeros);
  _hlego_test->Add(_hlego);
  _hlego_inefficiency->Divide(_hlego_test);
  */
    cout << "\ndetected: " << _detected_num << endl;
    cout << "\n in \'slice\' of beamcal: " << _test_num << endl;
    cout << "max radius: " << _max_radius << endl;
    
    // ------ clock end ------
    //    _t2 = Clock::now();
    //    cout << "*******************this is the end***********************" << endl;
    //    cout << "******************* time elapsed: " << (_end - _begin) << " ***********************" << endl;
    //    cout << "******************* time elapsed: " << std::chrono::duration_cast<std::chrono::nanoseconds>(_t2 - _t1).count() << " ***********************" << endl;
    _hcol1->Draw("COLZPOL");   
    _rootfile->Write();
}
