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
typedef std::chrono::high_resolution_clock Clock;

#include "BeamCalRecon_edit.h"
#include "scipp_ilc_utilities.h"
#include "polar_coords.h"
//#include "beamcal_reconstructor.h"
#include "include/beamcal_reconstructor_edit.h"

//#include "beamcal_scanner.h"
#include "include/beamcal_scanner_edit.h"

#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/MCParticle.h>


// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#include <TFile.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TH2D.h>


using namespace lcio;
using namespace marlin;
using namespace std;


BeamCalRecon_edit BeamCalRecon_edit;

static TFile* _rootfile;
static TProfile* _radeff;
static int _detected_num = 0;
static int _test_num = 0;
//static TH2F* _hitmap_bgd;
//static TH2D* _hitmap_bgd;
static TProfile2D* _hitmap_bgd;
static TProfile2D* _hitmap_zeros;
static TProfile2D* _test_slice;

vector<pair<float,float>> bgd_plot;


//auto _t1 = Clock::now();
//auto _t2 = Clock::now();


BeamCalRecon_edit::BeamCalRecon_edit() : Processor("BeamCalRecon_edit") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
    registerProcessorParameter( "BeamcalGeometryFile" , "input file"  , _beamcal_geometry_file_name , std::string("input.xml") ) ;
    registerProcessorParameter( "BackgroundEventList" , "input file"  , _background_event_list , std::string("input.xml") ) ;
    registerProcessorParameter( "BackgroundEventsToRead" , "number"  , _num_bgd_events_to_read , 10 ) ;
    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}



void BeamCalRecon_edit::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile(_root_file_name.c_str(),"RECREATE");
    _radeff = new TProfile("radeff","Radial Efficiency",14*2,0.0,140.0,0.0,1.0);
    _hitmap_bgd = new TProfile2D("hitmap_bgd","Hit Distribution",300.0,-150.0,150.0,300.0,-150.0,150.0);
    _hitmap_zeros = new TProfile2D("hitmap_zeros","Hit Distribution",300.0,-150.0,150.0,300.0,-150.0,150.0);
    _test_slice = new TProfile2D("hitmap_slice","Hit Distribution",300.0,-150.0,150.0,300.0,-150.0,150.0);

    // TH2F
    // TProfile2D
    //    _t1 = Clock::now();

    //Load up all the bgd events, and initialize the reconstruction algorithm.
    scipp_ilc::beamcal_recon_C::initialize_beamcal_reconstructor(_beamcal_geometry_file_name, _background_event_list, _num_bgd_events_to_read);

    _nRun = 0 ;
    _nEvt = 0 ;

}



void BeamCalRecon_edit::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void BeamCalRecon_edit::processEvent( LCEvent* signal_event ) {
    //Make sure we are using an electron that actually hits the Positive BeamCal

  //  _hitmap_bgd->Fill

    MCParticle* electron = NULL;
    bool detectable_electron = scipp_ilc::get_detectable_signal_event(signal_event,electron);
    if ( not detectable_electron ) return;

    //Get the radius at which the signal electron hit
    const double* endpoint = electron->getEndpoint();
    double end_x = (endpoint[0] - 0.007*endpoint[2]);
    double end_y = endpoint[1];
    double endx = end_x;
    double endy = end_y;

    double radius,phi;
    scipp_ilc::cartesian_to_polar(end_x,end_y,radius,phi);


    //Perform the reconstrunction algorithm, determine if the algorithm
    //detected the electron.
    scipp_ilc::beamcal_recon_C::beamcal_cluster* signal_cluster;
    signal_cluster = scipp_ilc::beamcal_recon_C::reconstruct_beamcal_event(signal_event);
    bool detected = signal_cluster->exceeds_sigma_cut;

    if(!detected){
      cout << "detected:  " << detected<< "        x-y: "<< endx << "\t" <<endy << endl;
    }
    //Plot our results with respect to the radius of the signal electron.
    _radeff->Fill(radius,detected); //bools and ints are basically interchangeable...
    _detected_num += detected;

    if(detected && endx > 0 && endy < 0){
      _test_num += detected;
      _test_slice->Fill(endx,endy,detected);
    }

    if(detected){
      _hitmap_bgd->Fill(endx,endy,detected);
      //      _hitmap_bgd->Fill(endx,endy);
    }else{
      _hitmap_zeros->Fill(endx,endy,true);
    }
    cout << _nEvt++ << endl;
}

void BeamCalRecon_edit::check( LCEvent * evt ){
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void BeamCalRecon_edit::end(){ 
    cout << "\ndetected: " << _detected_num << endl;
    cout << "\n in \'slice\' of beamcal: " << _test_num << endl;

    //    _t2 = Clock::now();
    //    cout << "*******************this is the end***********************" << endl;
    //    cout << "******************* time elapsed: " << (_end - _begin) << " ***********************" << endl;
    //    cout << "******************* time elapsed: " << std::chrono::duration_cast<std::chrono::nanoseconds>(_t2 - _t1).count() << " ***********************" << endl;
    _rootfile->Write();
}
