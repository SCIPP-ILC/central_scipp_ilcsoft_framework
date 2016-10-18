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

#include "BeamCalReconstruction.h"
#include "scipp_ilc_utilities.h"
#include "polar_coords.h"
#include "beamcal_reconstructor.h"
#include "beamcal_scanner.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/MCParticle.h>


// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#include <TFile.h>
#include <TProfile.h>



using namespace lcio;
using namespace marlin;
using namespace std;


BeamCalReconstruction BeamCalReconstruction;

static TFile* _rootfile;
static TProfile* _radeff;
static int _detected_num = 0;

BeamCalReconstruction::BeamCalReconstruction() : Processor("BeamCalReconstruction") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
    registerProcessorParameter( "BeamcalGeometryFile" , "input file"  , _beamcal_geometry_file_name , std::string("input.xml") ) ;
    registerProcessorParameter( "BackgroundEventList" , "input file"  , _background_event_list , std::string("input.xml") ) ;
    registerProcessorParameter( "BackgroundEventsToRead" , "number"  , _num_bgd_events_to_read , 10 ) ;
    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}



void BeamCalReconstruction::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile(_root_file_name.c_str(),"RECREATE");
    _radeff = new TProfile("radeff","Radial Efficiency",14,0.0,140.0,0.0,1.0);

    //Load up all the bgd events, and initialize the reconstruction algorithm.
    scipp_ilc::beamcal_recon::initialize_beamcal_reconstructor(_beamcal_geometry_file_name, _background_event_list, _num_bgd_events_to_read);

    _nRun = 0 ;
    _nEvt = 0 ;

}



void BeamCalReconstruction::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void BeamCalReconstruction::processEvent( LCEvent* signal_event ) { 
    //Make sure we are using an electron that actually hits the Positive BeamCal
    MCParticle* electron = NULL;
    bool detectable_electron = scipp_ilc::get_detectable_signal_event(signal_event,electron);
    if ( not detectable_electron ) return;

    //Get the radius at which the signal electron hit
    const double* endpoint = electron->getEndpoint();
    double end_x = (endpoint[0] - 0.007*endpoint[2]);
    double end_y = endpoint[1];
    double radius,phi;
    scipp_ilc::cartesian_to_polar(end_x,end_y,radius,phi);


    //Perform the reconstrunction algorithm, determine if the algorithm
    //detected the electron.
    scipp_ilc::beamcal_recon::beamcal_cluster* signal_cluster;
    signal_cluster = scipp_ilc::beamcal_recon::reconstruct_beamcal_event(signal_event);
    bool detected = signal_cluster->exceeds_sigma_cut;


    //Plot our results with respect to the radius of the signal electron.
    _radeff->Fill(radius,detected); //bools and ints are basically interchangeable...
    _detected_num += detected;
    

    cout << _nEvt++ << endl;;
}



void BeamCalReconstruction::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void BeamCalReconstruction::end(){ 
    cout << "\ndetected: " << _detected_num << endl;
    _rootfile->Write();
}
