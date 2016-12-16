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

#include "Thrust.h"
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

// #include <CLHEP/Vector/ThreeVector.h>
// #include <CLHEP/Random/RanluxEngine.h>

#include <EVENT/LCCollection.h>
#include <EVENT/LCIO.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/LCTOOLS.h>


using namespace lcio;
using namespace marlin;
using namespace std;


static TFile* _rootfile;

Thrust Thrust;

Thrust::Thrust() : Processor("Thrust") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
    registerProcessorParameter( "typeOfThrustFinder" ,
      "Type of thrust reconstruction algorithm to be used:\n#\t1 : Tasso algorithm\n#\t2 : JetSet algorithm"  ,
      _typeOfThrustFinder , 2 ) ;
}



void Thrust::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("thrust_test.root","RECREATE");
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



void Thrust::processRunHeader( LCRunHeader* run) { 
    run->parameters().setValue("thrust",12300321);
    //    _nRun++ ;

} 



void Thrust::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    _inParVec = evt->getCollection( _colName) ;

    if (!_partMom.empty()) _partMom.clear();

    for (int n=0;n<_inParVec->getNumberOfElements() ;n++)
    {
      MCParticle* aPart = dynamic_cast<MCParticle*>( _inParVec->getElementAt(n) );

      const double* partMom = aPart->getMomentum();
      _partMom.push_back( Hep3Vector(partMom[0], partMom[1], partMom[2]) ); 

    _nEvt ++ ;
    }

  //reset variables for output   
  _principleThrustValue = -1;
  _majorThrustValue     = -1;
  _minorThrustValue     = -1;
  _principleThrustAxis.set(0,0,0);
  _majorThrustAxis.set(0,0,0);
  _minorThrustAxis.set(0,0,0);

  // Switch to the desired type of thrust finder
  if (_typeOfThrustFinder == 1)
    { 
      TassoThrust();
    }
  else if (_partMom.size()<=1)
    { 
      TassoThrust();
    }
  else if (_typeOfThrustFinder == 2)
    {
      JetsetThrust();
    }
  // ###write
  //    evt->parameters().setValue("thrust",_principleThrustValue);

  FloatVec thrax;
  thrax.clear();
  thrax.push_back(_principleThrustAxis.x());
  thrax.push_back(_principleThrustAxis.y());
  thrax.push_back(_principleThrustAxis.z());

  _inParVec->parameters().setValue("principleThrustValue",_principleThrustValue);
  _inParVec->parameters().setValues("principleThrustAxis",thrax);

  if (_typeOfThrustFinder == 2)
    {
      thrax.clear();
      thrax.push_back(_majorThrustAxis.x());
      thrax.push_back(_majorThrustAxis.y());
      thrax.push_back(_majorThrustAxis.z());

      _inParVec->parameters().setValue("majorThrustValue",_majorThrustValue);
      _inParVec->parameters().setValues("majorThrustAxis",thrax);

      thrax.clear();
      thrax.push_back(_minorThrustAxis.x());
      thrax.push_back(_minorThrustAxis.y());
      thrax.push_back(_minorThrustAxis.z());

      _inParVec->parameters().setValue("minorThrustValue",_minorThrustValue);
      _inParVec->parameters().setValues("minorThrustAxis",thrax);

      float Oblateness;
      Oblateness = _majorThrustValue - _minorThrustValue;
      _inParVec->parameters().setValue("Oblateness",Oblateness);
      if ( (_majorThrustValue < 0) || (_minorThrustValue < 0) )
    {
      _inParVec->parameters().setValue("Oblateness",-1);
    }
    }

    streamlog_out( DEBUG4 ) << " thrust: " << _principleThrustValue << " TV: " << _principleThrustAxis << endl;
    streamlog_out( DEBUG4 ) << "  major: " << _majorThrustValue << " TV: " << _majorThrustAxis << endl;
    streamlog_out( DEBUG4 ) << "  minor: " << _minorThrustValue << " TV: " << _minorThrustAxis << endl;

  if (_principleThrustValue >= _max) _max = _principleThrustValue;
  if (_principleThrustValue <= _min) _min = _principleThrustValue;
}




void Thrust::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void Thrust::end(){ 
    _rootfile->Write();
}
