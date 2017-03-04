#ifndef ThrustSusy_h
#define ThrustSusy_h 1
#include <vector>
#include "marlin/Processor.h"
#include "lcio.h"
#include <iostream>
#include <string>
#include <IMPL/ReconstructedParticleImpl.h>

#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Random/RanluxEngine.h>

namespace CLHEP{}    // declare namespace CLHEP for backward compatibility
using namespace CLHEP ;
using namespace lcio ;
using namespace marlin ;


/**  Example processor for marlin.
 * 
 *  If compiled with MARLIN_USE_AIDA 
 *  it creates a histogram (cloud) of the MCParticle energies.
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Needs the collection of MCParticles.
 *
 *  <h4>Output</h4> 
 *  A histogram.
 * 
 * @param CollectionName Name of the MCParticle collection
 * 
 * @author F. Gaede, DESY
 * @version $Id: ThrustSusy.h,v 1.4 2005-10-11 12:57:39 gaede Exp $ 
 */

class ThrustSusy : public Processor {

    public:

        virtual Processor*  newProcessor() { return new ThrustSusy ; }


        ThrustSusy() ;

        /** Called at the begin of the job before anything is read.
         * Use to initialize the processor, e.g. book histograms.
         */
        virtual void init() ;

        /** Called for every run.
        */
        virtual void processRunHeader( LCRunHeader* run ) ;
        virtual void modifyRunHeader( LCRunHeader* run ) {}
        /** Called for every event - the working horse.
        */
        virtual void processEvent( LCEvent * evt ) ; 

        virtual void check( LCEvent * evt ) ; 


        /** Called after data processing for clean up.
        */
        virtual void end() ;


    protected:

        /** Input collection name.
        */
        std::string _colName ;
        std::string _root_file_name;
        int TassoThrustSusy();
        int JetsetThrustSusy();
        double sign(double a,double b);
        double min(double a,double b);

      /** Input collection name.
       */
        std::string partMomCheck;
        int _typeOfThrustSusyFinder;
        int _thrustDetectability;
        float _principleThrustSusyValue;
        float _majorThrustSusyValue;
        float _minorThrustSusyValue;
        Hep3Vector _principleThrustSusyAxis;
        Hep3Vector _majorThrustSusyAxis;
        Hep3Vector _minorThrustSusyAxis;
        float _min,_max;
        LCCollection* _inParVec;
        std::vector<Hep3Vector> _partMom;
        std::string filename;
        RanluxEngine myrnd;

        int _nRun ;
        int _nEvt ;
};

#endif



