#ifndef Pred_Compare_h
#define Pred_Compare_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>


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
 * @version $Id: Pred_Compare.h,v 1.4 2005-10-11 12:57:39 gaede Exp $ 
 */

class Pred_Compare : public Processor {

    public:

        virtual Processor*  newProcessor() { return new Pred_Compare ; }


        Pred_Compare() ;

        /** Called at the begin of the job before anything is read.
         * Use to initialize the processor, e.g. book histograms.
         */
        virtual void init() ;

        /** Called for every run.
        */
        virtual void processRunHeader( LCRunHeader* run ) ;

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
        
        std::string interest;
        int total;
        int hh_hh, hh_hm, hh_mm, hm_hh, hm_hm, hm_mm, mm_hh, mm_hm, mm_mm;

        int _nRun ;
        int _nEvt ;
};

#endif



