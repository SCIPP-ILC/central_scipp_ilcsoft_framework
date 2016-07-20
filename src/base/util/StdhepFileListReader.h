#ifndef StdhepFileListReader_h
#define StdhepFileListReader_h 1

#include "marlin/DataSourceProcessor.h"

using namespace lcio ;


namespace marlin{
  
  /** Reads binary StdHep files.
   *  Example processor for reading non-LCIO input files - creates events with
   *  MCParticle collections from binary StdHep files. Has to be the first active processor
   *  and requires that no LCIO input collection is used (parameter LCIOInputFiles).
   *
   *  <h4>Input - Prerequisites</h4>
   *  StdHep file.
   *
   *  <h4>Output</h4> 
   *  LCEvent with MCParticle collection.
   *
   * @param StdHepFileName   name of input file
   *
   * @author F. Gaede, DESY
   * @version $Id: StdhepFileListReader.h,v 1.3 2005-10-11 12:56:28 gaede Exp $ 
   */
  
  class StdhepFileListReader : public DataSourceProcessor {
    
  public:

    StdhepFileListReader() ;

    virtual StdhepFileListReader*  newProcessor() ;


    /** Creates events with MCParticle collections from the StdHep input file and
     *  calls all active processors' processEvent() and processRunHeader Method.
     *
     */
    virtual void readDataSource( int numEvents ) ;
    
    
    virtual void init() ;
    virtual void end() ;
    
  protected:
    
    std::string _fileName ;

  };
 
} // end namespace marlin 

#endif
