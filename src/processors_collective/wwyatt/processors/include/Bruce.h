#ifndef overlay
#define overlay

#include <Overlay/Overlay.h>
#include "lcio.h"
#include <string>
#include <EVENT/MCParticle.h>



using namespace lcio;
using namespace marlin;

class Bruce : public Overlay{
 public:
  virtual marlin::Processor*  newProcessor() { return new Overlay ; }
  Overlay() ;
  Overlay( Overlay const&) = delete;
  Overlay& operator=(Overlay const&) = delete;
  virtual const std::string & name() const { return Processor::name() ; }
  virtual void modifyEvent( EVENT::LCEvent * evt ) ; 
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void check( EVENT::LCEvent * evt ) ; 
  virtual void end() ;
 private:
  /** Helper method.
   */
  LCEvent* readNextEvent( ) ;

  /** Input file names.
   */
  StringVec _fileNames{};
  int       _numOverlay = 0;
  double    _expBG = 1;
  bool      _runOverlay = false;
  int       _nSkipEventsRandom = 0 ;

  StringVec _colVec{};
  std::map<std::string, std::string> _colMap{};

  IO::LCReader* _lcReader = NULL;
  EVENT::LCEvent* _overlayEvent = NULL;

  int _activeRunNumber = 0;
  int _nRun = 0;
  int _nEvt = 0;
  int _nOverlayEvt = 0;
  IntVec _events{};
} ;

#endif
