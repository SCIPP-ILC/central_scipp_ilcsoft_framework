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

#include "StdhepFileListReader.h"

#include "marlin/ProcessorMgr.h"

#include "IMPL/LCEventImpl.h"
#include "EVENT/LCEvent.h"
#include "IMPL/LCRunHeaderImpl.h"

#include "UTIL/LCStdHepRdr.h"
#include "UTIL/LCTOOLS.h"

#include <iostream>
#include <fstream>
using namespace std;


namespace marlin{



    StdhepFileListReader aStdhepFileListReader ;


    StdhepFileListReader::StdhepFileListReader() : DataSourceProcessor("StdhepFileListReader") {

        _description = "Reads stdhep filelists and extracts LCIO events with MCParticle collections."
            " Make sure to not specify any LCIOInputFiles in the steering in order to read StdHep files." ;

        registerProcessorParameter( "FileListName" , "input file"  , _fileName , std::string("input.list") ) ;

    }

    StdhepFileListReader*  StdhepFileListReader::newProcessor() { 
        return new StdhepFileListReader ;
    }

    void StdhepFileListReader::init() {    
        printParameters() ;    
    }


    void StdhepFileListReader::readDataSource( int numEvents ) {
        try { 
            ifstream filelist (_fileName, ifstream::in);
            string stdhepFile;
            
            cout << "file: " << _fileName << "    numEvents: " << numEvents << endl;
            while ( filelist >> stdhepFile ) {
                LCStdHepRdr* rdr = new  LCStdHepRdr( stdhepFile.c_str() ) ;

                LCCollection* col ;
                LCEventImpl* evt ;

                int evtNum = 0 ;
                int runNum = 0 ;

                while( ( col = rdr->readEvent() ) != 0 ) {

                    if ( numEvents > 0 && evtNum+1 > numEvents )
                    {
                        delete col;
                        break;
                    }

                    if ( isFirstEvent() ) {   // create run header

                        LCRunHeaderImpl* rHdr = new LCRunHeaderImpl ;

                        rHdr->setDescription( " Events read from stdhep input file: " + _fileName ) ; 
                        rHdr->setRunNumber( runNum ) ;

                        ProcessorMgr::instance()->processRunHeader( rHdr ) ;
                        _isFirstEvent = false ;	
                    }

                    evt = new LCEventImpl ;
                    evt->setRunNumber( runNum ) ;
                    evt->setEventNumber( evtNum++ ) ;


                    evt->addCollection(  col, "MCParticle"  ) ;

                    ProcessorMgr::instance()->processEvent( evt ) ;

                    delete evt ;
                }

                delete rdr;
            }
            filelist.close();
        } catch(IOException& e) {
            cout << " Unable to read and analyze the STDHEP file - " << e.what() << endl ;
        }
    }



    void StdhepFileListReader::end() {

    }


}
