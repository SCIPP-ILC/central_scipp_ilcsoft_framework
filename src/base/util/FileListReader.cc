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

#include "FileListReader.h"

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



    FileListReader aFileListReader ;


    FileListReader::FileListReader() : DataSourceProcessor("FileListReader") {

        _description = "Reads slcio filelists and extracts LCIO events with MCParticle collections."
            " Make sure to not specify any LCIOInputFiles in the steering in order to read StdHep files." ;

        registerProcessorParameter( "FileListName" , "input file"  , _fileName , std::string("input.list") ) ;

    }

    FileListReader*  FileListReader::newProcessor() { 
        return new FileListReader ;
    }

    void FileListReader::init() {    
        printParameters() ;    
    }


    void FileListReader::readDataSource( int numEvents ) {
        int numEventsRead = 0;
        try { 
            ifstream filelist (_fileName, ifstream::in);
            LCReader* lcReader = LCFactory::getInstance()->createLCReader() ;
            string slcioFile;
            LCEvent* event = NULL;

            while ( filelist >> slcioFile ) {
                LCRunHeaderImpl* runHeader = new LCRunHeaderImpl ;
                runHeader->setDescription( " Events read from input file list: " + _fileName ) ; 
                runHeader->setRunNumber( numEventsRead ) ;
                ProcessorMgr::instance()->processRunHeader( runHeader ) ;

                lcReader->open(slcioFile);
                while( (event=lcReader->readNextEvent()) ) {
                    ProcessorMgr::instance()->processEvent( event ) ;

                    numEventsRead++;
                    if ( numEventsRead >= numEvents ) {
                        delete event;
                        break;
                    }
                }
                lcReader->close();
                if ( numEventsRead >= numEvents ) break;
            }
            filelist.close();
        } catch(IOException& e) {
            cout << " Unable to read and analyze the LCIO file - " << e.what() << endl ;
        }
    }



    void FileListReader::end() {

    }


}
