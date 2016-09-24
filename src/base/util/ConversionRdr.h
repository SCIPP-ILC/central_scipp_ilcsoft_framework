#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0

#ifndef UTIL_LCStdHepReader_H
#define UTIL_LCStdHepReader_H 1

class MCParticleCollection ;

#include <iostream>

namespace plcio{

  class lStdHep ;

  
  /**Basic utility for reading a binary stdhep file and filling
   * a LCCollectionVec with MCParticles containing the stdhep
   * file information.
   * 
   * @author cassell, F Gaede
   * @version $Id: $
   */
  class LCStdHepReader{
    
  public:

	/** Open the stdhep input file in the constructer
	 */
    LCStdHepReader(const char* evfile) ;

	/** noop
	 */
	~LCStdHepReader() ;

    /** Get number of events in the stdhep file.
     *  This number is read from the file header (no guarantee that it is correct)
     */
    long getNumberOfEvents() const ;

    /** Get total number of expected events in the whole set of stdhep files
     *  from which this stdhep file belongs to.
     *  This number is read from the file header (no guarantee that it is correct)
     */
    long getNumberOfTotalEventsExpected() const ;

    /** Reads the next stdhep event and fills the MCParticleCollection
     *  Returns 0 if EOF.
     */
    int readEvent( MCParticleCollection* ) ;


    /** Print the file header to the given ostream.
     */
    void printHeader(std::ostream& os = std::cout ) ; 


    /** Return the charge of the particle times 3  - code copied from HepPDT package.
     */
    int threeCharge( int pdgID ) const ;

  private:
    
	lStdHep* _reader;
    

  }; // class

} // namespace UTIL

#endif /* ifndef UTIL_LCStdHepReader_H */
