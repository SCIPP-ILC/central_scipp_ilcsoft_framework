#ifndef JetRazor_h
#define JetRazor_h 1
#include <vector>
#include "marlin/Processor.h"
#include "lcio.h"
#include <iostream>
#include <string>
#include <IMPL/ReconstructedParticleImpl.h>

#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Random/RanluxEngine.h>

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TObjArray.h"
#include "TArrayI.h"
#include "TClass.h"


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
 * @version $Id: JetRazor.h,v 1.4 2005-10-11 12:57:39 gaede Exp $ 
 */

class JetRazor : public Processor {

    public:

        virtual Processor*  newProcessor() { return new JetRazor ; }


        JetRazor() ;

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
        JetFinder(Double_t ycut = 0.);
        virtual ~JetFinder();

        void setPartList(TObjArray* e);    // Input the particle 4(3)-vector list
        // e: 4-vector  TLorentzVector ..(px,py,pz,E) or
        //    3-vector  TVector3       ..(px,py,pz) 
        // If you input TVector3, the energy of particle
        // will be E = sqrt(px**2 + py**2 + pz**2) 
        void setYCut(Double_t ycut);       // Set the YCut value

        void doFindJets();                 // Clustering the particles into Jets  

        Int_t njets(){ return m_njets; };     // The number of jets found.  
        TLorentzVector* jet4vec(Int_t index); // Return the 4 vector of a jet (particle sum).
        // index: The index of the jet of interest  
        Int_t nParticlesPerJet(Int_t index);  // Number of particles in a particular jet
        // index: The index of the jet of interest  
        TArrayI getPartIndex(){ return m_ipart_jet_assoc; };    // Return the particle index array.
        // m_ipart_jet_assoc[i] = j means 
        // particle i is placed in jet j.  
        Int_t fewestParts(){ return m_ifewest_parts; }; // minimum number of particles to make a jet.  
        Double_t getYCut() { return m_dycut; }; // Obtain the current ycut

        void setDURHAM(); // Select DURHAM algorithm
        void setJADE();   //        JADE   algorithm
        void setJADEE();  //        JADE E algorithm

        Double_t calcinvmass(const TLorentzVector &jet1,
        const TLorentzVector &jet2);

    private:
        Int_t m_njets;     // Number of jets found  
        TObjArray m_jet;   // m_jet[i] is the 4 vector sum of all the particles in jet i.  
        TArrayI m_ipart_jet_assoc; // m_ipart_jet_assoc[i] = j means particle i was placed in jet j.  
        TArrayI m_inparts_per_jet; // m_inparts_per_jet[i] = j means jet i has j particles in it.  
        Int_t m_ifewest_parts; // m_ifewest_parts is number of particles in the jet 
        // with the fewest particles.  
        Double_t m_evis;
        Double_t m_dycut;
        Int_t m_algorithm; // Algorithm used in Jet clustering

        const static Int_t UNASSOC;
        const static Int_t DURHAM;
        const static Int_t JADE;
        const static Int_t JADEE;


    protected:

        /** Input collection name.
        */
        std::string _colName ;
        std::string _root_file_name;
        int TassoThrustRazor();
        int JetsetJetRazor();
        double sign(double a,double b);
        double min(double a,double b);

        /** Input collection name.
        */
        std::string partMomCheck;
        std::string betaCheck;

        bool partMom1;
        bool partMom0; 

        float _min,_max;
        LCCollection* _inParVec;
        std::vector<Hep3Vector> _partMom;
        std::string filename;
        RanluxEngine myrnd;

        int _nRun ;
        int _nEvt ;

        TObjArray* m_4vec;

        ClassDef(JetFinder,1) // Jetfinder base class   
};

#endif



