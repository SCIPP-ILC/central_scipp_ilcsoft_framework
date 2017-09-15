#ifndef BeamCalRecon_xy_h
#define BeamCalRecon_xy_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <TH2.h>
#include <vector>

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
 * @version $Id: BeamCalRecon_xy.h,v 1.4 2005-10-11 12:57:39 gaede Exp $ 
 */

class BeamCalRecon_xy : public Processor {

    public:

        virtual Processor*  newProcessor() { return new BeamCalRecon_xy ; }

        BeamCalRecon_xy() ;

	// Used to edit the LEGO plots
        virtual void PlotTH2F(TH2F* graph, std::stringstream& stream, std::string title, int bins);
	virtual void PlotTH1F(TH1F* graph, std::stringstream& stream, std::string title);

	//	virtual void FillRadiusThetaTable(vector<pair<float,float>>& thisTable, bool one);
	virtual void FillRadiusThetaTable(bool one, double radius, double phi_in_radian);

	//	virtual void PrintRadiusThetaTable();
	//	virtual void PrintRadiusThetaTable_two(int a[][13]);
	//	virtual void PrintRadiusThetaTable_two(int a[int][int]);
	virtual void PrintRadiusThetaTable(std::string key);    //this virtual works!!!!!1
	//	virtual void PrintRadiusThetaTable(std::string key, ofstream fout);

        /** Called at the begin of the job before anything is read.
         * Use to initialize the processor, e.g. book histograms.
         */
        virtual void init() ;

        /** Called for every run.        */
        virtual void processRunHeader( LCRunHeader* run ) ;

        /** Called for every event - the working horse.        */
        virtual void processEvent( LCEvent * evt ) ; 

        virtual void check( LCEvent * evt ) ; 

        /** Called after data processing for clean up.        */
        virtual void end() ;


    protected:

        /** Input collection name.        */
        std::string _colName ;
        std::string _beamcal_geometry_file_name;
        std::string _background_event_list;
        int _num_bgd_events_to_read;
        std::string _root_file_name;
	
	bool _test_bool;

        int _nRun ;
        int _nEvt ;
	double _max_radius;

} ;

#endif



