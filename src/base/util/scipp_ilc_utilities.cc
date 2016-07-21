#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0
#include <string>
#include <iostream>
#include <cmath>
#include "scipp_ilc_utilities.h"
#include "scipp_ilc_globals.h"
#include "polar_coords.h"

using namespace std;
using namespace lcio;

namespace scipp_ilc {
    static bool is_detectable(MCParticle* electron) {
        const double* endpoint = electron->getEndpoint();
        double end_z = endpoint[2];
        double end_y = endpoint[1];
        double end_x = endpoint[0] - abs(end_z)*_transform;

        double radius = hypot(end_x,end_y);


        double beamCalFront = 2000;
        double beamCalRear = 4000;
        return (beamCalFront < end_z and end_z < beamCalRear and radius < _radius_cut);


    }


    
    //Get electron/positron signal event and ensure it actually hits the detector
    bool get_detectable_signal_event(LCEvent* signal_event, MCParticle*& electron) {
        LCCollection* particles = signal_event->getCollection("MCParticle") ;
        if( particles == NULL ) return false;

        int num_particles = particles->getNumberOfElements();
        for( int particle_index = 0; particle_index < num_particles; particle_index++ ) {
            MCParticle* particle = dynamic_cast<MCParticle*>( particles->getElementAt( particle_index ) ) ;
            int pdgid = particle->getPDG();
            int status = particle->getGeneratorStatus(); //FINAL_STATE = 1
            if ( abs(pdgid) == 11 and status == 1 ) {
                electron = particle;
                return is_detectable(electron);
            }
        }

        return false;
    }


    void transform_to_cm(double pX, double E, double& pX_new, double& E_new){
	    double theta = 0.007;
	    double in_E = 250.0;
	    double beta = sin(theta);
	    double gamma = pow((1-pow(beta, 2)), -0.5);

	    // *
	    // *    |gamma         -gamma*beta| |p|                       |p'|
	    // *    |-gamma*beta         gamma|*|E| = TRANSOFORMED4vector |E'|
	    // *        
	    // *                                     *                                     *      

	    pX_new = pX*gamma - gamma*beta*E;
	    E_new = E*gamma - gamma*beta*pX;
    }    
    void transform_to_lab(double pX, double E, double& pX_new, double& E_new){
	    double theta = 0.007;
	    double in_E = 250.0;
	    double beta = sin(theta);
	    double gamma = pow((1-pow(beta, 2)), -0.5);

	     //*
	     //*    |gamma         gamma*beta| |p|                       |p'|
	     //*    |gamma*beta         gamma|*|E| = TRANSOFORMED4vector |E'|
	     //*        
	     //*                                     *                                     *      

	    pX_new = pX*gamma + gamma*beta*E;
	    E_new = E*gamma + gamma*beta*pX;
   }
}//END
