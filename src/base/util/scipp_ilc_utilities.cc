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


    /*double transform_to_cm(double pX, double pY, double pZ, double E){
	    double theta = inc_ang;
	    double in_E = E_e;
	    double beta = Math.sin(theta);
	    double gamma = Math.pow((1-Math.pow(beta, 2)), -0.5);
	    double* cm_Vec;

	     *
	     *    |gamma         -gamma*beta| |p|                       |p'|
	     *    |-gamma*beta         gamma|*|E| = TRANSOFORMED4vector |E'|
	     *        
	     *                                     *                                     *      

	    cm_Vec[0] = pX*gamma - gamma*beta*E;
	    cm_Vec[1] = pY;
	    cm_Vec[2] = pZ;
	    cm_Vec[3] = E*gamma - gamma*beta*pX;
	    return cm_Vec;
    }    
    double transform_to_lab(double pX, double pY, double pZ, double E){
	    double theta = inc_ang;
	    double in_E = E_e;
	    double beta = Math.sin(theta);
	    double gamma = Math.pow((1-Math.pow(beta, 2)), -0.5);
	    double* cm_Vec;

	     *
	     *    |gamma         -gamma*beta| |p|                       |p'|
	     *    |-gamma*beta         gamma|*|E| = TRANSOFORMED4vector |E'|
	     *    *  
	     *                                     *                                     *      

	    cm_Vec[0] = pX*(-gamma) + gamma*beta*E;
	    cm_Vec[1] = pY;
	    cm_Vec[2] = pZ;
	    cm_Vec[3] = E*(-gamma) + gamma*beta*pX;
	    return cm_Vec;
    }*/

}
