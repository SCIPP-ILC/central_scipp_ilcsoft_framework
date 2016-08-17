#include <string>
#include <EVENT/MCParticle.h>
#include <EVENT/LCCollection.h>
#include "lcio.h"

namespace scipp_ilc {
    bool get_detectable_signal_event(lcio::LCEvent* signal_event, lcio::MCParticle*& electron);

    void transform_to_cm(double pX, double E, double& pX_new, double& E_new);

    void transform_to_lab(double pX, double E, double& pX_new, double& E_new);

    int get_hitStatus(double x, double y);
}
