#include <string>
#include <EVENT/MCParticle.h>
#include <EVENT/LCCollection.h>
#include "lcio.h"

namespace scipp_ilc {
    bool get_detectable_signal_event(lcio::LCEvent* signal_event, lcio::MCParticle*& electron);
}
