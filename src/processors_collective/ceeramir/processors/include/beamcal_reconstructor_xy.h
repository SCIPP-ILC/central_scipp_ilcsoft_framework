#ifndef BEAMCAL_RECONSTRUCTOR_H
#define BEAMCAL_RECONSTRUCTOR_H
#include <vector>
#include <unordered_map>
#include "lcio.h"

namespace scipp_ilc {
    namespace beamcal_recon_xy {
        typedef std::unordered_map<int,float> pixel_map;

        void initialize_beamcal_reconstructor(std::string geom_file_name, std::string bgd_list_file_name, int bgd_events_to_be_read);
        extern std::vector<pixel_map*>* _database;

        struct beamcal_cluster;
        beamcal_cluster*  reconstruct_beamcal_event(lcio::LCEvent* signal_event);
    }
}
#endif
