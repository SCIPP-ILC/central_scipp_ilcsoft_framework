#ifndef BEAMCAL_SCANNER_H
#define BEAMCAL_SCANNER_H

#include <unordered_map>
#include <vector>

namespace scipp_ilc {
    namespace beamcal_recon_xy {

        struct beamcal_cluster {
            std::vector<int>* id_list;
            float significance;
            float energy;
            double background_average;

            //this value is left undefined until
            //the signal processing phase
            bool exceeds_sigma_cut;
        };



        beamcal_cluster* scan_beamcal(std::unordered_map<int,float>* pixels, std::unordered_map<int,double>* average_map, 
                                    std::unordered_map<int,double>* std_dev_map);
    }
}
#endif

