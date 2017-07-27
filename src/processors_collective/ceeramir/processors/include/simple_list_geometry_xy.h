#include <vector>
#include <unordered_map>

namespace scipp_ilc {
    namespace beamcal_recon_xy {
        struct surrounding_ids {
            std::vector<int>* list;

            //denotes the index of the element
            //half-way around the id
            int half_turn_index;
        };


        extern const int _IDlimit;
        extern int _LastRing;
        extern std::unordered_map<int,surrounding_ids*>* _pixel_graph;

	//	extern const bool _polar_coord_ID;
	extern bool _polar_coord_ID;

        int getID(double x, double y);
        void get_pixel_center(int ID, double& x, double& y);
        void initialize_geometry(std::string geom_file);
    }
}
