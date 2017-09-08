#include <iostream>
#include <utility>
#include <unordered_set>
#include <cmath>
#include <algorithm>

#include "test_beamcal_scanner.h"
#include "test_beamcal_reconstructor.h"
#include "list_geometry.h"


using namespace std;


namespace scipp_ilc {
    namespace beamcal_recon_test {



        static bool compare_pair( pair<int,float> pair1, pair<int,float> pair2 ) {
            //Sort from greatest to least energy.
            //Highest energy is zeroth element in list
            return ( pair1.second > pair2.second );
        }



        /*
         * Identifies the 50 highest (background-average subtracted) energy
         * layer-compressed pixels. To do this, it first takes all the pixels
         * from the pixel map, and loads them into a vector. The vector is
         * then sorted by their (background-average subtracted) energy, and
         * the top 50 highest are loaded into the seed map.
         */
        static void create_seed_list (pixel_map* pixels,
                                        unordered_map<int,double>* average_map,
                                        unordered_map<int,double>* std_dev_map,
                                        unordered_map<int,float>* seed_list) {

            //Read all pixels in the pixel map, get their bgd-avg subtracted
            //energy, create an [ID,bgd-sub energy] pair, and load those
            //pairs into a vector.
            vector< pair<int, float> > sorted_pixels;
            for (auto pixel : *pixels) { 
                int ID = pixel.first;
                float energy = pixel.second;
                float bgd_subtracted_energy = energy - (*average_map)[ID];
                pair<int,float> bgd_subtracted_pair( ID, bgd_subtracted_energy );
                sorted_pixels.push_back(bgd_subtracted_pair);
            }
            sort(sorted_pixels.begin(), sorted_pixels.end(), compare_pair);

            //Load the 50 highest bgd-sub energy pixels into the seed list,
            //storing the pixels' significance as the hashmap's values.
            int count = 0;
            int maximum = 50;
            for (auto bgd_subtracted_pixel : sorted_pixels) {
                int ID = bgd_subtracted_pixel.first;
                float bgd_subtracted_energy = bgd_subtracted_pixel.second;
                float std_dev = (*std_dev_map)[ID];
                if (std_dev == -1.0) { continue; }

                float significance = 0.0;
                significance = bgd_subtracted_energy / std_dev;
                (*seed_list)[ID] = significance;

                count++;
                if (count >= maximum) { break; }
            }
        }



        /*
         * Calculates the significance of a cluster made up of the IDs in ID_list. Every new cluster
         * requires a new average and standard deviation be calculated, by finding the energy sum
         * of the cluster pixles in every pixel map in the database.
         */
        float get_significance(vector<int>* ID_list, pixel_map* pixels, float& energy, double& average_background) {

            //Calculate average and std_dev for cluster
            double total_background_energy = 0.0;
            double total_squared_background_energy = 0.0;
            int weight = 0;
            for ( const pixel_map* stored_map : *_database ) {

                double map_background_energy = 0.0;
                for ( int ID : *ID_list ) {
                    if ( stored_map->find(ID) != stored_map->end() ) {
                        map_background_energy += stored_map->at(ID);
                    }
                }
                
                if (map_background_energy != 0.0) weight++;
                total_background_energy += map_background_energy;
                total_squared_background_energy += map_background_energy*map_background_energy;
            }

            if ( weight == 0 ) return 0.0;

            average_background = total_background_energy / weight;
            double square_of_averages = average_background*average_background;
            double average_of_squares = total_squared_background_energy / weight;
            double standard_deviation = sqrt(average_of_squares - square_of_averages);

	    /* double root_mean_square   = sqrt(square_of-averages / ) */


            //calculate significance
            for ( int ID : *ID_list ) energy += (*pixels)[ID];
            if ( energy <= 0.0 ) return 0;

            float significance = (float) ( (energy-average_background) / standard_deviation );

            return significance;
        }


      float root_mean_square(double *v, int n, vector<int>* ID_list, pixel_map* pixels, float& energy, double& average_background)
      {
	int i;
	double sum = 0.0;
	for(i = 0; i < n; i++)
	  sum += v[i] * v[i];
	return sqrt(sum / n);
      }

      int main(void)
      {
	double v[] = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
	printf("%f\n", rms(v, sizeof(v)/sizeof(double)));
	return 0;
      }


        /*
         * This function accrues all pixels adjacent to the pixel designated by "ID".
         * Everytime it adds an adjacent pixel, it tests to see if the significance of
         * the new, larger, cluster is larger than the significance before the it added
         * the adjacent cluster. If it is not, the extra pixel is removed from the cluster,
         * and the search around "ID" pixel continues.
         *
         * This function can be made recursive by switching on a line of code (see below).
         * Without recursion, the function will search the pixels immediately adjacent to it.
         * With recursion turned on, it will search the pixels immediately adjacent to it, and
         * upon succesfully adding a new adjacent pixel to the cluster, it will then see if it
         * can add any pixels adjacent to the newly added pixel. And it will check the pixels
         * adjacent to the pixel adjacent to the "ID" pixels, and so on.
         */
        static float cluster_seeker(int ID, float significance, vector<int>* ID_list, pixel_map* pixels,
                                    unordered_set<int>* searched_IDs, float& energy, double& bgd) {

            surrounding_ids* surroundings = (*_pixel_graph)[ID];
            
            int maximum_pixels = 4;
            int current_pixels = 1;
            for ( int neighbor_ID : *(surroundings->list) ) {
                //searched_IDs->end() means that neighbor_ID was not found
                if ( searched_IDs->find(neighbor_ID) != searched_IDs->end() ) continue;

                ID_list->push_back(neighbor_ID);
                float temp_energy = 0.0;
                double temp_bgd = 0.0;
                float new_significance = get_significance(ID_list,pixels,temp_energy,temp_bgd);

                if (new_significance > significance) {
                    searched_IDs->emplace(neighbor_ID);
                    energy = temp_energy;
                    bgd = temp_bgd;
                    significance = new_significance;
                    //To enable recursive clustering:
                    //comment out the above line, and uncomment the below line 
                    //significance = cluster_seeker(neighbor_ID,new_significance,ID_list,pixels,searched_IDs,energy,bgd);
                    
                    maximum_pixels++;
                } else {
                    ID_list->pop_back();
                }

                if (current_pixels > maximum_pixels) { break; }
            }
            return significance;
        }



        /*
         * Identify the most "significant" ( (energy - bgd_average) / standard_deviation ) cluster.
         *
         * This is done by iterating over the seed pixels, clustering around those pixels
         * (if you enable clustering that is), and selecting the cluster with the highest
         * significance value.
         */
        static beamcal_cluster* most_significant_cluster (pixel_map* pixels, unordered_map<int,float>* seed_list,
                                                            unordered_map<int,double>* average_map) {

            vector<int>* chosen_cluster = NULL;
            float chosen_significance = 0.0;
            float chosen_energy = 0.0;
            double chosen_bgd = 0.0;

            for( auto seed : *seed_list ) {
                int ID = seed.first;
                float significance = seed.second;
                float energy = (*pixels)[ID];
                double bgd = (*average_map)[ID];

                unordered_set<int>* searched_IDs = new unordered_set<int>();
                searched_IDs->emplace(ID);

                vector<int>* ID_list = new vector<int>();
                ID_list->push_back(ID);
                

                //uncomment the below line to enable pixel clustering
                //significance = cluster_seeker(ID,significance,ID_list,pixels,searched_IDs,energy,bgd);

                //choose the most significant cluster
                if ( significance > chosen_significance ) {
                    delete chosen_cluster;

                    chosen_significance = significance;
                    chosen_cluster = ID_list;
                    chosen_energy = energy;
                    chosen_bgd = bgd;

                } else {
                    delete ID_list;
                }
                delete searched_IDs;
            }

            
            //create cluster object, which contains useful information for
            //reconstructing the event later.
            beamcal_cluster* new_cluster;
            new_cluster = (beamcal_cluster*) malloc( sizeof(beamcal_cluster) );

            new_cluster->id_list = chosen_cluster;
            new_cluster->significance = chosen_significance;
            new_cluster->energy = chosen_energy;
            new_cluster->background_average = chosen_bgd;

            return new_cluster;
        }



        /*
         * Tries to identify the location of a signal event on the beamcal.
         * It requires two hashmaps: an one the averages for every pixel,
         * and one with the standard deviations for every pixel.
         * The scanning process is done in two steps: First, it identifies a
         * number of "seed" pixels using a simple algorithm. Second, it uses
         * a more rigorous clustering algorithm on the chosen seed pixels.
         * This second algorithm will determine if a signal event is present,
         * and return its location.
         */
        beamcal_cluster* scan_beamcal(pixel_map* pixels, unordered_map<int,double>* average_map, 
                                    unordered_map<int,double>* std_dev_map) {

            //Step 1: identify seed pixels
            unordered_map<int,float>* seed_list = new pixel_map();
            create_seed_list(pixels,average_map,std_dev_map,seed_list);

            //Step 2: use more advanced clustering algorithm to find
            //signal event amid seed pixels.
            return most_significant_cluster(pixels,seed_list,average_map);
        }
    }
}
