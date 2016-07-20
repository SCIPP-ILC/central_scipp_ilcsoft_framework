/**
 * PolarCoords.hpp
 *
 * Created on June 13
 *
 * @author Christopher Milke
 * @version 1.0 
 * 
 * This is a simple class designed to define polar coordinates 
 * for the use of the tile geometry. The coordinate "r" is (of course)
 * defined as 0 at the origin, and increases away from the origin as a
 * function of x and y ( r^2 = x^2 + y^2   as is standard). The coordinate
 * "phi" is defined as the counter-clockwise angle (in radians) with the 
 * x-axis as phi=0 (and 2pi,4pi,etc). 
 * 
 */

namespace scipp_ilc {
    extern const double _beam_out_angle;

    //converts x,y cartesian coordinates to r,phi polar coordinates
    void cartesian_to_polar(double x, double y, double& radius, double& phi);
    
    
    //converts r,phi polar coordinates to x,y cartesian coordinates
    void polar_to_cartesian(double radius, double phi, double& x, double& y);


    //sets the origin of the x axis to be the outgoing beampipe
    void z_to_beam_out(double& x, double& y, double& z);
}
