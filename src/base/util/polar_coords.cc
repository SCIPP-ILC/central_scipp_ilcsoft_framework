
/**
 * PolarCoords.cc
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

#define _USE_MATH_DEFINES
#include <cmath>
#include <stdlib.h>

#include <iostream>

#include "polar_coords.h"
#include "scipp_ilc_globals.h"

using namespace std;


namespace scipp_ilc {
    //converts x,y cartesian coordinates to r,phi polar coordinates
    void cartesian_to_polar(double x, double y, double& radius, double& phi) {
        radius = hypot(x,y);

        //There are so many else ifs here to handle the edge
        //cases where a point is on an axis
        if (y==0 && x>=0) phi = 0;
        else if (y==0 && x<0) phi = M_PI;
        
        else if (x==0 && y>0) phi = M_PI/2;
        else if (x==0 && y<0) phi = 3*M_PI/2;

        else if (x>0 && y>0) phi = atan(y/x);
        else if (x>0 && y<0) phi = atan(y/x) + 2*M_PI;
        else if (x<0 && y>0) phi = atan(y/x) + M_PI;
        else if (x<0 && y<0) phi = atan(y/x) + M_PI;
    
        else {
            cout << "Something has gone horribly wrong at polar_coords.cc" << endl;
            cout << "x = " << x << ", y = " << y << endl;
            exit(314);
        }    
    }
    
    
    //converts r,phi polar coordinates to x,y cartesian coordinates
    void polar_to_cartesian(double radius, double phi, double& x, double& y) {
        x = radius * cos(phi);
        y = radius * sin(phi);
    }


    //sets the origin of the x axis to be the outgoing beampipe
    void z_to_beam_out(double& x, double& y, double& z) {
        x = x - abs(z*_transform);
    }


}
