#ifndef FOURVEC
#define FOURVEC 1
/* ==== fourvec ====
 * I made fourvec to make all the code more readable.
 * To make a fourvec just pass it components
 * fourvec example(momentum[0], momentum[1], momentum[2], particle.getEnergy()); 
 * To get momentum:
 * example.x //return x momentum.
 * example.T //returns transverse momentum
 * Also it is case insensitive so:
 * example.y = example.Y
 */

struct fourvec{
  union{ double X; double x=0.0; };
  union{ double Y; double y=0.0; };
  union{ double Z; double z=0.0; };
  union{ double E; double e=0.0; };
  fourvec operator+(const fourvec& a) const; //vector addition
  double operator*(const fourvec& a) const; //dot product
  fourvec operator*(const double) const; //scalar multiplication
  fourvec operator+=(const fourvec& a); //vector addition
  fourvec operator*=(const double a); //scalar multiplication
  fourvec operator/=(const double a); //scalar division
  fourvec operator/(const double a)const; //scalar division

  fourvec operator*=(const fourvec& a); //WARINING this is not a normal operation.
  /* The way the operator *= works is it multiplies all of the components of the vectors together and returns a new vector.
   * Example: Let A=<1,2,3> and B=<4,5,6>; 
   * A*=B will return <4,10,18>;
   */

  fourvec();
  fourvec(const double x,const double y);
  fourvec(const double x,const double y,const double z);
  fourvec(const double x,const double y,const double z,const double e);

  fourvec(const double*,const unsigned short SIZE); // This is to initialize from an array.
  /* First parameter is the array, second parameter is the size of the array.
   * Example:
   * int a[3] = {1,2,3}; //Implicitly meaning x=1, y=2, z=3, e=0;
   * fourvec b(a, 3); //Returns <1,2,3,0> as a fourvec object.
   */
  double getTMag();
  double getMag();
};

#endif
