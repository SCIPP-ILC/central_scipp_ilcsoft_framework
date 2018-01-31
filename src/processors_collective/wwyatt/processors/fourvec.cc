#include <fourvec.h>
#include <cmath>
#include <Will.h>
fourvec fourvec::operator+(const fourvec& a) const{
  return fourvec(
		 a.x+x,
		 a.y+y,
		 a.z+z,
		 a.e+e);
}
double fourvec::operator*(const fourvec& a) const{
  return a.x*x + a.y*y + a.z*z;
}
fourvec fourvec::operator*(const double a) const{
  return fourvec(a*x,a*y,a*z);
}

fourvec fourvec::operator*=(const double a){
  this->x=this->x*a;
  this->y=this->y*a;
  this->z=this->z*a;
  this->e=this->e*a;
  return *this;
}
fourvec fourvec::operator*=(const fourvec& a){
  this->x=this->x*a.x;
  this->y=this->y*a.y;
  this->z=this->z*a.z;
  this->e=this->e*a.e;
  return *this;
}
fourvec fourvec::operator/=(const double a){
  this->x=this->x/a;
  this->y=this->y/a;
  this->z=this->z/a;
  this->e=this->e/a;
  return *this;
}
fourvec fourvec::operator/(const double a) const{
  return fourvec(a/(*this).x,a/(*this).y,a/(*this).z);
}

fourvec fourvec::operator+=(const fourvec& a){
  *this = *this+a;
  return *this;
}

fourvec::fourvec(const double _x,const double _y){x=_x;y=_y;}
fourvec::fourvec(const double _x,const double _y,const double _z):fourvec(_x,_y){z=_z;}
fourvec::fourvec(const double _x,const double _y,const double _z,const double _e):fourvec(_x,_y,_z){e=_e;}
fourvec::fourvec():fourvec(0,0,0,0){}
fourvec::fourvec(const double* input,const unsigned short SIZE){
  switch(SIZE){
  case 4:
    e=input[3];
  case 3:
    z=input[2];
  case 2:
    y=input[1];
  case 1:
    x=input[0];
  }
}
double fourvec::getTMag(){
  return sqrt(pow(this->x,2)+pow(this->y,2));
}
double fourvec::getMag(){
  return sqrt(pow(this->x,2)+pow(this->y,2)+pow(this->z,2));
}
