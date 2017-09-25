#include <Will.h>
#include "scipp_ilc_globals.h"
using namespace Will;
fourvec fourvec::operator+(const fourvec& a) const{
  return fourvec(
		 a.x+x,
		 a.y+y,
		 a.z+z,
		 a.e+e,
		 sqrt(pow(a.x+x,2)+pow(a.y+y,2))
		 );
}
double fourvec::operator*(const fourvec& a) const{
  return a.x*x + a.y*y + a.z*z;
}
fourvec fourvec::operator+=(const fourvec& a){
  *this=a+*this;
  return *this;
}
fourvec::fourvec(const double _x,const double _y){x=_x;y=_y;t=getTMag(new double[2]{_x,_y});}
fourvec::fourvec(const double _x,const double _y,const double _z):fourvec(_x,_y){z=_z;}
fourvec::fourvec(const double _x,const double _y,const double _z,const double _e):fourvec(_x,_y,_z){e=_e;}
fourvec::fourvec(const double _x,const double _y,const double _z,const double _e,const double _t):fourvec(_x,_y,_z,_e){t=_t;}
fourvec::fourvec():fourvec(0,0,0,0,0){}
fourvec::fourvec(const double* input,const unsigned short SIZE){
  switch(SIZE){
  case 4:
    e=input[3];
  case 3:
    z=input[2];
  case 2:
    t=getTMag(new double[2]{input[0], input[1]});
    y=input[1];
  case 1:
    x=input[0];
  }
}

prediction::prediction(double x,double y){
  electron.x=x;
  positron.y=y;
}



map<int, double> Will::maxEnergy(LCCollection* col, initializer_list<int> ids, vector<MCParticle*>& fs){
  const int SIZE=ids.size();
  //double* max=new double[SIZE];
  map<int, double>* max=new map<int,double>;
  //  map<int, int> max_id;
  for(int i=0; i < col->getNumberOfElements(); ++i){
    MCParticle* particle=dynamic_cast<MCParticle*>(col->getElementAt(i));
    if(particle->getGeneratorStatus()!=1)continue;
    fs.push_back(particle);
    int pid=particle->getPDG();
    for(auto id:ids){
      if(pid==id && particle->getEnergy() > (*max)[id]){
	(*max)[id]=particle->getEnergy();
	//	max_id[id]=i;
      }
    }
  }
  return *max;
}
map<int,MCParticle*> Will::maxParticle(LCCollection* col, initializer_list<int> ids){
  const int SIZE=ids.size();
  map<int, MCParticle*>max;
  //  map<int, int> max_id;
  for(int i=0; i < col->getNumberOfElements(); ++i){
    MCParticle* particle=dynamic_cast<MCParticle*>(col->getElementAt(i));
    if(particle->getGeneratorStatus()!=1)continue;
    int pid=particle->getPDG();
    for(auto id:ids){
      if(pid!=id)continue;
      double energy=0.0;
      if(max[id]==NULL) max[id]=particle;
      else energy=max[id]->getEnergy();
      if(particle->getEnergy()>energy) max[id]=particle;
    }
  }
  return max;
}

double* Will::getVector(MCParticle* particle){
  double* output=new double[4];
  const double* mom=particle->getMomentum();
  output[0]=mom[0];
  output[1]=mom[1];
  output[2]=mom[2];
  output[3]=particle->getEnergy();
  return output;
}

fourvec Will::getFourVector(MCParticle* particle){
  fourvec* output=new fourvec;
  const double* mom=particle->getMomentum();
  output->x=mom[0];
  output->y=mom[1];
  output->z=mom[2];
  output->E=particle->getEnergy();
  output->T=getTMag(mom);
  return *output;
}

  
double* Will::addVector(double* a, double* b, const int SIZE){
  double* output=new double[SIZE];
  for(int i=0; i<SIZE; ++i) output[i]=a[i]+b[i];
  return output;
}


double Will::getTMag(const double* input){
  return sqrt(pow(input[0], 2) + pow(input[1], 2));
}
double Will::getMag(const double* input){
  return sqrt(pow(input[0], 2) + pow(input[1], 2) + pow(input[2],2));
}

double Will::getTMag(fourvec input){
  return sqrt(pow(input.x, 2) + pow(input.y, 2));
}
double Will::getMag(fourvec input){
  return getMag(new double[3]{input.x,input.y,input.z});
}

double Will::getTheta(const double* input){
  return asin(getTMag(input)/getMag(input));
}
double Will::getTheta(fourvec input){
  return asin(getTMag(input)/getMag(input));
}

double Will::getDot(const double* a, const double* b){
  return getDot(fourvec(a,3),fourvec(b,3));
}
double Will::getDot(fourvec a, fourvec b){
  return a.x*b.x+a.y*b.y+a.z*b.z;
}

double Will::getTheta(const double* a, const double* b){
  return getTheta(fourvec(a,3),fourvec(b,3));
}
double Will::getTheta(fourvec a, fourvec b){
  return acos( a*b / (getMag(a)*getMag(b)) );
}


measure Will::getMeasure(LCCollection* col){
  measure out;
  vector<MCParticle*> hadronic_system;
  map<int, double> max=maxEnergy(col, {11, -11}, hadronic_system);
  double mag=0.0;
  //Checks for scatter in electron or positron.
  for(auto particle: hadronic_system){
    int id = particle->getPDG();
    if(particle->getEnergy()==max[11]){
      //ELECTRON
      out.electron=getFourVector(particle);
      if(out.electron.T!=0.0){
	out.e_scatter=out.scattered=true;
	out.electronic=out.electron;
      }
    }else if(particle->getEnergy()==max[-11]){
      //POSITRON
      out.positron=getFourVector(particle);
      if(out.positron.T!=0.0){
	out.p_scatter=out.scattered=true;
	out.electronic=out.positron;
      }    
    }else{
      //HADRONIC
      fourvec hadron=getFourVector(particle);
      out.hadronic+=hadron;
      out.mag+=hadron.T;
    }
    get_hitStatus(particle);
  }
  out.scattered ? meta.SCATTERS++ : meta.NOSCATTERS++;
  return out;
}

fourvec Will::transform_to_lab(fourvec input){
  //  cout << "start: " << input.x << endl;
  double theta = 0.007;
  double beta = sin(theta);
  double gamma = pow((1-pow(beta, 2)), -0.5);
  double e_temp = input.E;
  double E_new = input.E*gamma + gamma*beta*input.x;
  double pX_new = input.x*gamma + gamma*beta*e_temp;
  input.x=pX_new;
  input.e=E_new;
  input.t=getTMag(input);
  fourvec out(pX_new,input.y,input.z,E_new);
  //  cout << "end: " << pX_new << endl;
  return out;
}

fourvec Will::transform_to_lab(MCParticle* input){
  return transform_to_lab(getFourVector(input));
}

int Will::get_hitStatus(fourvec input, bool lab){
  if(lab){
    input=transform_to_lab(input);
    meta.MSC++;
    //    cout << "POST: " << input.x << endl;
  }
  double x=input.x, y=input.y, z=input.z;
  double rad = sqrt(pow(x, 2)+pow(y, 2));
  double x_shift = x - scipp_ilc::_BeamCal_zmin*tan(scipp_ilc::_crossing_angle);
  double rad_shift = sqrt(pow(x_shift, 2)+pow(y, 2));
  //cout << "T: " << input.Z << endl;
  //cout << "Rad: " << rad << endl;
  int stat = scipp_ilc::get_hitStatus(input.x, input.y, input.z);
  meta.STATS[stat]++;
  //  cout << " 
  //cout << "s: " <<STATS[1] << "-"<<STATS[2]<<"-"<<STATS[3]<<"-"<<STATS[4]<<endl;
  return stat;
}

int Will::get_hitStatus(MCParticle* input){
  const double* mom = input->getMomentum();
  double x=mom[0],y=mom[1],z=mom[2];
  int stat = scipp_ilc::get_hitStatus(x,y,z);
  return stat;
}

void Will::print(string input){
  cout << input << endl;
}

void Will::print(string input, string input2){
  cout << input << " : " << input2 << endl;
}

META Will::getMETA(){return meta;}
