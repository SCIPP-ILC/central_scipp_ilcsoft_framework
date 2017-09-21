#include <Will.h>
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
  return getTMag(new double[2]{input.x,input.y});
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
	out.scattered=true;
	out.electronic=out.electron;
      }
    }else if(particle->getEnergy()==max[-11]){
      //POSITRON
      out.positron=getFourVector(particle);
      if(out.positron.T!=0.0){
	out.scattered=true;
	out.electronic=out.positron;
      }    
    }else{
      //HADRONIC
      fourvec hadron=getFourVector(particle);
      out.hadronic+=hadron;
      out.mag+=hadron.T;
    }
  }
  return out;
}


void Will::print(string input){
  cout << input << endl;
}

void Will::print(string input, string input2){
  cout << input << " : " << input2 << endl;
}
