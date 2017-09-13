#include <Will.h>
using namespace Will;
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

double* Will::legacy(fourvec input){
  double *out = new double[4];
  out[0]=input.x;
  out[1]=input.y;
  out[2]=input.z;
  out[3]=input.E;
  return out;
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
double Will::getMag(const fourvec input){
  return sqrt(pow(input.x, 2) + pow(input.y, 2) + pow(input.z,2));
}


prediction Will::getPrediction(LCCollection* col){
  prediction out;
  vector<MCParticle*> hadronic_system;
  map<int, double> max=maxEnergy(col, {11, -11}, hadronic_system);
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


