#include <Will.h>
#include "scipp_ilc_globals.h"
#include <sstream>
#include <iostream>
#include <TFile.h>
#include <TH1F.h>

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
fourvec fourvec::operator*(const double a) const{
  return fourvec(a*x,a*y,a*z);
}

fourvec fourvec::operator*=(const double a){
  this->x=this->x*a;
  this->y=this->y*a;
  this->z=this->z*a;
  this->e=this->e*a;
  this->t=this->t*a;
  return *this;
}
fourvec fourvec::operator*=(const fourvec& a){
  this->x=this->x*a.x;
  this->y=this->y*a.y;
  this->z=this->z*a.z;
  this->e=this->e*a.e;
  this->t=this->t*a.t;
  return *this;
}
fourvec fourvec::operator/=(const double a){
  this->x=this->x/a;
  this->y=this->y/a;
  this->z=this->z/a;
  this->e=this->e/a;
  this->t=this->t/a;
  return *this;
}
fourvec fourvec::operator/(const double a) const{
  return fourvec(a/(*this).x,a/(*this).y,a/(*this).z);
}

fourvec fourvec::operator+=(const fourvec& a){
  *this = *this+a;
  return *this;
}
/*string Will::str(fourvec in){
  ostringstream strs;
  strs << in.x;
  string x_str = strs.str();
  strs.str(string());
  strs << in.y;
  string y_str=strs.str();
  strs.str(string());
  strs << in.z;
  string z_str=strs.str();
  return "x:y:z " + x_str +":"+ y_str +":"+ z_str + " ";
  }*/


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
prediction::prediction(measure input){
  //const double ENERGY=input.electron.e+input.positron.e+input.hadronic.e;
  const double ENERGY=500;
  alpha=(ENERGY-input.hadronic.E - input.hadronic.z);
  beta=(ENERGY-input.hadronic.E + input.hadronic.z);
  electron.x = -input.hadronic.x;
  electron.y = -input.hadronic.y;
  electron.z = -(pow(input.electron.T, 2)-pow(alpha, 2))/(2*alpha);
  electron.e = ENERGY - input.hadronic.e - electron.z - input.hadronic.z;
  electron.t = getTMag(electron);

  positron.x = -input.hadronic.x;
  positron.y = -input.hadronic.y;
  positron.z = (pow(input.positron.T, 2)-pow(beta, 2))/(2*beta);
  positron.e = ENERGY - input.hadronic.e + positron.z + input.hadronic.z;
  positron.t = getTMag(positron);
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
      }
    }
  }
  return *max;
}

map<int,MCParticle*> Will::maxParticle(LCCollection* col, initializer_list<int> ids){
  const int SIZE=ids.size();
  map<int, MCParticle*>max;
  //map<int, int> max_id;
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
  input.T = sqrt(pow(input.x, 2) + pow(input.y, 2));
  return input.T;
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
  return a*b;
}

double Will::getTheta(const double* a, const double* b){
  return getTheta(fourvec(a,3),fourvec(b,3));
}
double Will::getTheta(fourvec a, fourvec b){
  //If acos(-1) it returns 0;
  double theta=a*b/(getMag(a)*getMag(b));
  double val=acos(theta);
  if(val!=val) return 0;
  return val;
}

measure Will::getMeasure(LCCollection* col){
  measure out;
  vector<MCParticle*> all_particles;
  vector<MCParticle*> hadronic_system;

  map<int, double> max=maxEnergy(col, {11, -11}, all_particles);
  /* Max energy is a big function. It finds the max energy of the given particle ids.
   * The sketchy thing I did was have it return a vector with final state particles (genStat==1)
   * That is what the all_particles thing is.
   */
  int hits = 0;
  out.mag=0;
  //Checks for scatter in electron or positron.
  for(auto particle: all_particles){
    if(particle->getGeneratorStatus()!=1) continue; //Redundant, maxEnergy does this.
    int id = particle->getPDG();
    if(particle->getEnergy()==max[11]){
      ++hits;
      //ELECTRON
      out.electron=getFourVector(particle);
      if(out.electron.T!=0.0){
	out.e_scatter=out.scattered=true;
	out.electronic=out.electron;
      }
    }else if(particle->getEnergy()==max[-11]){
      ++hits;
      //POSITRON
      out.positron=getFourVector(particle);
      if(out.positron.T!=0.0){
	out.p_scatter=out.scattered=true;
	out.electronic=out.positron;
      }    
    }else{
      //HADRONIC
      fourvec hadron=getFourVector(particle);
      hadronic_system.push_back(particle);
      out.hadronic+=hadron;
      out.mag += getTMag(hadron);
    }
  }
  out.hadronic_nopseudo=out.hadronic;
 
  //PSEUDO PARTICLE
  out.pseudo=*(new fourvec(
			   -(out.hadronic.x+out.positron.x+out.electron.y),
			   -(out.hadronic.y+out.positron.y+out.electron.y) ));
  out.hadronic.x=0;
  out.hadronic.y=0;
  double total_energy=out.hadronic.E;
  out.hadronic.z=0;
  out.hadronic.e=0;
  out.hadronic.t=0;
  int count = 0;
  short n_hadrons = hadronic_system.size();
  out.mag += getTMag(out.pseudo);    
  for(auto particle: hadronic_system){
    fourvec hadron = getFourVector(particle);
    ++count;
    out.hadronic += (hadron+=out.pseudo*(hadron.E/total_energy));
    //int id=particle->getPDG();
    //const double angle = 0.5; //Cutting angle
    //if(id==12||id==14||id==16||id==18|| (id>=1000001 && id<=1000039))continue;
    //      if(true || abs(hadron.z/getMag(hadron)) < angle ){
    //out.hadronic += hadron;
    //}
  }
  //out.hadronic+=out.pseudo;
 
  out.scattered ? meta.SCATTERS++ : meta.NOSCATTERS++; //Accounting for later statistical use.
  return out;
}

fourvec Will::transform_to_lab(fourvec input){
  double px=input.x;
  double e = input.e;
  scipp_ilc::transform_to_lab(px,e,px,e);
  fourvec out(px,input.y,input.z,e);
  return out;
}


fourvec Will::transform_to_lab(MCParticle* input){
  return transform_to_lab(getFourVector(input));
}

int Will::get_hitStatus(fourvec input, bool lab){
  if(lab) input=transform_to_lab(input);
  return scipp_ilc::get_hitStatus(input.x, input.y, input.z);
}

int Will::get_hitStatus(MCParticle* input){
  const double* mom = input->getMomentum();
  double x=mom[0],y=mom[1],z=mom[2];
  int stat = scipp_ilc::get_hitStatus(x,y,z);
  return stat;
}

fourvec Will::getBeamcalPosition(const fourvec input, signed short dir){
  fourvec lab = transform_to_lab(input);
  fourvec pos;
  //Positron moves in -z direction
  double direction = lab.z / abs(lab.z);
  if(dir != 0 && dir != direction) meta.err_direction++;
  pos.z = META::BEAMCAL * direction;
  pos.x = lab.x * pos.z / lab.z + pos.z * .007 * (-direction);
  pos.y = lab.y * pos.z / lab.z;
  return pos;
}

void Will::print(string input){
  cout << input << endl;
}

void Will::print(string input, string input2){
  cout << input << " : " << input2 << endl;
}

META Will::getMETA(){return meta;}

hmgrid Will::getHMGrid(vector<fourvec> predicted, vector<fourvec> actual){
  hmgrid output;
  //Create position vectors
  for(unsigned int i=0; i < predicted.size(); ++i){
    recordHMValue(output, predicted[i], actual[i]);
  }
  return output;
}

//Calculates a HM grid with the option of an energy cut.
hmgrid Will::getHMGrid(vector<Bundle> input, double energy_cut){
  hmgrid output;
  for(Bundle bundle: input){
    if(bundle.system_energy>=energy_cut)
      recordHMValue(output, bundle.predicted, bundle.actual);
  }
  return output;
}
//Sees if the predicted and actual vector hit the beamcal, records the results in a hmgrid object.
void Will::recordHMValue(hmgrid &output, fourvec predicted, fourvec actual){
  fourvec real=getBeamcalPosition(actual);
  fourvec pred=getBeamcalPosition(predicted);
  bool hit_real=get_hitStatus(real)<3;
  bool hit_pred=get_hitStatus(pred)<3;
  if     (  hit_pred &&  hit_real )output.hh++;
  else if( !hit_pred &&  hit_real )output.mh++;
  else if(  hit_pred && !hit_real )output.hm++;
  else if( !hit_pred && !hit_real )output.mm++;
}
void Will::printHMGrid(vector<fourvec> pred, vector<fourvec> actual){
  printHMGrid(getHMGrid(pred,actual));
}
void Will::printHMGrid(vector<Bundle> input, double energy_cut){
  printHMGrid(getHMGrid(input, energy_cut));
}

void Will::printHMGrid(hmgrid input){
  double sum=input.hh+input.hm+input.mh+input.mm;
  cout << "Total Events in HM Grid: " << sum << endl;
  cout << "          | Truth Hit | Truth Miss" << endl;
  cout << "Pred Hit  | " << input.hh/sum << "   |  " << input.hm/sum << endl;
  cout << "Pred Miss | " << input.mh/sum << "   |  " << input.mm/sum << endl;
}

TH1F* Will::getDistribution(string name,vector<Bundle> input, double energy_cut, double upper_bound){
  stringstream strs;
  strs << energy_cut;
  string cut = strs.str();
  string title="Theta Distribution Above "+cut;
  TH1F* output=new TH1F(name.c_str(),title.c_str(),300,0.0,upper_bound);
  for(Bundle bundle: input){
    if(bundle.system_energy > energy_cut)
      output->Fill(getTheta(bundle.actual, bundle.predicted));
  }
  return output;
}
