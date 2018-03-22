#include <TwoPhoton.h>
#include <iostream>
#include <iomanip>
using namespace TwoPhoton;
fourvec TwoPhoton::transform_to_lab(fourvec input){
  double px=input.x;
  double e=input.e;
  scipp_ilc::transform_to_lab(px,e,px,e);
  fourvec out(px,input.y,input.z,e);
  return out;
}


fourvec TwoPhoton::transform_to_lab(MCParticle* input){
  return transform_to_lab(getFourVector(input));
}

int TwoPhoton::get_hitStatus(fourvec input, bool lab){
  if(lab) input=transform_to_lab(input);
  return scipp_ilc::get_hitStatus(input.x, input.y, input.z);
}

int TwoPhoton::get_hitStatus(MCParticle* input){
  const double* mom = input->getMomentum();
  double x=mom[0],y=mom[1],z=mom[2];
  int stat = scipp_ilc::get_hitStatus(x,y,z);
  return stat;
}

prediction::prediction(bundle input){
  //const double ENERGY=input.electron.e+input.positron.e+input.hadronic.e;
  const double ENERGY=500;
  alpha=(ENERGY-input.hadronic.E - input.hadronic.z);
  beta=(ENERGY-input.hadronic.E + input.hadronic.z);
  electron.x = -input.hadronic.x;
  electron.y = -input.hadronic.y;
  electron.z = -(pow(getTMag(input.electron), 2)-pow(alpha, 2))/(2*alpha);
  electron.e = ENERGY - input.hadronic.e - electron.z - input.hadronic.z;

  positron.x = -input.hadronic.x;
  positron.y = -input.hadronic.y;
  positron.z = (pow(getTMag(input.positron), 2)-pow(beta, 2))/(2*beta);
  positron.e = ENERGY - input.hadronic.e + positron.z + input.hadronic.z;
}

map<int, double> TwoPhoton::maxEnergy(LCCollection* col, initializer_list<int> ids, vector<MCParticle*>& fs){
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

map<int,MCParticle*> TwoPhoton::maxParticle(LCCollection* col, initializer_list<int> ids){
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

bundle TwoPhoton::getHadronicSystem(LCCollection* col, bool ptKick){
  bundle out;
  vector<MCParticle*> all_particles;
  vector<MCParticle*> hadronic_system;

  map<int, double> max=maxEnergy(col, {11, -11}, all_particles);
  /* Max energy is a big function. It finds the max energy of the given particle ids.
   * The sketchy thing I did was have it return a vector with final state particles (genStat==1)
   * That is what the all_particles thing is.
   */
  //I am not proud of using the all_particles vector, but I did. It was efficient code, but confusing to read.

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
      if(getTMag(out.electron)!=0.0){
	out.e_scatter=out.scattered=true;
      }
    }else if(particle->getEnergy()==max[-11]){
      ++hits;
      //POSITRON
      out.positron=getFourVector(particle);
      if(getTMag(out.positron)!=0.0){
	out.p_scatter=out.scattered=true;
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
  if(ptKick){
  //PSEUDO PARTICLE
    out.pseudo=*(new fourvec(
			     -(out.hadronic.x+out.positron.x+out.electron.x),
			     -(out.hadronic.y+out.positron.y+out.electron.y) ));
    double total_energy=out.hadronic.E;
    out.hadronic.x=0;
    out.hadronic.y=0;
    out.hadronic.z=0;
    out.hadronic.e=0;
    int count = 0;
    short n_hadrons = hadronic_system.size();
    out.mag += getTMag(out.pseudo);    
    for(auto particle: hadronic_system){
      fourvec hadron = getFourVector(particle);
      ++count;
      out.hadronic += (hadron+=out.pseudo*(hadron.E/total_energy));
      //To remove unobservable particles.
      //int id=particle->getPDG();
      //const double angle = 0.5; //Cutting angle
      //if(id==12||id==14||id==16||id==18|| (id>=1000001 && id<=1000039))continue;
      //      if(abs(hadron.z/getMag(hadron)) < angle ){
      //out.hadronic += hadron;
      //}
    }
  }
  return out;
}

fourvec TwoPhoton::getBeamcalPosition(const fourvec input, signed short dir){
  fourvec lab = transform_to_lab(input);
  fourvec pos;
  //Positron moves in -z direction
  double direction = lab.z / abs(lab.z);
  pos.z = 3265 * direction;
  pos.x = lab.x * pos.z / lab.z + pos.z * .007 * (-direction);
  pos.y = lab.y * pos.z / lab.z;
  return pos;
}



hmgrid TwoPhoton::getHMGrid(vector<fourvec> predicted, vector<fourvec> actual){
  hmgrid output;
  //Create position vectors
  for(unsigned int i=0; i < predicted.size(); ++i){
    recordHMValue(output, predicted[i], actual[i]);
  }
  return output;
}

//Calculates a HM grid with the option of an energy cut.
hmgrid TwoPhoton::getHMGrid(vector<Result> input, double energy_cut){
  hmgrid output;
  for(Result result: input){
    if(result.system_energy>=energy_cut)
      recordHMValue(output, result.predicted, result.actual);
  }
  return output;
}
//Sees if the predicted and actual vector hit the beamcal, records the results in a hmgrid object.
void TwoPhoton::recordHMValue(hmgrid &output, fourvec predicted, fourvec actual){
  fourvec real=getBeamcalPosition(actual);
  fourvec pred=getBeamcalPosition(predicted);

  bool hit_real=get_hitStatus(real)<3;
  bool hit_pred=get_hitStatus(pred)<3;
  if     (  hit_pred &&  hit_real )output.hh++;
  else if( !hit_pred &&  hit_real )output.mh++;
  else if(  hit_pred && !hit_real )output.hm++;
  else if( !hit_pred && !hit_real )output.mm++;
}
string TwoPhoton::getCombinationFromIndex(unsigned int dec){
  string out="";
  for(int i = 0; i <= 3; i++){
    int val=1<<i;
    if((dec & val)>0)out="1"+out;
    else out="O"+out;
  }
  return out;
}
string TwoPhoton::getGoodOrBad(int index){
  if(index==2||index==8||index==10)
    return "!";
  else return "-";
}
void TwoPhoton::printGuessTable(vector<Result> positron, vector<Result> electron){ 
  int array[16]={0};
  for(unsigned int i = 0; i < positron.size(); ++i){
    fourvec preal=getBeamcalPosition(positron[i].actual);
    fourvec ppred=getBeamcalPosition(positron[i].predicted);
    fourvec ereal=getBeamcalPosition(electron[i].actual);
    fourvec epred=getBeamcalPosition(electron[i].predicted);
    
    unsigned int ehit_pred=(get_hitStatus(epred)<3) << 3;
    unsigned int ehit_real=(get_hitStatus(ereal)<3) << 2;
    unsigned int phit_pred=(get_hitStatus(ppred)<3) << 1;
    unsigned int phit_real=(get_hitStatus(preal)<3) << 0;
    
    unsigned int out=phit_real|phit_pred|ehit_real|ehit_pred;
    ++array[out];
    
  }
  cout << "Printing hit table " << endl;
  cout << "#e-e+ " << endl << "#PTPT\t\t[Quantity]" << endl;
  double total=positron.size();
  long unsigned int bad=0;
  for(unsigned int i = 0; i < 16; ++i){
    if(i==2||i==8||i==10){
      bad += array[i];
    }
    cout << setprecision(4) << fixed << getGoodOrBad(i) << getCombinationFromIndex(i) << ": " << 100*array[i]/total << "%\t[" << array[i] << "]" << endl;
  }
  cout << "total: " << positron.size() << endl;
  cout << "bad events(!): " << 100*bad/total <<"%" << endl;
}
void TwoPhoton::printHMGrid(vector<fourvec> pred, vector<fourvec> actual){
  printHMGrid(getHMGrid(pred,actual));
}
void TwoPhoton::printHMGrid(vector<Result> input, double energy_cut){
  printHMGrid(getHMGrid(input, energy_cut));
}

void TwoPhoton::printHMGrid(hmgrid input){
  double sum=input.hh+input.hm+input.mh+input.mm;
  cout << "Total Events in HM Grid: " << sum << endl;
  cout << "          | Truth Hit | Truth Miss" << endl;
  cout << "Pred Hit  | " << input.hh/sum << "   |  " << input.hm/sum << endl;
  cout << "Pred Miss | " << input.mh/sum << "   |  " << input.mm/sum << endl;
}
double* TwoPhoton::getVector(MCParticle* particle){
  double* output=new double[4];
  const double* mom=particle->getMomentum();
  output[0]=mom[0];
  output[1]=mom[1];
  output[2]=mom[2];
  output[3]=particle->getEnergy();
  return output;
}

fourvec TwoPhoton::getFourVector(MCParticle* particle){
  fourvec* output=new fourvec;
  const double* mom=particle->getMomentum();
  output->x=mom[0];
  output->y=mom[1];
  output->z=mom[2];
  output->E=particle->getEnergy();
  return *output;
}

  
double* TwoPhoton::addVector(double* a, double* b, const int SIZE){
  double* output=new double[SIZE];
  for(int i=0; i<SIZE; ++i) output[i]=a[i]+b[i];
  return output;
}


double TwoPhoton::getTMag(const double* input){
  return sqrt(pow(input[0], 2) + pow(input[1], 2));
}
double TwoPhoton::getMag(const double* input){
  return sqrt(pow(input[0], 2) + pow(input[1], 2) + pow(input[2],2));
}

double TwoPhoton::getTMag(fourvec input){
  return sqrt(pow(input.x, 2) + pow(input.y, 2));
}
double TwoPhoton::getMag(fourvec input){
  return getMag(new double[3]{input.x,input.y,input.z});
}

double TwoPhoton::getTheta(const double* input){
  return asin(getTMag(input)/getMag(input));
}
double TwoPhoton::getTheta(fourvec input){
  return asin(getTMag(input)/getMag(input));
}

double TwoPhoton::getTheta(const double* a, const double* b){
  return getTheta(fourvec(a,3),fourvec(b,3));
}
double TwoPhoton::getTheta(fourvec a, fourvec b){
  //If acos(-1) it returns 0;
  double theta=a*b/(getMag(a)*getMag(b));
  double val=acos(theta);
  if(val!=val) return 0;
  return val;
}

double TwoPhoton::getPhi(fourvec input){
  return atan(input.y/input.x);
}
