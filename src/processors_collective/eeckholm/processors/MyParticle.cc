//created by Jane Shtalenkova
//January  10, 2016
//extends LCIO MCParticle class


#include "MyParticle.h"

//MyParticle::MyParticle(MCParticle* in_particle) : source{in_particle}{cout << "MyParticle object created" << endl;}
MyParticle::MyParticle(){cout << "MyParticle object created" << endl;}

/*MyParticle::~MyParticle(){ 
    cout << "MyParticle object destroyed" << endl;
    delete source; 
    delete momentum; 
}*/

void MyParticle::setSource(MCParticle* in_part){source=in_part;}

//create editable momentum (non-constant)
double* MyParticle::mom(){
    const double* temp_mom = source->getMomentum();
    momentum[0] = temp_mom[0];    
    momentum[1] = temp_mom[1];    
    momentum[2] = temp_mom[2];
    momentum[3] = source->getEnergy();
    return momentum;    
}

double MyParticle::energy(){
    en = source->getEnergy();
    return en;    
}

int MyParticle::id(){
    return source->getPDG();    
}

int MyParticle::stat(){
    return source->getGeneratorStatus();
}

//booleans
bool MyParticle::isHadronic(){return had;}

bool MyParticle::isElectronic(){return elec;}

bool MyParticle::isDetectable(){return able;}

bool MyParticle::isDetected(){return detect;}

//setters
void MyParticle::setHadronic(bool set){had=set;}

void MyParticle::setElectronic(bool set){elec=set;}

void MyParticle::setDetectable(bool set){able=set;}

void MyParticle::setDetected(bool set){detect=set;}

//transform momentum to lab frame from CM frame
double* MyParticle::getLorentzMom(){
    scipp_ilc::transform_to_lab(momentum[0], momentum[3], momentum[0], momentum[3]);   
    return momentum;
}



