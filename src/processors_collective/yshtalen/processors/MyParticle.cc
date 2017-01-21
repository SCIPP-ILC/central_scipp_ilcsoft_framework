//created by Jane Shtalenkova
//January  10, 2016
//extends LCIO MCParticle class


#include "MyParticle.h"

MyParticle::MyParticle(MCParticle* in_particle) : source{in_particle}{}

MyParticle::~MyParticle(){ 
    delete source; 
    delete momentum; 
}

//create editable momentum (non-constant)
MyParticle::mom(){
    const double* temp_mom = source->getMomentum();
    momentum[0] = temp_mom[0];    
    momentum[1] = temp_mom[1];    
    momentum[2] = temp_mom[2];
    momentum[3] = source->getEnergy();
    return momentum;    
}

MyParticle::energy(){
    en = source->getEnergy();
    return en;    
}

MyParticle::id(){
    return source->getPDG();    
}

MyParticle::stat(){
    return source->getGeneratorStatus();
}

//booleans
MyParticle::isHadronic(){return had;}

MyParticle::isElectronic(){return elec;}

MyParticle::isDetectable(){return able;}

MyParticle::isDetected(){return detect;}

//setters
MyParticle::setHadronic(set){had=set;}

MyParticle::setElectronic(set){elec=set;}

MyParticle::setDetectable(set){able=set;}

MyParticle::setDetected(set){detect=set;}

//transform momentum to lab frame from CM frame
MyParticle::getLorentzMom(){
    scipp_ilc::transform_to_lab(double momentum[0], double momentum[3], double& momentum[0], double& momentum[3])   
    return momentum;
}



