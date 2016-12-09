#ifndef BHABAH_UTIL
#define BHABAH_UTIL

#include<string>


//Namespace for all my bahbah processing. NOTE: PUT INTO SEPERATE FILE.
namespace BB{
  //Setting PDG constants to particle names (DID NOT VERIFY)
  const static int PHOTON = 22;
  const static int POSITRON = 11;
  const static int ELECTRON = -11;

  //Setting errors static so I can spell check them easily.
  const static string ERROR_NOT_BAHBAH_PARTICLE = "Not a Bahbah particle User-Error: Trying to add a particle that is not accounted for.";
  const static string ERROR_MAX_PHOTONS = "Max Photon User-Error: Trying to add photon but it is already full.";
  const static string ERROR_ALREADY_PARTICLE = "Particle Already There User-Error: Trying to add a paticle but the space is not NULL.";

  //Making a bundle to store and process my event because the data is being given to me nicly.
  class Bundle{
  public:
    int count = 0;
    MCParticle* photonA = NULL;
    MCParticle* photonB = NULL;
    MCParticle* positron = NULL;
    MCParticle* electron = NULL;
    MCParticle* photons[2] = {NULL, NULL};

    Bundle(){}//No constructor
    virtual void initializ();

    //Returns the total number of particles added. MAX is 4 right now.
    int getCount(){
      return count;
    }

    //-- Physics Section --\\
    //returns the angle from the x-y axis
    double getTheta(MCParticle* _input){
      return atan(_input->getMomentum()[0]/_input->getMomentum()[1]);
    }
    //returns the angle from the z axis
    double getPhi(MCParticle* _input){
      //Get the norm in the x,y plane:
      const double m_x = (_input->getMomentum())[0];
      const double m_y = (_input->getMomentum())[1];
      const double m_z = _input->getMomentum()[2];
      double norm = sqrt(m_x*m_x + m_y*m_y);
      double phi = atan(norm/m_z);

      //The norm in the opposite side, and the z axis is the same side.
      return phi;
    }

    //-- CompSci Section --\\
    //Helper function, can be set to the namespace for general use. I could not find println(), so I made it.
    void err(string _input){
      cout << _input << endl;
    }

    //Checks to see if particle is already there, then adds it for processing. Once full it runs initialize.
    void addParticle(MCParticle* _input){
      switch(_input->getPDG()){
      case(PHOTON):
	addPhoton(_input);
	break;
      case(POSITRON):
	if(positron == NULL){
	  positron = _input;
	}else err(ERROR_ALREADY_PARTICLE);
	break;
      case(ELECTRON):
	if(electron == NULL){
	  electron = _input;
	}else err(ERROR_ALREADY_PARTICLE);
	break;
      default:
	err(ERROR_NOT_BAHBAH_PARTICLE);
      } 
      if(++count == 4){
	initialize();
      }
    }
    
    //Preparing for updating class for n photons. Currently max is 2 photons.
    void addPhoton(MCParticle* _input){
      if(photonA == NULL){
	photonA = _input;
	photons[0] = photonA;
      }else if(photonB ==NULL){
	photonB = _input;
	photons[1] = photonB;
      }else err(ERROR_MAX_PHOTONS);
    }

  private:
  };
}

#endif
