#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0

/*
 * Copy of short-example.cc for figuring out how to use FastJet
 */

#include <vector>
//#include "marlin/Processor.h"
//#include "marlin/VerbosityLevels.h"

#include "fastjet/ClusterSequence.hh"
//#include "lcio.h"
#include </cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-09/CLHEP/2.1.4.1/include/CLHEP/Vector/ThreeVector.h>
#include <iostream>
namespace CLHEP{}
using namespace CLHEP ;
//using namespace lcio ; 
using namespace fastjet;
using namespace std;

vector<PseudoJet> makeJets(vector<Hep3Vector> _parp){

    vector<PseudoJet> particles;
    // an event with three particles: px py pz E
    for (unsigned n=0; n<_parp.size(); n++){
        cout << "entered particle loop" << endl; 
        particles.push_back(   PseudoJet(_parp.at(n)[0],_parp.at(n)[1], _parp.at(n)[2], _parp.at(n)[3]));
        //particles.push_back(   PseudoJet(_parp.at(n)));
        //particles.push_back( PseudoJet( 99.0, 0.1, 0, 100.0) );
        //particles.push_back( PseudoJet( 4.0, -0.1, 0, 5.0) );
        //particles.push_back( PseudoJet( -99.0, 0, 0, 99.0) );
    }
    // choose a jet definition
    double R = 0.7;
    JetDefinition jet_def(antikt_algorithm, R);
    // run the clustering, extract the jets
    ClusterSequence cs(particles, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
    return jets;
}

vector<double> getJetEnergy(vector<Hep3Vector> _parp){
    vector<PseudoJet> jets = makeJets(_parp);
    vector<double> jetEnergy;
    for (unsigned i = 0; i<jets.size(); i++){
        cout << " energy of jet: " << jets[i].e();
        jetEnergy.push_back(jets[i].e());
    }
    return jetEnergy;
}

int main () {
    vector<Hep3Vector> _example;

    _example.push_back( Hep3Vector( 99.0, 0.1, 0) ); // do i need to give particle energy? 
    _example.push_back( Hep3Vector( 4.0, -0.1, 0) );
    _example.push_back( Hep3Vector( -99.0, 0, 0) );
    vector<PseudoJet> jets = makeJets(_example) ; 
    vector<double> e_vec = getJetEnergy(_example);
    cout << "size of e_vec"<< e_vec.size() << endl;     
    cout << "energies: "  << endl; 
    for (unsigned i = 0; i<e_vec.size() ; i++){
        cout << e_vec[i] << endl; ;
    }
    // print out some info
    //cout << "Clustered with " << jet_def.description() << endl;
    // print the jets
    cout << " pt y phi" << endl;
    for (unsigned i = 0; i < jets.size(); i++) {
        cout << "jet " << i << ": "<< jets[i].perp() << " "<< jets[i].rap() << " " << jets[i].phi() << endl;
        vector<PseudoJet> constituents = jets[i].constituents();
        for (unsigned j = 0; j < constituents.size(); j++) {
            cout << " constituent " << j << "’s pt: "<< constituents[j].perp() << endl;
        }
    }

    return 0;
}


