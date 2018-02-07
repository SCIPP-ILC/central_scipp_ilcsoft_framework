#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0

//#include <vector>
//#include "marlin/Processor.h"
//#include "marlin/VerbosityLevels.h"

#include "fastjet/ClusterSequence.hh"
//#include "lcio.h"
#include <iostream>
using namespace fastjet;
using namespace std;

int main () {
    vector<PseudoJet> particles;
    // an event with three particles: px py pz E
    particles.push_back( PseudoJet( 99.0, 0.1, 0, 100.0) );
    particles.push_back( PseudoJet( 4.0, -0.1, 0, 5.0) );
    particles.push_back( PseudoJet( -99.0, 0, 0, 99.0) );
    // choose a jet definition
    double R = 0.7;
    JetDefinition jet_def(antikt_algorithm, R);
    // run the clustering, extract the jets
    ClusterSequence cs(particles, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
    // print out some info
    cout << "Clustered with " << jet_def.description() << endl;
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

