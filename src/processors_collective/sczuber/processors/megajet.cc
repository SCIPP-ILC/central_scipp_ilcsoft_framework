#include <iostream>
//#include <ostream>
//#include <vector>
#include <fstream>
#include <cmath>
//#include "fastjet/ClusterSequence.hh"
#include "megajet.h"

using namespace std; 

/*int main(){
    
    int jets[6] = {0,1,2,3,4,5};
    for (int i = 0; i<6; i++){
        cout << jets[i] << endl; 
    }

    int subset[10] = {};
    return 0;


}
*/

vector<PseudoJet> findMegaJets(vector<PseudoJet> jets){
  
  
    //ostream& operator<<(ostream&);  
    for (unsigned int i = 0; i< pow(2,6); i++){
        cout << "NEW SUBSET   " << i << endl; 
        std::vector<PseudoJet> subset;
        for (unsigned int j = 0; j<6; j++){
            unsigned int bit = pow(2,j);  
            if ((i  &  bit) != 0){
                subset.push_back(jets[j]);
            }
        }
        for (unsigned int k = 0; k<subset.size(); k++){
            cout << subset[k]; 
        }
        cout << endl;
         
        subset.clear();
        
    }
    vector<PseudoJet> megajets; 
    return megajets;
}
