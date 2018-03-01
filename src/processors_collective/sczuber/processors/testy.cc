#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std; 

int main(){
    
    int jets[6] = {0,1,2,3,4,5};
    for (int i = 0; i<6; i++){
        cout << jets[i] << endl; 
    }
    //char counter[] = new char[6];
    for (unsigned int i = 0; i< pow(2,6); i++){
        cout << "NEW SUBSET   " << i << endl; 
        std::vector<int> subet;

        for (unsigned int j = 0; j<6; j++){
            unsigned int bit = pow(2,j); 
            //cout << "bit " << bit << endl; 
            if ((i  &  bit) != 0){
                subet.push_back(jets[j]);
            }
        }
        for (int k = 0; k<subet.size(); k++){
            cout << subet[k]; 
        }
        cout << endl; 
        subet.clear();
        
    }
    int subset[10] = {};
    return 0;


}
