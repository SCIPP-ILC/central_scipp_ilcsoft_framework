#include <Will.h>
#include "scipp_ilc_globals.h"
#include <sstream>
#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <fourvec.h>
using namespace Will;
using namespace TwoPhoton;









void Will::print(string input){
  cout << input << endl;
}

void Will::print(string input, string input2){
  cout << input << " : " << input2 << endl;
}


TH1F* Will::getDistribution(string name,vector<Result> input, double energy_cut, double upper_bound){
  stringstream strs;
  strs << energy_cut;
  string cut = strs.str();
  string title="Theta Distribution Above "+cut;
  TH1F* output=new TH1F(name.c_str(),title.c_str(),300,0.0,upper_bound);
  for(Result result: input){
    if(result.system_energy > energy_cut)
      output->Fill(getTheta(result.actual, result.predicted));
  }
  return output;
}

META Will::getMETA(){return meta;}
