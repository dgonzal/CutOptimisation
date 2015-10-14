#include "CutManager.h"


using namespace std;

bool CutManager::add_Cut(string nick, string cutname, vector<string> parameters){
  cutNicks.push_back(nick);
  if(cutname=="CutSmaller"){
    double step = stod(parameters[0]);
    double min = stod(parameters[1]);
    double max = stod(parameters[2]);
    double numberOfSteps;
    double frac = modf((max-min)/step, &numberOfSteps);
    if(frac!=0){
      numberOfSteps+=1;
      cout<<"number of steps is not an int, going one step beyond"<<endl;
    }   
    cutStore.push_back(new CutSmaller(step,min,max,int(numberOfSteps)));
    scan_permutations*=numberOfSteps;
    return true;
  }
  else if(cutname=="CutGreater"){
    double step = stod(parameters[0]);
    double min = stod(parameters[1]);
    double max = stod(parameters[2]);
    double numberOfSteps;
    double frac = modf((max-min)/step, &numberOfSteps);
    if(frac!=0){
      numberOfSteps+=1;
      cout<<"number of steps is not an int, going one step beyond"<<endl;
    }
    cutStore.push_back(new CutGreater(step,min,max,int(numberOfSteps)));
    scan_permutations*=numberOfSteps;
    return true;
  }
  else
    cout<<" cut seems not to be implemented, prefere to exit "<<endl;
  abort();
  return false;
}
bool CutManager::check_Cut(int cuti, double val){
  return cutStore[cuti]->passes(val);
}
void CutManager::permute_Cuts(int cuti){
  for(unsigned int p=0; p < cutStore.size(); p++){
    int result = cuti;
    for(unsigned int i=0; i <= p; i++){
      if(i<p) result/=cutStore[p]->get_numberOfSteps();
      if(p==i)result = int(result)%cutStore[p]->get_numberOfSteps();
    }
    //cout<< "Cut "<<p+1<< " Result: " <<result<<" ";
    cutStore[p]->set_counter(result);
  }
  //cout<<endl;
}
string CutManager::get_OuputName(){
  string result;
  for(unsigned int i = 0; i<cutNicks.size(); i++){
    //cout<<cutNicks[i] +"_"+ cutStore[i]->get_CutValues()<<endl;;
    if(i>0) result+= "_";
    result+= cutNicks[i] +"_"+ cutStore[i]->get_CutValues();
  }
  return result;
}