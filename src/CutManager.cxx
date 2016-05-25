#include "CutManager.h"

#include <math.h>       
#include <stdlib.h>     
using namespace std;

bool CutManager::add_Cut(string nick, string cutname, vector<string> parameters){
  double numberOfSteps;
  if(cutname=="CutSmaller" || cutname=="CutGreater" || cutname=="CutExact"){
    double step = stod(parameters[0]);
    double min = stod(parameters[1]);
    double max = stod(parameters[2]); 
    //cout<<"step "<<step<<" min "<<min<<" max "<<max<<endl;
    double frac = modf((max-min)/step, &numberOfSteps);
    if(frac!=0){
      numberOfSteps++;
      cout<<"number of steps is not an int, going one step beyond"<<endl;
    }   
    numberOfSteps++;

    int obs_num = -1;
    for(auto num : order_cuts)
      if(num>obs_num) obs_num =num;
    obs_num ++;
    order_cuts.insert(order_cuts.begin()+find_position(int(numberOfSteps)),obs_num);

    if(cutname=="CutSmaller"){
      cutStore.insert(cutStore.begin()+find_position(int(numberOfSteps)), new CutSmaller(step,min,max,int(numberOfSteps)));
    }
    else if(cutname=="CutGreater"){
      cutStore.insert(cutStore.begin()+find_position(int(numberOfSteps)),new CutGreater(step,min,max,int(numberOfSteps)));
      
    }
    else if(cutname=="CutExact"){
      cutStore.insert(cutStore.begin()+find_position(int(numberOfSteps)),new CutExact(step,min,max,int(numberOfSteps)));
    }
    //cout<< nick <<" "<<find_position(int(numberOfSteps))<<endl;
    cutNicks.insert(cutNicks.begin()+find_position(int(numberOfSteps)),nick);
  }
  else{
    cout<<"CutClass seems not to be implemented in CutManager, prefere to exit!"<<endl;
    abort();
  }

  scan_permutations*=numberOfSteps;
  return true;
}
bool CutManager::check_Cut(int cuti, double val){
  int reorder = 0;
  for(unsigned int i =0; i < order_cuts.size(); i++)
    if(cuti == order_cuts[i]) reorder= i;
  return cutStore[reorder]->passes(val);
}
void CutManager::permute_Cuts(int cuti){
  if(PermutVec.size()>2 || PermutVec.size()==1)
    cuti = PermutVec[cuti];
  else if(PermutVec.size()==2)
    cuti = PermutVec[0]<PermutVec[1] ? PermutVec[0]+cuti : PermutVec[1]+cuti; 
  for(unsigned int p=0; p < cutStore.size(); p++){
    int result = cuti;
    //cout<<result<< " ";
    for(unsigned int i=0; i <= p; i++){
      if(i<p) result/=cutStore[i]->get_numberOfSteps();
      if(p==i)result = result%cutStore[i]->get_numberOfSteps();
    }
    //cout<< "Cut "<<p+1<< " Result: " <<result<<" ";
    cutStore[p]->set_counter(int(result));
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


unsigned int CutManager::find_position(unsigned int NumberOfSteps){
  unsigned int result =0;
  for(unsigned int i=1; i<cutStore.size()+1; i++){
    if(cutStore.size()>i && cutStore[i]->get_numberOfSteps() == NumberOfSteps) result = i;
    if(cutStore.size()>i+1 && cutStore[i-1]->get_numberOfSteps() > NumberOfSteps && cutStore[i+1]->get_numberOfSteps() < NumberOfSteps) result = i;
    if(cutStore.size()<i+1 && cutStore[i-1]->get_numberOfSteps() > NumberOfSteps) result = i;
  }
  return result;
}

int CutManager::get_NumberScanPoints(){
  if(PermutVec.size() ==0)
    return scan_permutations;
  else if(PermutVec.size()==1 || PermutVec.size()>2)
    return PermutVec.size();
  else
    return (fabs(PermutVec[0]-PermutVec[1]));
}
