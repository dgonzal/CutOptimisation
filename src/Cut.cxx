# include "Cut.h"

using namespace std;


Cut::Cut(int numberOfSteps_):numberOfSteps(numberOfSteps_){
  counter=0;
}
bool Cut::passes(double eventValue){
  cout<<" calling the interface function passes() , something went really wrong"<<endl;
  std::abort();
  return false;
}
const string Cut::get_CutValues(){
  cout<<" calling the interface function get_Cut(), something went really wrong"<<endl;
  std::abort();
  return " ";
}

CutGreater::CutGreater(double step_, double min_, double max_, int numberOfSteps_):Cut(numberOfSteps_){
  step=step_;
  min=min_;
  max=max_;
}
bool CutGreater::passes(double eventValue){
  return step*Cut::get_counter()+min > eventValue;
}
const string CutGreater::get_CutValues(){
  if((step*Cut::get_counter()+min) > max) 
    return  to_string(int(step*Cut::get_counter()+min))+" bigger then MAX value";
  return to_string(int(step*Cut::get_counter()+min));
}
CutSmaller::CutSmaller(double step_, double min_, double max_, int numberOfSteps_):Cut(numberOfSteps_){
  step=step_;
  min=min_;
  max=max_;
}
bool CutSmaller::passes(double eventValue){
  //cout<<"cut val "<<step*Cut::get_counter()+min<<endl;
  return step*Cut::get_counter()+min < eventValue;
}
const string CutSmaller::get_CutValues(){
  if((step*Cut::get_counter()+min) > max) 
    return  to_string(int(step*Cut::get_counter()+min))+" bigger then MAX value";
  return to_string(int(step*Cut::get_counter()+min));
}

CutExact::CutExact(double step_, double min_, double max_, double numberOfSteps_):Cut(numberOfSteps_){
  step=step_;
  min=min_;
  max=max_;
}
bool CutExact::passes(double eventValue){
  return step*Cut::get_counter()+min == eventValue;
}
const string CutExact::get_CutValues(){
  if((step*Cut::get_counter()+min) > max) 
    return  to_string(int(step*Cut::get_counter()+min))+" bigger then MAX value";
  return to_string(int(step*Cut::get_counter()+min));
}
