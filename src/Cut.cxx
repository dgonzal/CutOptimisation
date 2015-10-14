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
  return to_string(int(step*Cut::get_counter()+min));
}
CutSmaller::CutSmaller(double step_, double min_, double max_, int numberOfSteps_):Cut(numberOfSteps_){
  step=step_;
  min=min_;
  max=max_;
}
bool CutSmaller::passes(double eventValue){
  return step*Cut::get_counter()+min < eventValue;
}
const string CutSmaller::get_CutValues(){
  return to_string(int(step*Cut::get_counter()+min));
}
