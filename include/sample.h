#pragma once

#include <string>
#include <vector>
#include <iostream>

#include <boost/shared_array.hpp>

class sample{
 public:
  sample(const std::string & fname_pattern, const std::vector<std::string> & observables_, const std::string & TreeName= "AnalysisTree", const std::vector<std::string> & sys=std::vector<std::string>());
  int get_NumberEvents(){return n_events;}
  double get_Value(int event, int observable_pos){return observablesValues[event][observable_pos];}
  std::vector<double> get_SysWeight(int event);
  const std::string get_Name(){return sample_name;}
 private:
  size_t n_events;
  std::vector<std::string> observables; 
  std::vector<std::string> sys; 
  std::vector<std::vector<double>> observablesValues;
  std::vector<std::vector<double>> systematicsWeights;
  std::string sample_name;
  std::vector<std::string> operators_used = {"-","+","*","/"};
};
