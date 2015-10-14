#pragma once

#include <string>
#include <vector>
#include <iostream>

#include <boost/shared_array.hpp>

class sample{
 public:
  sample(const std::string & fname_pattern, const std::vector<std::string> & observables_, const std::string & TreeName= "AnalysisTree");
  int get_NumberEvents(){return n_events;}
  double get_Value(int event, int observable_pos){return observablesValues[observable_pos][event];}
  const std::string get_Name(){return sample_name;}
 private:
  size_t n_events;
  std::vector<std::string> observables; 
  std::vector<std::vector<double>> observablesValues;
  std::string sample_name;
};
