#pragma once

#include <vector>
#include <string>
#include <iostream>


// all samples used within one run of the program must be constructed with the same obs_selection.
struct obs_selection{
  enum calc_type{multiply,divide,add, substract};
  std::vector<std::string> observable_names;
  std::vector<std::string> added_names;//variables derivated from obs_names
  std::vector<operands_type> operands;
  std::vector<calc_type> operation;
  size_t indexof_obs(const string & observable_name) const {
    size_t index = std::find(observable_names.begin(), observable_names.end(), observable_name) - observable_names.begin();
    if(index == observable_names.size()){
      throw std::string("sample_info::indexof: did not find observable '" + observable_name + "'");
    }
    return index;
  }
  size_t indexof(const string & observable_name) const {
    size_t index = std::find(observable_names.begin(), observable_names.end(), observable_name) - observable_names.begin();
    size_t add_index = std::find(added_names.begin(), added_names.end(), observable_name) - added_names.begin();
    if(index == observable_names.size() && add_index == added_names.size()){
      throw std::string("sample_info::indexof: did not find observable '" + observable_name + "'");
    }
    if(add_index == added_names.size()) add_index =0;

    return index+add_index;
  }
};
