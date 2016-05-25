#pragma once 

#include "Cut.h"

#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <stdlib.h>

/*
**
** Cut Factory Class
** In this class the strings from the config are transformed into c++ classes
**
**
*/

class CutManager {
 public:
  CutManager(){scan_permutations=1;}
  bool add_Cut(std::string nick, std::string cutname, std::vector<std::string> parameters);
  bool check_Cut(int cuti, double val);
  void set_PermutationVector(std::vector<int> & PermutVec_){PermutVec=PermutVec_;}
  void permute_Cuts(int cuti);
  int number_Cuts(){return cutStore.size();}
  int get_NumberScanPoints();
  std::vector<int> & get_PermutationVector(){return PermutVec;}
  std::string get_OuputName();
 private:
  std::vector<int> PermutVec;
  std::vector<Cut*> cutStore;
  std::vector<int> order_cuts;
  std::vector<std::string> cutNicks;
  unsigned int scan_permutations;
  unsigned int find_position(unsigned int NumberOfSteps);
};

