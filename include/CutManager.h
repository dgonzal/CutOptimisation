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
  CutManager(){ scan_permutations=1;}
  bool add_Cut(std::string nick, std::string cutname, std::vector<std::string> parameters);
  bool check_Cut(int cuti, double val);
  void permute_Cuts(int cuti);
  int number_Cuts(){return cutStore.size();}
  int get_NumberScanPoints(){return scan_permutations;}
  std::string get_OuputName();
 private:
  std::vector<Cut*> cutStore;
  std::vector<std::string> cutNicks;
  int scan_permutations;
  int find_position(int NumberOfSteps);
};

