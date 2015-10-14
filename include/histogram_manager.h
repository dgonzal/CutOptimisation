#pragma once 

#include "TH1D.h"

#include <string>

class HistoManager {
 public: 
  HistoManager();
  void write_HistsToFile(const std::string & filename);
  void fill_ithHist(int i, double val, double weight);
  void reset_Hists();
  std::vector<TH1D*> hists;
};
