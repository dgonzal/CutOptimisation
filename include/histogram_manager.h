#pragma once 

#include "TH1D.h"

#include <string>
#include <vector>

#include "config_parser.h"

class HistoManager {
 public: 
  HistoManager();
  void create_Hists(const obsInfo & PlotInfo,const std::string & sampleName,const std::vector<std::string> & SysNames);
  void write_HistsToFile(const std::string & filename);
  void fill_ithHist(int i, double val, double weight, std::vector<double> sys_weights);
  void reset_Hists();
  std::vector<TH1D*> hists;
  std::vector<std::vector<TH1D*>> sys;
};
