#pragma once 

#include <boost/algorithm/string.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/filesystem.hpp>

#include <string>
#include <vector>

struct sampleInfo{
  std::string dir;
  std::string nick;
  enum sampleTypes{signal, background} sampleType;
};
struct cutInfo{
  std::vector<std::string> info;
  std::string nick;
  std::string cutname(){return info[0];}
  std::string variable(){return info[1];}
  std::vector<std::string> cutInput(){
    std::vector<std::string> myInput;
    for(unsigned int i=2;i<info.size();i++)
      myInput.push_back(info[i]);
    return myInput;
  }
};
struct obsInfo{
  std::string name;
  int bins;
  double min;
  double max;
};

class configInfo {
 public:
  configInfo(const boost::property_tree::ptree & configTree);
  enum optimisation_method{Theta};
  void printConfig();
  const std::string & get_Workdir(){return workdir;}
  const std::vector<obsInfo> & get_PlotInfo(){return observablePlotInfo;}
  const std::string & get_ThetaDir(){return ThetaDir;}
  const std::string & get_ThetaModelDir(){return ThetaModelDir;}
  std::vector<sampleInfo> get_Samples(){return samples;}
  std::vector<std::string> get_SampleInformation(int option =0);
  std::vector<std::string> get_Variables();
  std::vector<cutInfo> get_Cuts(){return cuts;}
  const std::string & get_Path(){return path;}
  const std::vector<std::string> & get_Sys(){return sys;}
  const std::string & get_Rebinning(){return rebinning_script;} 
  //optimisation_method get_Method(){return optiMethod}
 private:
  optimisation_method optiMethod;
  std::string path;
  std::string workdir;
  std::vector<obsInfo> observablePlotInfo;
  std::string weight;
  std::string ThetaDir;
  std::string ThetaModelDir;
  std::vector<sampleInfo> samples;
  std::vector<cutInfo> cuts;
  std::vector<std::string> sys; 
  std::string rebinning_script="";
}; 

