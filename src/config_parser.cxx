#include "config_parser.h"

#include <iostream>

using namespace std;

configInfo::configInfo(const boost::property_tree::ptree & configTree){
  std::cout<<"Current Configuration:"<<std::endl;
  workdir = "./workdir";
  path = "";
  weight = "";
  for (auto& section : configTree){
    //std::cout << '[' << section.first << "]\n";
    for (auto& key : section.second){
      if(section.first == "CutScan"){
	cutInfo new_cut;
	new_cut.nick = key.first;
	string completeInfo =  key.second.get_value<std::string>();
	boost::split(new_cut.info,completeInfo,boost::is_any_of("\t "));
	//std::cout<<new_cut.info[0]<<std::endl;
	cuts.push_back(new_cut);
	}
      if(section.first == "Signal" || section.first == "Background"){
	sampleInfo new_sample;
	new_sample.nick = key.first;
	new_sample.dir = key.second.get_value<std::string>();
	if(section.first == "Signal")
	  new_sample.sampleType = sampleInfo::sampleTypes::signal;
	else if(section.first == "Background")
	  new_sample.sampleType = sampleInfo::sampleTypes::background;
	samples.push_back(new_sample);
      }
      if(section.first == "Setup"){
	if(key.first=="SamplePath") path =  key.second.get_value<std::string>();
	if(key.first=="Model") ThetaModelDir = key.second.get_value<std::string>();
	if(key.first=="Observable"){ 
	  std::vector<std::string> info;
	  string mykey = key.second.get_value<std::string>();
	  boost::split(info,mykey,boost::is_any_of("\t "));
	  observable = info[0];
	  if(info.size()>3){
	    observablePlotInfo.bins = stoi(info[1]);
	    observablePlotInfo.min = stod(info[2]);
	    observablePlotInfo.max = stod(info[3]);
	  }
	  else{
	    observablePlotInfo.bins = 0;
	    observablePlotInfo.min = 0;
	    observablePlotInfo.max = 0;
	  }
	}
	if(key.first=="Theta") ThetaDir = key.second.get_value<std::string>();
	if(key.first=="Workdir") workdir = key.second.get_value<std::string>();
	if(key.first=="Weight") weight = key.second.get_value<std::string>();
      }

      //std::cout << key.first << " = " << key.second.get_value<std::string>() << "\n";
    }
  }
  std::cout<<"Config read in finished"<<std::endl;
  std::cout<<"=========================="<<std::endl;
}

void configInfo::printConfig(){
  std::cout<<"Following Setup was read in:"<<endl;
  std::cout<<"Observable: "<<observable<<" Histogram Binning: "<<observablePlotInfo.bins<< " Min/Max "<<observablePlotInfo.min<<"/"<<observablePlotInfo.max<<endl;
  std::cout<<"Weight: "<<weight<<endl;
  std::cout<<"Workdir: "<<workdir<<endl;
  std::cout<<"Thera Dir: "<<ThetaDir<<endl;
  std::cout<<"Thera Model: "<<ThetaModelDir<<endl;
  for(auto& sample : samples)
    std::cout<<sample.nick<<" "<<sample.dir<<" Type "<<sample.sampleType<<std::endl;
  for(auto& cut : cuts)
    std::cout<<"Cut Nickname "<<cut.nick<<" Cut Name "<<cut.cutname()<<" Variable "<< cut.variable()<<std::endl;
  std::cout<<"================="<<std::endl;
}

std::vector<std::string> configInfo::get_SampleInformation(int option){
  vector<string> names;
  if(option==0)
    for(auto & sample : samples)
      names.push_back(sample.dir);
  else if(option==1)
    for(auto & sample : samples)
      names.push_back(sample.nick);
  else if(option==2)
    for(auto & sample : samples)
      if(sample.sampleType==0)names.push_back(sample.nick);
  return names;
}

std::vector<std::string> configInfo::get_Variables(){
  vector<string> var;
  string myobs = observable; string myweight = weight;
  var.push_back(myobs);
  var.push_back(weight);
  for(auto & cut : cuts){
    var.push_back(cut.variable());
  }
  return var;
}

