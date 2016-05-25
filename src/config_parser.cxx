#include "config_parser.h"

#include <iostream>

using namespace std;

configInfo::configInfo(const boost::property_tree::ptree & configTree){
  auto cwd  = boost::filesystem::current_path();
  //std::cout<<"Current Configuration:"<<std::endl;
  //cout<<cwd.string()<<endl;
  workdir = "workdir";
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
	if(key.first=="Model"){
	  ThetaModelDir = key.second.get_value<std::string>();
	  if(!boost::filesystem::exists(ThetaModelDir)){ 
	    ThetaModelDir = cwd.string()+"/"+ThetaModelDir;
	    /*//This does not work on the bird cluster, why not ????
	    if(boost::filesystem::exists(cwd.string()+"/"+ThetaModelDir))
	      
	    else{
	      cerr<<"ThetaModelDir not found "<<ThetaModelDir<<endl;
	      cerr<<"ThetaModelDir not found "<<cwd.string()+"/"+ThetaModelDir<<endl;
	      abort();
	      }*/
	     //cout<<ThetaModelDir<<endl;
	  }
	}
	if(key.first=="RebinningScript"){
	  rebinning_script = key.second.get_value<std::string>();
	  if(!boost::filesystem::exists(rebinning_script)){ 
	    rebinning_script = cwd.string()+"/"+rebinning_script;
	    /*//This does not work on the bird cluster, why not ????
	    if(boost::filesystem::exists(cwd.string()+"/"+rebinning_script))
	      rebinning_script = cwd.string()+"/"+rebinning_script;
	    else{
	      cerr<<"RebinningScript not found "<<rebinning_script<<endl;
	      abort();
	      }*/
	  }
	}
	if(key.first=="Observable"){ 
	  std::vector<std::string> info;
	  string mykey = key.second.get_value<std::string>();
	  boost::split(info,mykey,boost::is_any_of("\t "));
	  if(info.size()%4!=0){
	    cout<<"wrong number of arguments provided for the plots"<<endl;
	    assert(1==0);
	  }
	  for(unsigned int i=0; i<info.size(); i = i+4){
	    obsInfo myStorageInfo;
	    myStorageInfo.name = info[i];
	    myStorageInfo.bins = stoi(info[i+1]);
	    myStorageInfo.min = stod(info[i+2]);
	    myStorageInfo.max = stod(info[i+3]);
	    observablePlotInfo.push_back(myStorageInfo);
	  }
	}
	if(key.first=="Theta") ThetaDir = key.second.get_value<std::string>();
	if(key.first=="Workdir") workdir = key.second.get_value<std::string>();
	if(key.first=="Weight") weight = key.second.get_value<std::string>();
	if(key.first=="Systematics"){
	  string sysinfo = key.second.get_value<std::string>();
	  boost::split(sys,sysinfo,boost::is_any_of("\t "));
	}
      }
      //std::cout << key.first << " = " << key.second.get_value<std::string>() << "\n";
    }
  }
  if(boost::contains(ThetaModelDir,workdir)){
    ThetaModelDir.erase(ThetaModelDir.find(workdir),workdir.size());
  }
   if(boost::contains(rebinning_script,workdir)){
    rebinning_script.erase(rebinning_script.find(workdir),workdir.size());
  }
  

  std::cout<<"Config read in finished"<<std::endl;
  std::cout<<"=========================="<<std::endl;
}

void configInfo::printConfig(){
  std::cout<<"Following Setup was read in:"<<endl;
  for(auto & info : observablePlotInfo)
    std::cout<<"Observable: "<<info.name<<" Histogram Binning: "<<info.bins<< " Min/Max "<<info.min<<"/"<<info.max<<endl;
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
  else if(option==2){
    for(auto & sample : samples)
      if(sample.sampleType==0)names.push_back(sample.nick);
  }
  else if(option==3)
    for(auto & sample : samples)
      if(sample.sampleType==1)names.push_back(sample.nick);
  return names;
}

std::vector<std::string> configInfo::get_Variables(){
  vector<string> var;
  for(auto & info : observablePlotInfo)
    var.push_back(info.name);
  string myweight = weight;
  var.push_back(weight);
  for(auto & cut : cuts){
    var.push_back(cut.variable());
  }
  return var;
}
