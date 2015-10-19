
#include "config_parser.h"
#include "sample.h"
#include "histogram_manager.h"
#include "CutManager.h"
#include "Cut.h"
#include "PythonWrapper.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>

#include <iostream>
#include <fstream>

using namespace std;

bool run_theta(configInfo & config, vector<sample> & samples);


int main(int argc, char* argv[]){
  /**
   ** TO DO: write help command & securit check
   **
   */
  bool debug = false;
  if(argc < 2 || argc > 2){
    cout<<"One INI File should be given as Argument"<<endl;
    return 1;
  }
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini(argv[1], pt);
  configInfo currentConfig(pt);
  if(debug)currentConfig.printConfig();
  if(!boost::filesystem::exists(currentConfig.get_Workdir()))
    boost::filesystem::create_directories(currentConfig.get_Workdir());
  else
    cout<<"Workdir already exists, I will ignore this and write stuff in there for now"<<endl;
  chdir(currentConfig.get_Workdir().c_str());
  //load sample trees  
  vector<sample> sample_vec;
  vector<string> observables = currentConfig.get_Variables();
  for(auto & data :  currentConfig.get_Samples()){
    sample mysample(currentConfig.get_Path()+data.dir,observables);
    sample_vec.push_back(mysample);
  }
  run_theta(currentConfig,sample_vec);

  return 0;
}

struct resultInfo{
  string cuts;
  vector<double> limits;
  double compareValue;
};

bool run_theta(configInfo & config, vector<sample> & samples){
  ofstream out("result.txt");
  out<<"CutName ";
  for(auto & sampleName : config.get_SampleInformation(2))
    out<<sampleName<<" ";
  out<<endl;
  vector<resultInfo> best_results;
  CutManager ScanCuts;
  for(auto & cut : config.get_Cuts())
    ScanCuts.add_Cut(cut.nick, cut.cutname(),cut.cutInput());
  unsigned int numberOfCuts = ScanCuts.number_Cuts();
  HistoManager ThetaHists;
  for(auto & sampleName : config.get_SampleInformation(1)){
    ThetaHists.hists.push_back(new TH1D((config.get_Observable()+"__"+sampleName).c_str(),sampleName.c_str(),config.get_PlotInfo().bins,config.get_PlotInfo().min,config.get_PlotInfo().max));
  }
  cout<<" number of scanned cuts "<<numberOfCuts<<" total number of points "<<ScanCuts.get_NumberScanPoints()<<endl;
  boost::progress_display show_progress(ScanCuts.get_NumberScanPoints());
  for(int p =0; p<ScanCuts.get_NumberScanPoints(); p++){
    ScanCuts.permute_Cuts(p);
    int i =0;
    //cout<<int(p%10)<<" "<<int((p/10)%20)<<endl;
    for(auto & sample : samples){
      for(int m = 0; m< sample.get_NumberEvents(); m++){
	bool passCuts = true;
	//int frac_samp = int(double(m)/double(sample.get_NumberEvents())*100);
	//if(frac_samp%10)cout<<"\r"<<"Scan Progress [%] "<<int(double(p)/double(ScanCuts.get_NumberScanPoints())*100) <<" "<<ScanCuts.get_OuputName() <<" "<<config.get_Samples()[i].nick<<" [%]: "<<frac_samp<<flush;// " pass "<< pass << " fail "<<fail<<flush;
	for(unsigned int it =0; it <numberOfCuts; it++){
	  passCuts = ScanCuts.check_Cut(it, sample.get_Value(m,it+2));
	  if(!passCuts) break;
	  if(it+1==numberOfCuts)ThetaHists.hists[i]->Fill(sample.get_Value(m,0),sample.get_Value(m,1));
	}
      }
      i++; 
    }
    resultInfo myInfo;
    myInfo.cuts = ScanCuts.get_OuputName();
    myInfo.compareValue =1;
    ThetaHists.write_HistsToFile(ScanCuts.get_OuputName()+".root");
    map<string, double> limits = expected_limits(config.get_ThetaDir(),ScanCuts.get_OuputName()+".root");
    out<<ScanCuts.get_OuputName()<<" ";
    for(map<string, double>::const_iterator it=limits.begin(); it!=limits.end(); ++it){
      out << it->second<<" ";
      myInfo.limits.push_back(it->second);
      myInfo.compareValue *=it->second;
    }
    out<<endl;
    if(p<5) best_results.push_back(myInfo);
    else{
      int position = -1;
      double distance = 0;
      for(unsigned int ip = 0; ip < 5;ip++){
	if(myInfo.compareValue<best_results[ip].compareValue && distance < fabs(myInfo.compareValue-best_results[ip].compareValue)){
	  position = ip;
	  distance = fabs(myInfo.compareValue-best_results[ip].compareValue);
	}
      }
      if(position != -1)best_results[position]= myInfo;
    }		   
    ThetaHists.reset_Hists();
    ++show_progress;
  }
  out<<"================================================================="<<endl;
  for(auto & result: best_results){
    out<<result.cuts<<" ";
    for(auto & number : result.limits)
      out<<number<<" ";
    out<<endl;
  }
  
  out.close();
  return true;
}

/*
** Ideas:
** -> Parallel reading of samples, don't forget to block mem if writting to it
** -> Split work load into making histograms & rebinning them and running theta
**
** -> Put stuff into a class that manage the histo gram production
** -> most of the stuff in in run_theta is rather generic and could go in main as part of a class
**
** -> In principale one could also disintangle theta computation and hist production 
**    if something goes wrong would could restart it easily at any given point
**
*/
