#include "config_parser.h"
#include "sample.h"
#include "histogram_manager.h"
#include "CutManager.h"
#include "Cut.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/filesystem.hpp>

#include <iostream>

using namespace std;

bool run_theta(configInfo & config, vector<sample> & samples);


int main(int argc, char* argv[]){
  /**
   ** TO DO: write help command & securit check
   **
  if(argv.size()<2 || argv.size()>2){
    cout<<"One INI File should be given as Argument"<<endl;
    return 1;
  }
  */
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini(argv[1], pt);
  configInfo currentConfig(pt);
  currentConfig.printConfig();
  //load sample trees  
  vector<sample> sample_vec;
  vector<string> observables = currentConfig.get_Variables();
  for(auto & data :  currentConfig.get_Samples()){
    sample mysample(currentConfig.get_Path()+data.dir,observables);
    sample_vec.push_back(mysample);
  }
  if(!boost::filesystem::exists(currentConfig.get_Workdir()))
    boost::filesystem::create_directories(currentConfig.get_Workdir());
  else
    cout<<"Workdir already exists, I will ignore this and write stuff in there for now"<<endl;

  run_theta(currentConfig,sample_vec);

  return 0;
}

bool run_theta(configInfo & config, vector<sample> & samples){
  CutManager ScanCuts;
  for(auto & cut : config.get_Cuts())
    ScanCuts.add_Cut(cut.nick, cut.cutname(),cut.cutInput());
  unsigned int numberOfCuts = ScanCuts.number_Cuts();
  HistoManager ThetaHists;
  for(auto & sampleName : config.get_SampleInformation(1)){
    ThetaHists.hists.push_back(new TH1D(sampleName.c_str(),sampleName.c_str(),config.get_PlotInfo().bins,config.get_PlotInfo().min,config.get_PlotInfo().max));
  }
  
  cout<<" number of scanned cuts "<<numberOfCuts<<" total number of points "<<ScanCuts.get_NumberScanPoints()<<endl;
  for(int p =0; p<ScanCuts.get_NumberScanPoints(); p++){
    int i =0;
    //cout<<int(p%10)<<" "<<int((p/10)%20)<<endl;
    for(auto & sample : samples){
      int fail=0 ,pass=0;
      for(int m = 0; m< sample.get_NumberEvents(); m++){
	bool passCuts = true;
	cout<<"\r"<<" "<<config.get_Samples()[i].nick<<" [%]: "<<(double(m)/double(sample.get_NumberEvents())*100)<< " pass "<< pass << " fail "<<fail<<flush;
	unsigned int it = 0;
	do{
	  passCuts = ScanCuts.check_Cut(it, sample.get_Value(m,it+2));
	  it++;
	}while(it<numberOfCuts && passCuts==true);
	if(passCuts) 
	  pass++;
	else
	  fail++;
	ThetaHists.hists[i]->Fill(sample.get_Value(m,0),sample.get_Value(m,1));
      }
      i++; 
    }
    ThetaHists.write_HistsToFile(config.get_Workdir()+"/"+ScanCuts.get_OuputName()+".root");
    ThetaHists.reset_Hists();
    ScanCuts.permute_Cuts(p);
  }
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
