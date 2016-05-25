#include "histogram_manager.h"

#include "TFile.h"

using namespace std;


HistoManager::HistoManager(){}

void HistoManager::write_HistsToFile(const string & filename){
  TFile * file = new TFile(filename.c_str(), "recreate");
  file->cd();
  for(auto & hist : hists)
    hist->Write();
  for(auto & vec : sys)
    for(auto & hist :vec)
      hist->Write();
  delete file;
}

void HistoManager::reset_Hists(){
  for(auto & hist : hists)
    hist->Reset();
  for(auto & vec : sys)
    for(auto & hist :vec)
      hist->Reset();
}

void HistoManager::fill_ithHist(int i, double val, double weight, vector<double> sys_weights){
  hists[i]->Fill(val,weight);
  for(unsigned int p =0; p < sys_weights.size(); p++)
    sys[i][p]->Fill(val,sys_weights[p]);
}

void HistoManager::create_Hists(const obsInfo & PlotInfo ,const std::string & sampleName,const std::vector<std::string> & SysNames){
  hists.push_back(new TH1D((PlotInfo.name+"__"+sampleName).c_str(),(PlotInfo.name+"__"+sampleName).c_str(),PlotInfo.bins,PlotInfo.min,PlotInfo.max));
  vector<TH1D*> sys_hists;
  for(auto systematic : SysNames){
    if(boost::algorithm::ends_with(systematic,"__plus") || boost::algorithm::ends_with(systematic,"__minus"))
      sys_hists.push_back(new TH1D((PlotInfo.name+"__"+sampleName+"__"+systematic).c_str(),(PlotInfo.name+"__"+sampleName+"__"+systematic).c_str(),PlotInfo.bins,PlotInfo.min,PlotInfo.max));
    else if(boost::algorithm::contains(systematic,"Up") || boost::algorithm::contains(systematic,"up") ||  boost::algorithm::contains(systematic,"plus"))
      sys_hists.push_back(new TH1D((PlotInfo.name+"__"+sampleName+"__"+systematic+"__plus").c_str(),(PlotInfo.name+"__"+sampleName+"__"+systematic).c_str(),PlotInfo.bins,PlotInfo.min,PlotInfo.max));
    else if(boost::algorithm::contains(systematic,"down") || boost::algorithm::contains(systematic,"Down") ||  boost::algorithm::contains(systematic,"minus"))
      sys_hists.push_back(new TH1D((PlotInfo.name+"__"+sampleName+"__"+systematic+"__minus").c_str(),(PlotInfo.name+"__"+sampleName+"__"+systematic).c_str(),PlotInfo.bins,PlotInfo.min,PlotInfo.max));
  }
  sys.push_back(sys_hists);
}
