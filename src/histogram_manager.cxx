#include "histogram_manager.h"

#include "TFile.h"

using namespace std;


HistoManager::HistoManager(){}

void HistoManager::write_HistsToFile(const string & filename){
  TFile * file = new TFile(filename.c_str(), "recreate");
  file->cd();
  for(auto & hist : hists)
    hist->Write();
  delete file;
}

void HistoManager::reset_Hists(){
  for(auto & hist : hists)
    hist->Reset();
}

void HistoManager::fill_ithHist(int i, double val, double weight){
  hists[i]->Fill(val,weight);
}
