#pragma once

#include <map>
#include <string>
#include <vector>

#include "include/obs_selection.h"
#include "include/sample.h"

#include "TH1D.h"



// all info to create a histogram from a sample for theta.
struct theta_histo_spec{
  int nbins;
  double xmin, xmax;
  boost::shared_ptr<OptCut> cut;
  std::string theta_obsname;
  std::string observable_name; // observable name as in opttuple
  theta_histo_spec(int nbins_, double xmin_, double xmax_, const shared_ptr<OptCut> & cut_, const string & theta_obsname_, const string & obsname_):
  nbins(nbins_), xmin(xmin_), xmax(xmax_), cut(cut_), theta_obsname(theta_obsname_), observable_name(obsname_){}
};

// write the histograms for all histogram specifications into the (newly created) root file fname (which will be overwritten, if it exists).
// As name, <obsname>__<sample.theta_name> is used for all samples. Samples with the same theta_name are added, resulting in only one histogram.
void create_theta_histfile(const std::string & fname, const obs_selection & obs_sel, const std::vector<sample> & samples, const std::vector<theta_histo_spec> & specs){
    TFile * file = new TFile(fname.c_str(), "recreate");
    file->cd();
    std::map<string, TH1D*> theta_histos;
    std::vector<size_t> specs_indices;
    for(size_t j=0; j<specs.size(); ++j){
        specs_indices.push_back(obs_sel.indexof(specs[j].observable_name));
    }
    for(size_t j=0; j<specs.size(); ++j){
        for(size_t isample=0; isample<samples.size(); ++isample){
            string full_histogram_name = specs[j].theta_obsname + "__";
            full_histogram_name += samples[isample].get_theta_procname();
            if(samples[isample].theta_systname_and_direction != ""){
                full_histogram_name += "__" + samples[isample].theta_systname_and_direction;
            }
            TH1D *& histo = theta_histos[full_histogram_name];
            if(histo==0){
                histo = new TH1D(full_histogram_name.c_str(), full_histogram_name.c_str(), specs[j].nbins, specs[j].xmin, specs[j].xmax);
            }
            //int entries = 0;
            for(size_t ievent=0; ievent < samples[isample].size(); ++ievent){
                const opt_event & evt = samples[isample][ievent];
                double sel_weight = (*specs[j].cut)(evt);
                if(sel_weight==0) continue;
                float obs_value = evt.data[specs_indices[j]];
                histo->Fill(obs_value, evt.weight * sel_weight);
                //entries++;
            }
            cout << "histogram " << full_histogram_name << ": " << histo->Integral() << endl;
        }
    }
    file->cd();
    file->Write();
    delete file;
}

// return the expected limits by running theta, as map from * to limit in pb. 
// needs to be updated!!
std::map<std::string, double> expected_limits(const std::string & theta_rootfile, const std::vector<std::string> & signal_processes){
  std::map<string, double> result;
  ofstream ta("analysis.py");
  ta << "execfile(\"cutopt-include.py\")" << endl;
  ta << "run_cutopt(\"" << theta_rootfile << "\")" << endl;
  ta.close();
    
  //run theta-auto:
  char * const args[] = {"analysis.py", 0};
  ProgramWrapper theta_auto((theta_dir + "utils2/theta-auto.py").c_str(), args);
  theta_auto.read_until_line("expected limit");
  for(size_t i=0; i<100; ++i){
    std::string nextline = theta_auto("");
    if(nextline=="end") break;
    size_t p = nextline.find(" ");
    if(p==std::string::npos){
      throw "unexpected line from theta-auto: " + nextline;
    }
    std::string process = nextline.substr(0, p);
    double limit = lexical_cast<double>(nextline.substr(p+1));
    result[process] = limit;
  }
  return result;
}
