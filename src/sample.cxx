#include "include/sample.h"

#include <sstream>      // std::stringstream

using namespace std;
using namespace boost;


sample::sample(const string & fname_pattern, const string & nickname_, const string & theta_procname_, const obs_selection & obs_sel_): nickname(nickname_), theta_procname(theta_procname_), obs_sel(obs_sel_){
    string path;
    size_t p = fname_pattern.rfind('/');
    if(p!=string::npos){
      path = fname_pattern.substr(0, p);
    }
    else p=0;
    string pattern = fname_pattern.substr(p+1);
    cout << "sample: using path '" << path << "', pattern '" << pattern << "'. " << flush;
    vector<string> files = getMatchingFiles(path, pattern);
    cout << "Found " << files.size() << " file" << (files.size() > 1?"s":"") << flush;
    TChain * chain = new TChain("CutTree");
    if(files.size()==0) return;
    for(size_t i=0; i<files.size(); ++i){
      chain->AddFile(files[i].c_str());
    }
    n_events = chain->GetEntries();
    cout << " containing " << n_events << " events." << endl;
    TObjArray* br_list = chain->GetListOfBranches();
    // build a map (observable index according to obs_sel) -> (index into struct from input file)
    map<size_t, size_t> obsselindex_to_fileindex;
    set<string> observables_left;
    observables_left.insert(obs_sel.observable_names.begin(), obs_sel.observable_names.end());
    for(int i=0; i<br_list->GetEntries(); ++i){
      string br_name = br_list->At(i)->GetName();
      if(observables_left.find(br_name) == observables_left.end()) continue;
      size_t obsselindex = obs_sel.indexof_obs(br_name);
      obsselindex_to_fileindex[obsselindex] = i;
      observables_left.erase(br_name);
    }
    if(observables_left.size()>0){
      stringstream ss;
      ss << "sample: did not find observables '";
      for(set<string>::const_iterator it = observables_left.begin(); it!=observables_left.end(); ++it){
	ss << *it << "' ";
      }
      ss << "in fnames " << fname_pattern;
      throw ss.str();
    }
    // read the data from the input file and "reshuffle" it into opt_event structure according to the
    // indices from obs_selection.
    // Set Directory to all given Variables and Weight.
    scoped_array<double_t> data(new double_t[obsselindex_to_fileindex.size()]);	
    float sample_weight = 1.0; //get_sample_weight(sinfo_path, nickname, lumi);
    events.reset(new opt_event[n_events]);
    float sum_of_weights = 0;
    for(unsigned int i=0; i<obsselindex_to_fileindex.size(); ++i){
      data[i]=98265;
      chain->SetBranchAddress(br_list->At(obsselindex_to_fileindex[i])->GetName(), &data[i]);
    }
    //Fill Variables and weight into memory
    opt_event::e_type typ = opt_event::other_bkg;
    if(is_signal()) typ = opt_event::signal;
    for(size_t ievent=0; ievent < n_events; ++ievent){
      opt_event evt(obsselindex_to_fileindex.size()+obs_sel.added_names.size());
      chain->GetEntry(ievent);
      for(map<size_t, size_t>::const_iterator it=obsselindex_to_fileindex.begin(); it!=obsselindex_to_fileindex.end(); ++it){
	//cout<<it->first<<" "<<it->second<<" "<<data[it->second]<<" "<<br_list->At(it->second)->GetName() <<endl;
	evt.data[it->first] = data[it->second];
      }
      for(unsigned int m=0; m<obs_sel.added_names.size();++m){
	evt.data[obsselindex_to_fileindex.size()+m] = case_operation(obs_sel.operation[m],evt.data[obs_sel.indexof(obs_sel.operands[m].first)],evt.data[obs_sel.indexof(obs_sel.operands[m].second)]);
	//cout<< obs_sel.added_names[m]<<" "<<obs_sel.operands[m].first<<" "<<obs_sel.operands[m].second<<" "<<evt.data[obsselindex_to_fileindex.size()+m]<<endl;
      }
      evt.weight = evt.data[obs_sel.indexof("weight")];//w * sample_weight;
      sum_of_weights += evt.weight;
      events[ievent] = evt;
      //print_event(evt,obs_sel );
    }
    delete chain;
  }
