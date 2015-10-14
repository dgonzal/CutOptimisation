#include "sample.h"

#include "TChain.h"

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>


#include <sstream>

using namespace std;
using namespace boost;

std::vector<std::string> getMatchingFiles(const std::string & dir, const std::string & pattern){
  namespace fs = boost::filesystem;
  boost::regex exp(pattern);
  boost::regex postfix("root");
  std::vector<std::string> result;
  fs::path full_path = fs::complete(fs::path(dir));
  fs::directory_iterator end_iter;
  for(fs::directory_iterator dir_itr(full_path); dir_itr != end_iter; ++dir_itr){
    if(fs::is_regular_file(dir_itr->status())){
      std::string filename = dir_itr->path().filename().string();
      if(boost::regex_search(filename, exp) && boost::regex_search(filename, postfix)){
	result.push_back(dir_itr->path().string());
      }
    }
  }
  return result;
}

sample::sample(const string & fname_pattern, const std::vector<std::string> & observables_, const string & TreeName): observables(observables_){
  string path;
  size_t p = fname_pattern.rfind('/');
  if(p!=string::npos){
    path = fname_pattern.substr(0, p);
  }
  else p=0;
  string pattern = fname_pattern.substr(p+1);
  std::cout << "sample: using path '" << path << "', pattern '" << pattern << "'. " << std::flush;
  std::vector<std::string> files = getMatchingFiles(path, pattern);
  std::cout << "Found " << files.size() << " file" << (files.size() > 1?"s":"") << std::flush;
  TChain * chain = new TChain(TreeName.c_str());
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
  observables_left.insert(observables.begin(), observables.end());
  for(int i=0; i<br_list->GetEntries(); ++i){
    string br_name = br_list->At(i)->GetName();
    if(observables_left.find(br_name) == observables_left.end()) continue;
    vector<string>::iterator p = find(observables.begin(),observables.end(),br_name);
    obsselindex_to_fileindex[distance(observables.begin(),p)] = i;
    //cout<<br_name<<" "<<i<<" distance "<<distance(observables.begin(),p) <<endl;
    observables_left.erase(br_name);
  }
  // weight is stored in observables, don't forget the case with unityweight!
  if((observables_left.size()>0 && !observables[1].empty()) || observables_left.size()>1){
    stringstream ss;
    cout << "sample: did not find observables '"<<endl;;
    for(set<string>::const_iterator it = observables_left.begin(); it!=observables_left.end(); ++it){
      cout << *it << "' "<<endl;;
    }
    cout << "in File " << fname_pattern;
    throw ss.str();
  }

  // read the data from the input file and "reshuffle" it data structure according to the
  // Set Directory to all given Variables and Weight.
  scoped_array<double_t> data(new double_t[obsselindex_to_fileindex.size()]);	
  for(unsigned int i=0; i<obsselindex_to_fileindex.size(); ++i){
    data[i]=98265;
    if(observables[i].empty())continue;
    chain->SetBranchAddress(br_list->At(obsselindex_to_fileindex[i])->GetName(), &data[i]);
  }
  cout<<" starting to read in"<<endl;
  //Fill Variables into memory
  // resize to the needed size not to waste time allocating new memory several times
  // first to the number of needed observables and then to the entries
  observablesValues.resize(observables.size(),vector<double>(n_events));
  int percent_fraction = n_events/10;
  for(size_t ievent=0; ievent < n_events; ++ievent){
    chain->GetEntry(ievent);
    //Where the data gets into my memory ?
    for(map<size_t, size_t>::const_iterator it=obsselindex_to_fileindex.begin(); it!=obsselindex_to_fileindex.end(); ++it){
      //cout<<it->first<<" second "<<it->second<<endl;
      if(!observables[it->first].empty())
	observablesValues[it->first][ievent] = data[it->first];
      else
	observablesValues[it->first][ievent] = 1;
      if(ievent%percent_fraction==0){
	//cout<<it->first<<" "<<it->second<<" "<<data[it->first]<<" "<<br_list->At(it->second)->GetName() <<" event # "<< ievent<<endl;
	cout<<"\r"<<"Percent done: "<<(ievent*100)/(n_events-1)<<flush;
      }
    }
  }
  cout << std::endl; // all done
  delete chain;
}
