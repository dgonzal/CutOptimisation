#pragma once 

#include <vector>
#include <string>

#include <boost/regex.hpp>




std::vector<std::string> getMatchingFiles(const std::string & dir, const std::string & pattern){
  namespace fs = boost::filesystem3;
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


double get_sample_weight(const std::string & path_to_sampleinfo_cfg, const std::string & nickname, double lumi_pb){ 
  return 1.0; //should be changed!
}


void print_event(const opt_event & evt, const obs_selection & obssel){
    for(size_t i=0; i<obssel.observable_names.size(); ++i){
        cout << obssel.observable_names[i] << " = " << evt.data[i] << endl;
    }
    for(size_t i=0; i<obssel.added_names.size(); ++i){
        cout << obssel.added_names[i] << " = " << evt.data[obssel.observable_names.size()+i] << endl;
    }
};

// gauss distribution around 0 with given width ...
double gauss(double width,double center=0){
    static TRandom3 rnd;
    return rnd.Gaus(center, width);
}

void randomize_cuts(vector<cutinfo> & cuts){
    for(size_t i=0; i < cuts.size(); ++i){
        double new_value = *(cuts[i].cut_threshold) + gauss(cuts[i].stepsize);
        if(new_value < cuts[i].range.first) new_value = cuts[i].range.first;
        if(new_value > cuts[i].range.second) new_value = cuts[i].range.second;
        *(cuts[i].cut_threshold) = new_value;
    }
}

vector<double> vary_bounds(vector<double> up_bounds, vector<double> down_bounds, int which_bound, double width, double down= std::numeric_limits<double>::min(), double up= std::numeric_limits<double>::max()){
  vector<double> result_vector; 	
  double result;
  if(which_bound==0){
    for(unsigned int i=0; i<up_bounds.size();++i){
      if(up_bounds[i]>=down && up_bounds[i]<=up){
	int count = 0;
	do{
	  result = up_bounds[i]+gauss(width);
	  count+=1;
	}while(count<100 &&( result<down || result>up || result<down_bounds[i]));
	if(result<down || result>up || result< down_bounds[i]) result=up_bounds[i];
      }
      else
	result = up_bounds[i];
      result_vector.push_back(result);
    }
  }
  else if(which_bound==1){
    for(unsigned int i=0; i<down_bounds.size();++i){
      if(down_bounds[i]>=down && down_bounds[i]<=up){
	int count = 0;
	do{
	  result = down_bounds[i]+gauss(width);
	  count+=1;
	}while(count<100 &&( result<down || result>up || result> up_bounds[i]));
	if(result<down || result>up || result> up_bounds[i]) result=down_bounds[i];
      }
      else
	result = down_bounds[i];
      result_vector.push_back(result);
    }
    
  }
  return result_vector;
}
