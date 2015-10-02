#pragma once 

#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>

#include <vector>
#include <string>

#include "include/obs_selection.h"

struct opt_event{
  float weight;
  enum e_type{
    // we need to distinguish the event type for some reweighting techniques ...
    signal, ttbar_bkg, other_bkg
  };
  e_type event_type;
  boost::shared_array<float> data; // make it shared, not scoped so copying of opt_event is possible ...
  opt_event(size_t n = 0): data(new float[n]){}
};

// cuts on the opt_event level; apply reweight instead of cut.
// 
class OptCut{
public:
  virtual double operator()(const opt_event & evt) =0;
  virtual double passed_events() =0;
};

class PassAll: public OptCut{
public:
  virtual double passed_events(){return 0;};
    double operator()(const opt_event & evt) {
        return true;
    }
};
class AndCut: public OptCut{
public:
  virtual double passed_events(){return 0;};
  void add(const boost::shared_ptr<OptCut> & cut){
    cuts.push_back(cut);
  }
   
    virtual double operator()(const opt_event & evt) {
      double result = 1.0;
      for(size_t i=0; i<cuts.size(); ++i){
	result *= (*cuts[i])(evt);
      }
      return result;
    }
  private:
    std::vector<boost::shared_ptr<OptCut> > cuts;
};
class MultiCut: public OptCut{
public:
  MultiCut(const multi_cutinfo & info_, const obs_selection & obs_sel_): info(info_), obs_sel(obs_sel_) {  
    pass_number = 0.0;
  }
  double passed_events(){return pass_number;}
  double operator()(const opt_event & evt) {
    bool result = false;
    size_t index;
    for(unsigned int i=0; i<info.mode.size();++i){
      index = obs_sel.indexof(info.obsname[i]);
      switch(info.mode[i]){
      case multi_cutinfo::cut_gtr: 
	result = evt.data[index] > *(info.cut_threshold[i]);
	if(!result) return result?1.0:0.0;
	break;
      case multi_cutinfo::cut_lt: 
	result = evt.data[index] < *(info.cut_threshold[i]);
	if(!result) return result?1.0:0.0;
	break;
      default:
	throw string("logic error in MultiCut::operator()");
      }
    }
    pass_number+= evt.weight;
    return result?1.0:0.0;
  }

private:
  double pass_number;
  multi_cutinfo info;
  obs_selection  obs_sel;
};


// It is a bit complicated to use because of the logic one has to consider! 
class CombiCut: public OptCut{
public:
  enum cut_type{OrCut,AndCut};
  virtual double passed_events(){return 0;};
  void add(const boost::shared_ptr<OptCut> & cut,const cut_type type){
    cuts.push_back(cut);
    my_cut.push_back(type);
  }

  void delete_last(int num){
    for(int i=0; i<num; ++i){
      cuts.pop_back();
      my_cut.pop_back();
    }
  }
  std::vector<boost::shared_ptr<OptCut> > cutinfo(){return cuts;} 
  virtual double operator()(const opt_event & evt) {
    double result = 0.0;
    for(size_t i=0; i<cuts.size(); ++i){
      if(my_cut[i]==OrCut)
	result = max(result, (*cuts[i])(evt));
      else if(my_cut[i]==AndCut){
	if(i==0) result = 1.0;
	result *= (*cuts[i])(evt);
      }
      else
	cout <<"error: no cut type in Combi-Cut declared"<<endl;
    }
    return result;
    }
private:
  std::vector<boost::shared_ptr<OptCut> > cuts;
  std::vector<cut_type> my_cut;
};

class OrCut: public OptCut{
public:
  virtual double passed_events(){return 0;};
  void add(const boost::shared_ptr<OptCut> & cut){
    cuts.push_back(cut);
  }
  
  std::vector<boost::shared_ptr<OptCut> > cutinfo(){return cuts;}
    
  virtual double operator()(const opt_event & evt) {
    double result = 0.0;
    for(size_t i=0; i<cuts.size(); ++i){
      result = max(result, (*cuts[i])(evt));
    }
    return result;
    }
private:
  std::vector<boost::shared_ptr<OptCut> > cuts;
};

class CutInfoCut: public OptCut{
public:
  virtual double passed_events(){return 0;};
    CutInfoCut(const cutinfo & info_, const obs_selection & obs_sel): info(info_) {
       index = obs_sel.indexof(info.obsname);
    }  
    virtual double operator()(const opt_event & evt) {
        bool result = false;
        switch(info.mode){
	case cutinfo::cut_gtr: result = evt.data[index] > *(info.cut_threshold); break;
	case cutinfo::cut_lt: result = evt.data[index] < *(info.cut_threshold); break;
	default:
	  throw string("logic error in CutInfoCut::operator()");
        }
        return result?1.0:0.0;
    }
private:
    cutinfo info;
    size_t index;
};

struct cutinfo{
  enum e_mode{ cut_gtr, cut_lt, not_used };
  string obsname;
  e_mode mode;
  pair<double, double> range;
  double stepsize;
  std::vector<e_mode> multi_mode;
  std::vector<pair<double, double> > multi_range;
  std::vector<double> multi_stepsize;
  boost::shared_ptr<double> cut_threshold;
cutinfo(const std::string & obsname_, const e_mode & mode_, const std::pair<double, double> & range_, double stepsize_, const boost::shared_ptr<double> & threshold_):
  obsname(obsname_), mode(mode_), range(range_), stepsize(stepsize_), cut_threshold(threshold_){
  }
};

struct multi_cutinfo{
  enum e_mode{ cut_gtr, cut_lt };
  std::vector<string> obsname;
  std::vector<e_mode> mode;
  std::vector<pair<double, double> > range;
  std::vector<double> stepsize;
  std::vector<boost::shared_ptr<double> > cut_threshold;
multi_cutinfo(const std::vector<string> & obsname_, const std::vector<e_mode> & mode_, const std::vector<std::pair<double, double> > & range_, std::vector<double> stepsize_, const std::vector<shared_ptr<double> > & threshold_):
    obsname(obsname_), mode(mode_), range(range_), stepsize(stepsize_), cut_threshold(threshold_){
  }
};
