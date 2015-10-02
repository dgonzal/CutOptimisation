#pragma once

#include <string>

#include <boost/shared_array.hpp>

#include "include/obs_selection"
#include "include/OptCut.h"

class sample{
 public:
  std::string theta_systname_and_direction; // empty for nominal, otherwise <syst>__plus or <syst>__minus
  bool is_data() const {
    return nickname.find("MC_") != 0;
  }
  bool is_signal() const {
    return nickname.find("MC_") == 0 && nickname.find("_B") != string::npos;
  }
  bool is_ttbar_bkg() const{
    return nickname.find("TTbar")!=string::npos;
  }
 sample(const std::string & fname_pattern, const std::string & nickname_, const std::string & theta_procname_, const obs_selection & obs_sel_): nickname(nickname_), theta_procname(theta_procname_), obs_sel(obs_sel_)
  const opt_event & operator[](size_t i) const{
    return events[i];
  }
  const string & get_theta_procname() const{
    return theta_procname;
  }  
  const string & get_nickname() const{
    return nickname;
  }
 private:
  size_t n_events;
  boost::shared_array<opt_event> events; //make a shared_array instead of vector to make "copies" of sample cheap
  obs_selection obs_sel; // should be the same for all samples within a program run ... (TODO: enforce ?!)
  std::string nickname; // sample nickname used by gc, etc.
  std::string theta_procname; // theta process name
  enum calc_type {add,multiply,divide};


};
