// This is the cut optimizer for the Z'->ttbar analysis. It uses ntuples produced by highmass_tuple
// as input.
// Thanks to Jochen Ott for the code.
#include "ProgramWrapper.hpp"
#include "PythonFunction.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TRandom3.h"

#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/scoped_array.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <string>
#include <stdexcept>
#include <limits>


using namespace std;
using namespace boost;


//const string sinfo_path = "/";
const string theta_dir = "/scratch/hh/dust/naf/cms/user/gonzalez/CMSSW/theta/";

// rules for samples:
// * for processing, each ("grid") sample is identified by the sample nickname. Each sample on this level
//  corresponds to one sample instance.
// * for the statistical evaluation, several samples might be added together (such as t and tbar in single top, or even different single top processes);
//   this is done via the theta_name.
//
// all samples used within one run of the program must be constructed with the same obs_selection.
struct obs_selection{
  enum calc_type{multiply,divide,add, substract};
  vector<string> observable_names;
  vector<string> added_names;//variables derivated from obs_names
  vector<operands_type> operands;
  vector<calc_type> operation;
  size_t indexof_obs(const string & observable_name) const {
    size_t index = find(observable_names.begin(), observable_names.end(), observable_name) - observable_names.begin();
    if(index == observable_names.size()){
      throw string("sample_info::indexof: did not find observable '" + observable_name + "'");
    }
    return index;
  }
  size_t indexof(const string & observable_name) const {
    size_t index = find(observable_names.begin(), observable_names.end(), observable_name) - observable_names.begin();
    size_t add_index = find(added_names.begin(), added_names.end(), observable_name) - added_names.begin();
    if(index == observable_names.size() && add_index == added_names.size()){
      throw string("sample_info::indexof: did not find observable '" + observable_name + "'");
    }
    if(add_index == added_names.size()) add_index =0;

    return index+add_index;
  }
};

double case_operation(obs_selection::calc_type operation, double operand1, double operand2){
  double result =0.;
  switch(operation){
  case obs_selection::multiply: 
    result = operand1*operand2;
    break;  
  case obs_selection::divide: 
    result = operand1/operand2;
    //cout<<result<<" = "<<operand1<<" / "<<operand2<<endl;
    break;
  case obs_selection::add: 
    result = operand1+operand2;
    break;
  case obs_selection::substract: 
    result = operand1-operand2;
    break;
  default:
    throw string("logic error in case_operation()");
  }
  return result;
}

void highmass_multi_srs(const string & theta_fname, const vector<sample> & samples, const obs_selection & obs_sel, double htlepcut,
                       const vector<boost::shared_ptr<OptCut> > & sr_regions, const boost::shared_ptr<OptCut> & htlepsb){
    vector<theta_histo_spec> specs;
    for(size_t i=0; i<sr_regions.size(); ++i){
      stringstream ss;
      ss << "mu_mtt" << i;
      specs.push_back(theta_histo_spec(80, 0, 4000, sr_regions[i], ss.str(), "mttbar"));
    }
    specs.push_back(theta_histo_spec(20, 0, htlepcut, htlepsb, "mu_htlep", "htlep"));
    create_theta_histfile(theta_fname, obs_sel, samples, specs);
}

void selection_multi_srs(const string & theta_fname, const vector<sample> & samples, const obs_selection & obs_sel, 
                       const vector<boost::shared_ptr<OptCut> > & sr_regions){
    vector<theta_histo_spec> specs;
    for(size_t i=0; i<sr_regions.size(); ++i){
      stringstream ss;
      ss << "mu_mtt" << i;
      specs.push_back(theta_histo_spec(100, 0, 4000, sr_regions[i], ss.str(), "mtt_clean"));
    }
    create_theta_histfile(theta_fname, obs_sel, samples, specs);
    PythonFunction rebinning;
    rebinning.rebinning("dsss",theta_fname.c_str());
}

int run(){
    obs_selection obs_sel;
    obs_sel.observable_names.push_back("weight");
    obs_sel.observable_names.push_back("pT_mu");
    obs_sel.observable_names.push_back("HT");
    obs_sel.observable_names.push_back("2D");
    obs_sel.observable_names.push_back("mtt");
    obs_sel.observable_names.push_back("met_pt");
    obs_sel.observable_names.push_back("Jet_pt_max");
    obs_sel.observable_names.push_back("btag");

    obs_sel.observable_names.push_back("Iso02");
    obs_sel.observable_names.push_back("Iso018");
    obs_sel.observable_names.push_back("Iso016");
    obs_sel.observable_names.push_back("Iso014");
    obs_sel.observable_names.push_back("Iso012");
    obs_sel.observable_names.push_back("Iso01");
    obs_sel.observable_names.push_back("Iso008");
    obs_sel.observable_names.push_back("Iso006");
    obs_sel.observable_names.push_back("Iso004");
    obs_sel.observable_names.push_back("Iso002");
    obs_sel.observable_names.push_back("Iso04");
    obs_sel.observable_names.push_back("Iso05");

    obs_sel.observable_names.push_back("toptag");
    obs_sel.observable_names.push_back("chi2");
    obs_sel.observable_names.push_back("mtt_clean");
    obs_sel.observable_names.push_back("heptoptag");
   
    obs_sel.added_names.push_back("HTl");
    obs_sel.operands.push_back(make_pair("pT_mu","met_pt"));
    obs_sel.operation.push_back(obs_selection::add);	

    obs_sel.added_names.push_back("relIso02");
    obs_sel.operands.push_back(make_pair("Iso02","pT_mu"));
    obs_sel.operation.push_back(obs_selection::divide);
    
    obs_sel.added_names.push_back("relIso018");
    obs_sel.operands.push_back(make_pair("Iso018","pT_mu"));
    obs_sel.operation.push_back(obs_selection::divide);
    obs_sel.added_names.push_back("relIso016");
    obs_sel.operands.push_back(make_pair("Iso016","pT_mu"));
    obs_sel.operation.push_back(obs_selection::divide);
    obs_sel.added_names.push_back("relIso014");
    obs_sel.operands.push_back(make_pair("Iso014","pT_mu"));
    obs_sel.operation.push_back(obs_selection::divide);
    obs_sel.added_names.push_back("relIso012");
    obs_sel.operands.push_back(make_pair("Iso012","pT_mu"));
    obs_sel.operation.push_back(obs_selection::divide);
   
    obs_sel.added_names.push_back("relIso01");
    obs_sel.operands.push_back(make_pair("Iso01","pT_mu"));
    obs_sel.operation.push_back(obs_selection::divide);
    obs_sel.added_names.push_back("relIso008");
    obs_sel.operands.push_back(make_pair("Iso008","pT_mu"));
    obs_sel.operation.push_back(obs_selection::divide);
    obs_sel.added_names.push_back("relIso006");
    obs_sel.operands.push_back(make_pair("Iso006","pT_mu"));
    obs_sel.operation.push_back(obs_selection::divide);
    obs_sel.added_names.push_back("relIso004");
    obs_sel.operands.push_back(make_pair("Iso004","pT_mu"));
    obs_sel.operation.push_back(obs_selection::divide);
    obs_sel.added_names.push_back("relIso002");
    obs_sel.operands.push_back(make_pair("Iso002","pT_mu"));
    obs_sel.operation.push_back(obs_selection::divide);
   
    obs_sel.added_names.push_back("relIso04");
    obs_sel.operands.push_back(make_pair("Iso04","pT_mu"));
    obs_sel.operation.push_back(obs_selection::divide);

    obs_sel.added_names.push_back("relIso05");
    obs_sel.operands.push_back(make_pair("Iso05","pT_mu"));
    obs_sel.operation.push_back(obs_selection::divide);

    string path = "/scratch/hh/dust/naf/cms/user/gonzalez/Analysis53X_v3/heptest/";
    vector<sample> samples;
    samples.push_back(sample(path + "QCDCycle.MC.TTbar", "MC_TTbar", "ttbar", obs_sel));
    samples.push_back(sample(path + "QCDCycle.MC.QCDMu", "MC_QCD", "qcd", obs_sel));
    samples.push_back(sample(path + "QCDCycle.MC.W", "MC_WJets", "wjets", obs_sel));
    samples.push_back(sample(path + "QCDCycle.MC.DY", "MC_Zjets", "zjets", obs_sel));

    samples.push_back(sample(path + "QCDCycle.MC.ZP1000w10.root", "MC_ZP1000w10", "zp1000", obs_sel));
    samples.push_back(sample(path + "QCDCycle.MC.ZP1250w12p5.root", "MC_ZP1250w15", "zp1250", obs_sel));
    samples.push_back(sample(path + "QCDCycle.MC.ZP1500w15.root", "MC_ZP1500w15", "zp1500", obs_sel));
    samples.push_back(sample(path + "QCDCycle.MC.ZP2000w20.root", "MC_ZP2000w20", "zp2000", obs_sel));    
    samples.push_back(sample(path + "QCDCycle.MC.ZP3000w30.root", "MC_ZP3000w30", "zp3000", obs_sel));
    samples.push_back(sample(path + "QCDCycle.MC.ZP4000w40.root", "MC_ZP4000w40", "zp4000", obs_sel));
    //define cuts:
    vector<cutinfo> cuts;    
    //2D-Cut
    shared_ptr<double> TwoD_cut_bool(new double(0.5));
    cutinfo TwoD_cutinfo("2D", cutinfo::cut_gtr, make_pair(0.0,1.), .1, TwoD_cut_bool);
    cuts.push_back(TwoD_cutinfo); 
    
    //cone-cuts
    shared_ptr<double> iso05_double(new double(5));
    cutinfo Iso05cutinfo("Iso05", cutinfo::cut_lt, make_pair(0,50), .5, iso05_double);
    cuts.push_back(Iso05cutinfo); 

    shared_ptr<double> reliso05_double(new double(0.2));
    cutinfo relIso05cutinfo("relIso05", cutinfo::cut_lt, make_pair(0,10), .01, reliso05_double);
    cuts.push_back(relIso05cutinfo); 
    
    shared_ptr<double> reliso02_double(new double(0.2));
    cutinfo relIso02cutinfo("relIso02", cutinfo::cut_lt, make_pair(0,10), .01, reliso02_double);
    cuts.push_back(relIso02cutinfo); 

    shared_ptr<double> reliso008_double(new double(0.2));
    cutinfo relIso008cutinfo("relIso008", cutinfo::cut_lt, make_pair(0,10), .01, reliso008_double);
    cuts.push_back(relIso008cutinfo); 

    shared_ptr<double> reliso04_double(new double(0.2));
    cutinfo relIso04cutinfo("relIso04", cutinfo::cut_lt, make_pair(0,10), .01, reliso04_double);
    cuts.push_back(relIso04cutinfo); 

    //vanilla cuts
    shared_ptr<double> met_cut(new double(50));
    cutinfo metcutinfo("met_pt", cutinfo::cut_gtr, make_pair(0,500), 10, met_cut);
    cuts.push_back(metcutinfo);
    shared_ptr<double> ht_cut(new double(150));
    cutinfo htcutinfo("HTl", cutinfo::cut_gtr, make_pair(0,3500), 5, ht_cut);
    cuts.push_back(htcutinfo);
    shared_ptr<double> Jet_pt_max_cut(new double(150));
    cutinfo Jet_pt_maxcutinfo("Jet_pt_max", cutinfo::cut_gtr, make_pair(0,900), 50,Jet_pt_max_cut);
    cuts.push_back(Jet_pt_maxcutinfo);
    shared_ptr<double> Jet_pt_max210_cut(new double(210));
    cutinfo Jet_pt_210maxcutinfo("Jet_pt_max", cutinfo::cut_gtr, make_pair(0,900), 50,Jet_pt_max210_cut);
    cuts.push_back(Jet_pt_210maxcutinfo);

    shared_ptr<double> Btag_cut(new double(0.5));
    cutinfo Btag_cutinfo("btag", cutinfo::cut_gtr, make_pair(0,10), 0.5,Btag_cut);
    cuts.push_back(Btag_cutinfo);
    cutinfo NoBtag_cutinfo("btag", cutinfo::cut_lt, make_pair(0,10), 0.5,Btag_cut);
    cuts.push_back(NoBtag_cutinfo);

    shared_ptr<double> Toptag_cut(new double(0.5));
    cutinfo Toptag_cutinfo("toptag", cutinfo::cut_gtr, make_pair(0,10), 0.5,Toptag_cut);
    cuts.push_back(Toptag_cutinfo);
    cutinfo NoToptag_cutinfo("toptag", cutinfo::cut_lt, make_pair(0,10), 0.5,Toptag_cut);
    cuts.push_back(NoToptag_cutinfo);

    shared_ptr<double> Chi2_cut(new double(30));
    cutinfo Chi2_cutinfo("chi2", cutinfo::cut_lt, make_pair(0,1000), 0.5,Chi2_cut);
    cuts.push_back(Chi2_cutinfo);

    shared_ptr<double> Chi2_45cut(new double(45));
    cutinfo Chi2_45cutinfo("chi2", cutinfo::cut_lt, make_pair(0,1000), 0.5,Chi2_45cut);
    cuts.push_back(Chi2_45cutinfo);
    shared_ptr<double> Chi2_35cut(new double(35));
    cutinfo Chi2_35cutinfo("chi2", cutinfo::cut_lt, make_pair(0,1000), 0.5,Chi2_35cut);
    cuts.push_back(Chi2_35cutinfo);

    shared_ptr<double> Chi2_cut10(new double(10));
    cutinfo Chi2_10cutinfo("chi2", cutinfo::cut_lt, make_pair(0,1000), 0.5,Chi2_cut10);
    cuts.push_back(Chi2_10cutinfo);

    shared_ptr<PassAll> allevents(new PassAll());
    
    shared_ptr<AndCut> twodcut_btag(new AndCut());
    twodcut_btag->add(shared_ptr<OptCut>(new CutInfoCut(TwoD_cutinfo, obs_sel)));
    twodcut_btag->add(shared_ptr<OptCut>(new CutInfoCut(metcutinfo, obs_sel)));	     
    twodcut_btag->add(shared_ptr<OptCut>(new CutInfoCut(htcutinfo, obs_sel)));	     
    twodcut_btag->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_maxcutinfo, obs_sel)));
    twodcut_btag->add(shared_ptr<OptCut>(new CutInfoCut(Btag_cutinfo, obs_sel)));
    //twodcut_btag->add(shared_ptr<OptCut>(new CutInfoCut(Toptag_cutinfo, obs_sel)));
    twodcut_btag->add(shared_ptr<OptCut>(new CutInfoCut(Chi2_10cutinfo, obs_sel)));

    shared_ptr<AndCut> twodcut_nobtag(new AndCut());
    twodcut_nobtag->add(shared_ptr<OptCut>(new CutInfoCut(TwoD_cutinfo, obs_sel)));
    twodcut_nobtag->add(shared_ptr<OptCut>(new CutInfoCut(metcutinfo, obs_sel)));	     
    twodcut_nobtag->add(shared_ptr<OptCut>(new CutInfoCut(htcutinfo, obs_sel)));	     
    twodcut_nobtag->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_maxcutinfo, obs_sel)));
    twodcut_nobtag->add(shared_ptr<OptCut>(new CutInfoCut(NoBtag_cutinfo, obs_sel)));
    //twodcut_nobtag->add(shared_ptr<OptCut>(new CutInfoCut(NoToptag_cutinfo, obs_sel)));
    twodcut_nobtag->add(shared_ptr<OptCut>(new CutInfoCut(Chi2_10cutinfo, obs_sel)));

    shared_ptr<AndCut> iso05_cut(new AndCut());
    iso05_cut->add(shared_ptr<OptCut>(new CutInfoCut(Iso05cutinfo, obs_sel)));
    iso05_cut->add(shared_ptr<OptCut>(new CutInfoCut(metcutinfo, obs_sel)));	     
    iso05_cut->add(shared_ptr<OptCut>(new CutInfoCut(htcutinfo, obs_sel)));	     
    iso05_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_maxcutinfo, obs_sel)));

    shared_ptr<AndCut> reliso04_cut(new AndCut());
    reliso04_cut->add(shared_ptr<OptCut>(new CutInfoCut(relIso04cutinfo, obs_sel)));
    reliso04_cut->add(shared_ptr<OptCut>(new CutInfoCut(metcutinfo, obs_sel)));	     
    reliso04_cut->add(shared_ptr<OptCut>(new CutInfoCut(htcutinfo, obs_sel)));	     
    reliso04_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_maxcutinfo, obs_sel)));

    shared_ptr<AndCut> reliso04_btag(new AndCut());
    reliso04_btag->add(shared_ptr<OptCut>(new CutInfoCut(relIso04cutinfo, obs_sel)));
    reliso04_btag->add(shared_ptr<OptCut>(new CutInfoCut(metcutinfo, obs_sel)));	     
    reliso04_btag->add(shared_ptr<OptCut>(new CutInfoCut(htcutinfo, obs_sel)));	     
    reliso04_btag->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_maxcutinfo, obs_sel)));
    reliso04_btag->add(shared_ptr<OptCut>(new CutInfoCut(Btag_cutinfo, obs_sel)));
    reliso04_btag->add(shared_ptr<OptCut>(new CutInfoCut(Chi2_cutinfo, obs_sel)));

    shared_ptr<AndCut> reliso04_nobtag(new AndCut());
    reliso04_nobtag->add(shared_ptr<OptCut>(new CutInfoCut(relIso04cutinfo, obs_sel)));
    reliso04_nobtag->add(shared_ptr<OptCut>(new CutInfoCut(metcutinfo, obs_sel)));	     
    reliso04_nobtag->add(shared_ptr<OptCut>(new CutInfoCut(htcutinfo, obs_sel)));	     
    reliso04_nobtag->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_maxcutinfo, obs_sel)));
    reliso04_nobtag->add(shared_ptr<OptCut>(new CutInfoCut(NoBtag_cutinfo, obs_sel)));
    reliso04_nobtag->add(shared_ptr<OptCut>(new CutInfoCut(Chi2_cutinfo, obs_sel)));

    shared_ptr<AndCut> reliso05_cut(new AndCut());
    reliso05_cut->add(shared_ptr<OptCut>(new CutInfoCut(relIso05cutinfo, obs_sel)));
    reliso05_cut->add(shared_ptr<OptCut>(new CutInfoCut(metcutinfo, obs_sel)));	     
    reliso05_cut->add(shared_ptr<OptCut>(new CutInfoCut(htcutinfo, obs_sel)));	     
    reliso05_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_maxcutinfo, obs_sel)));

    shared_ptr<AndCut> reliso02_cut(new AndCut());
    reliso02_cut->add(shared_ptr<OptCut>(new CutInfoCut(relIso04cutinfo, obs_sel)));
    reliso02_cut->add(shared_ptr<OptCut>(new CutInfoCut(relIso02cutinfo, obs_sel)));
    reliso02_cut->add(shared_ptr<OptCut>(new CutInfoCut(metcutinfo, obs_sel)));	     
    reliso02_cut->add(shared_ptr<OptCut>(new CutInfoCut(htcutinfo, obs_sel)));	     
    reliso02_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_maxcutinfo, obs_sel)));

    shared_ptr<CombiCut> twoD_btag_cut(new CombiCut());
    shared_ptr<CombiCut> twoD_nobtag_cut(new CombiCut());
    shared_ptr<CombiCut> twoD_nobtag2_cut(new CombiCut());

    //twoD_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(relIso04cutinfo,obs_sel)),CombiCut::OrCut);
    twoD_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(TwoD_cutinfo, obs_sel)),CombiCut::AndCut);
    twoD_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(metcutinfo, obs_sel)),CombiCut::AndCut);
    twoD_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(htcutinfo, obs_sel)),CombiCut::AndCut);
    twoD_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_maxcutinfo, obs_sel)),CombiCut::AndCut);
    //twoD_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_210maxcutinfo, obs_sel)),CombiCut::AndCut);
    twoD_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Btag_cutinfo, obs_sel)),CombiCut::AndCut);
    twoD_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Toptag_cutinfo, obs_sel)),CombiCut::AndCut);
    twoD_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Chi2_10cutinfo, obs_sel)),CombiCut::AndCut);

    //twoD_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(relIso04cutinfo,obs_sel)),CombiCut::OrCut);
    twoD_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(TwoD_cutinfo, obs_sel)),CombiCut::AndCut);
    twoD_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(metcutinfo, obs_sel)),CombiCut::AndCut);
    twoD_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(htcutinfo, obs_sel)),CombiCut::AndCut);
    twoD_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_maxcutinfo, obs_sel)),CombiCut::AndCut);
    //twoD_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_maxcutinfo, obs_sel)),CombiCut::AndCut);
    twoD_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(NoBtag_cutinfo, obs_sel)),CombiCut::AndCut);
    twoD_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(NoToptag_cutinfo, obs_sel)),CombiCut::AndCut);
    twoD_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Chi2_10cutinfo, obs_sel)),CombiCut::AndCut);

    //twoD_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(relIso04cutinfo,obs_sel)),CombiCut::OrCut);
    twoD_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(TwoD_cutinfo, obs_sel)),CombiCut::AndCut);
    twoD_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(metcutinfo, obs_sel)),CombiCut::AndCut);
    twoD_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(htcutinfo, obs_sel)),CombiCut::AndCut);
    twoD_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_maxcutinfo, obs_sel)),CombiCut::AndCut);
    //twoD_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_210maxcutinfo, obs_sel)),CombiCut::AndCut);
    twoD_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(Btag_cutinfo, obs_sel)),CombiCut::AndCut);
    twoD_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(NoToptag_cutinfo, obs_sel)),CombiCut::AndCut);
    twoD_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(Chi2_10cutinfo, obs_sel)),CombiCut::AndCut);

    shared_ptr<CombiCut> start_btag_cut(new CombiCut());
    shared_ptr<CombiCut> start_nobtag_cut(new CombiCut());
    shared_ptr<CombiCut> start_nobtag2_cut(new CombiCut());
  
    //testing and improvments
    vector<double> pt_bounds_up;
    vector<double> pt_bounds_down;
    vector<double> iso_vec;
    vector<double> relIso_vec;
    
    relIso_vec.push_back(0.2);
   
    relIso_vec.push_back(0.18);
    relIso_vec.push_back(0.16);
    relIso_vec.push_back(0.14);
    relIso_vec.push_back(0.12);
    
    relIso_vec.push_back(0.1);
    relIso_vec.push_back(0.08);
    relIso_vec.push_back(0.06);
    relIso_vec.push_back(0.04);
    relIso_vec.push_back(0.02);
    //relIso_vec.push_back(0.01);
    

    for(unsigned int i=0; i<relIso_vec.size(); ++i){
      //double limit = 47.8514/(relIso_vec[i]-0.0250832)-200.501;
      //double limit = 24.578/(relIso_vec[i]-0.0374231)-97.4084;
      double limit =  29.1356/(relIso_vec[i]- 0.023111)-164.383;
      //double limit =45/(relIso_vec[i]-0.0206224)-210;
      cout<<limit<<endl;
    
      if(limit>0 && limit <2000){
	pt_bounds_up.push_back(limit);
	pt_bounds_down.push_back(limit);
      }
      else if(limit >=2000){
	pt_bounds_up.push_back(1900);
	pt_bounds_down.push_back(1900);

      }

    }
    pt_bounds_up.push_back(2000);
    pt_bounds_down.push_back(2000);


    for(unsigned int i=0; i<obs_sel.added_names.size()-2;++i){

      //if(obs_sel.added_names.size()-1!=i)pt_bounds.push_back(100+(double)i*1000);
      //if(obs_sel.added_names.size()-1==i)pt_bounds.push_back(10000);

      iso_vec.push_back(0.2);

      vector<shared_ptr<double> > threshold;
      vector<multi_cutinfo::e_mode> multi_mode;
      vector<pair<double, double> > range_pairs;
      vector<double> stepsize;
      vector<string> obsnames;
      
      cout<<pt_bounds_up[i]<<" "<<pt_bounds_up[i+1]<<endl;

      pair<double, double> pt_range(0,2000);
      pair<double, double> iso_range(0,5);
      shared_ptr<double> iso_threshold (new double(iso_vec[i]));
      shared_ptr<double> pt_threshold_down (new double(pt_bounds_up[i]));
      shared_ptr<double> pt_threshold_up (new double(pt_bounds_up[i+1]));


      threshold.push_back(pt_threshold_down);
      threshold.push_back(pt_threshold_up);
      threshold.push_back(iso_threshold);
 
      multi_mode.push_back(multi_cutinfo::cut_gtr);
      multi_mode.push_back(multi_cutinfo::cut_lt);
      multi_mode.push_back(multi_cutinfo::cut_lt);
      
      range_pairs.push_back(pt_range);
      range_pairs.push_back(pt_range);
      range_pairs.push_back(iso_range);

      
      stepsize.push_back(0.01);
      stepsize.push_back(0.01);
      stepsize.push_back(0.01);
      
      obsnames.push_back("pT_mu");
      obsnames.push_back("pT_mu");
      obsnames.push_back(obs_sel.added_names[i]);


      multi_cutinfo test1(obsnames, multi_mode, range_pairs, stepsize, threshold);
      start_btag_cut->add(shared_ptr<OptCut>(new MultiCut(test1, obs_sel)),CombiCut::OrCut);
      start_nobtag_cut->add(shared_ptr<OptCut>(new MultiCut(test1, obs_sel)),CombiCut::OrCut);
      start_nobtag2_cut->add(shared_ptr<OptCut>(new MultiCut(test1, obs_sel)),CombiCut::OrCut);

    }

    //start_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(relIso05cutinfo,obs_sel)),CombiCut::OrCut);
    //start_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(relIso008cutinfo,obs_sel)),CombiCut::OrCut);
    //start_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(relIso02cutinfo,obs_sel)),CombiCut::OrCut);
    start_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(relIso04cutinfo,obs_sel)),CombiCut::OrCut);
    //start_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(metcutinfo, obs_sel)),CombiCut::AndCut);
    //start_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(htcutinfo, obs_sel)),CombiCut::AndCut);
    //start_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_maxcutinfo, obs_sel)),CombiCut::AndCut);
    //start_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_210maxcutinfo, obs_sel)),CombiCut::AndCut);
    start_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Btag_cutinfo, obs_sel)),CombiCut::AndCut);
    start_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Toptag_cutinfo, obs_sel)),CombiCut::AndCut);
    //start_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Chi2_10cutinfo, obs_sel)),CombiCut::AndCut);

    start_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(relIso04cutinfo,obs_sel)),CombiCut::OrCut);
    //start_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(metcutinfo, obs_sel)),CombiCut::AndCut);
    //start_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(htcutinfo, obs_sel)),CombiCut::AndCut);
    //start_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_maxcutinfo, obs_sel)),CombiCut::AndCut);
    //start_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_210maxcutinfo, obs_sel)),CombiCut::AndCut);
    start_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(NoBtag_cutinfo, obs_sel)),CombiCut::AndCut);
    start_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(NoToptag_cutinfo, obs_sel)),CombiCut::AndCut);
    //start_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Chi2_10cutinfo, obs_sel)),CombiCut::AndCut);


    start_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(relIso04cutinfo,obs_sel)),CombiCut::OrCut);
    //start_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(metcutinfo, obs_sel)),CombiCut::AndCut);
    //start_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(htcutinfo, obs_sel)),CombiCut::AndCut);
    //start_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_maxcutinfo, obs_sel)),CombiCut::AndCut);
    //start_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_210maxcutinfo, obs_sel)),CombiCut::AndCut);
    start_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(Btag_cutinfo, obs_sel)),CombiCut::AndCut);
    start_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(NoToptag_cutinfo, obs_sel)),CombiCut::AndCut);
    //start_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(Chi2_10cutinfo, obs_sel)),CombiCut::AndCut);

    // scan through cuts; each time create the theta histfile and re-calculate expected limit.
    vector<string> signal_processes;
    for(size_t i=0; i<samples.size(); ++i){
        if(samples[i].is_signal()) signal_processes.push_back(samples[i].get_theta_procname());
    }

    //"classical" analysis:
    ofstream out("result/2D_CMSTag.txt");
    cout<<"2D"<<endl;
    vector<shared_ptr<OptCut> > signal_regions;
    signal_regions.push_back(twodcut_nobtag);
    signal_regions.push_back(twodcut_btag);
    //signal_regions.push_back(vanilla_cut);
    selection_multi_srs("rootfiles/theta_histos2D.root", samples, obs_sel, signal_regions);
    map<string, double> limits = expected_limits("rootfiles/theta_histos2D_rebinned.root", signal_processes);
    for(map<string, double>::const_iterator it=limits.begin(); it!=limits.end(); ++it){
      out <<  it->first << " " << it->first.substr(it->first.find("p")+1,it->first.size()-it->first.find("p")-1)  << " " << it->second << endl;
    }
    out.close();


    signal_regions.clear();
    signal_regions.push_back(twoD_btag_cut);
    signal_regions.push_back(twoD_nobtag_cut);
    signal_regions.push_back(twoD_nobtag2_cut);
    out.open("result/test1.txt");
    cout<<"test"<<endl;
    selection_multi_srs("rootfiles/theta_histos_test1.root", samples, obs_sel, signal_regions);
    limits = expected_limits("rootfiles/theta_histos_test1_rebinned.root", signal_processes);
    for(map<string, double>::const_iterator it=limits.begin(); it!=limits.end(); ++it){
      out <<  it->first << " " << it->first.substr(it->first.find("p")+1,it->first.size()-it->first.find("p")-1)  << " " << it->second << endl; 
    }
    out.close();

    //ofstream comp_out("result/3D.txt", ios::out | ios::app );
    ofstream comp_out("result/testscan.txt", ios::out | ios::app );
    vector<double> best_comb;

    for(int q = 0; q<1; ++q){
      double chi2_const = 45;
      shared_ptr<double> Chi2_cut_const(new double(chi2_const));
      cutinfo Chi2_cutinfo_const("chi2", cutinfo::cut_lt, make_pair(0,1000), 0.5,Chi2_cut_const);
      
      
      for(int i= 0; i<1; ++i){
	
	double chi2_tag = 50;
	shared_ptr<double> Chi2_cut_tagvary(new double(chi2_tag));
	cutinfo Chi2_cutinfo_tagvary("chi2", cutinfo::cut_lt, make_pair(0,1000), 0.5,Chi2_cut_tagvary);
	
	for(int l= 0; l<1; ++l){
	  
	  double chi2_notag = 30;
	  shared_ptr<double> Chi2_cut_notagvary(new double(chi2_notag));
	  cutinfo Chi2_cutinfo_notagvary("chi2", cutinfo::cut_lt, make_pair(0,1000), 0.5,Chi2_cut_notagvary);
	
	  
	  for(int p= 0; p<1; ++p){
	    double met_num = 50;
	    shared_ptr<double> Met_cut_vary(new double(met_num));
	    cutinfo Met_cutinfo_vary("met_pt", cutinfo::cut_gtr, make_pair(0,1000), 0.5,Met_cut_vary);
	    
	    for(int m= 0; m<1; ++m){
	      double jet_num = 230;
	      shared_ptr<double> Jet_cut_vary(new double(jet_num));
	      cutinfo Jet_cutinfo_vary("Jet_pt_max", cutinfo::cut_gtr, make_pair(0,1000), 0.5,Jet_cut_vary);
	      
	      for(int k= 0; k<1; ++k){
		double ht_num = 130;
		shared_ptr<double> HT_cut_vary(new double(ht_num));
		cutinfo HT_cutinfo_vary("HTl", cutinfo::cut_gtr, make_pair(0,1000), 5,HT_cut_vary);
		
		
		start_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Chi2_cutinfo_tagvary, obs_sel)),CombiCut::AndCut);
		start_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Chi2_cutinfo_notagvary, obs_sel)),CombiCut::AndCut);
		start_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(Chi2_cutinfo_const, obs_sel)),CombiCut::AndCut);
		
		start_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Met_cutinfo_vary, obs_sel)),CombiCut::AndCut);
		start_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Met_cutinfo_vary, obs_sel)),CombiCut::AndCut);
		start_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(Met_cutinfo_vary, obs_sel)),CombiCut::AndCut);
		
		start_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_cutinfo_vary, obs_sel)),CombiCut::AndCut);
		start_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_cutinfo_vary, obs_sel)),CombiCut::AndCut);
		start_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_cutinfo_vary, obs_sel)),CombiCut::AndCut);
		
		start_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(HT_cutinfo_vary, obs_sel)),CombiCut::AndCut);
		start_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(HT_cutinfo_vary, obs_sel)),CombiCut::AndCut);
		start_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(HT_cutinfo_vary, obs_sel)),CombiCut::AndCut);
		
   	
		signal_regions.clear();
		
		signal_regions.push_back(start_btag_cut);
		signal_regions.push_back(start_nobtag_cut);
		signal_regions.push_back(start_nobtag2_cut);
      
		cout<<"chi2_tag: "<<chi2_tag<<" chi2_notag: "<<chi2_notag<<" chi2_notop_btag: "<<chi2_const <<" met: "<<met_num<<" Jet_pT: "<<jet_num<<" HT: "<<ht_num<<endl;
		
		vector<double> compare_limits;
		
		selection_multi_srs("rootfiles/theta_histos_global_opt.root", samples, obs_sel, signal_regions);
		limits = expected_limits("rootfiles/theta_histos_global_opt_rebinned.root", signal_processes);
		for(map<string, double>::const_iterator it=limits.begin(); it!=limits.end(); ++it){
		  compare_limits.push_back(it->second);
		}
		
		
		
		if(best_comb.size()==0){
		  best_comb.swap(compare_limits);
		  out.open("result/test2.txt");
		  out<<"#chi2 tag: "<<chi2_tag<<" chi2 notag: "<<chi2_notag<<" met: "<<met_num<<" Jet pT: "<<jet_num<<" HT: "<<ht_num<<endl;
		  for(map<string, double>::const_iterator it=limits.begin(); it!=limits.end(); ++it)
		  out <<  it->first << " " << it->first.substr(it->first.find("p")+1,it->first.size()-it->first.find("p")-1)  << " " << it->second << endl;
		  out.close();
		}
		else{
		  double status=1;
		
		  for(int u =4; u<compare_limits.size()-1;++u){
		    status *=  best_comb[u]/compare_limits[u];
		  }
		  
		  if(status>1){
		    out.open("result/test2.txt");
		    out<<"chi2_tag: "<<chi2_tag<<" chi2_notag: "<<chi2_notag<<" chi2_notop_btag: "<<chi2_const <<" met: "<<met_num<<" Jet_pT: "<<jet_num<<" HT: "<<ht_num<<" discr "<< status<<endl;
		    for(map<string, double>::const_iterator it=limits.begin(); it!=limits.end(); ++it)
		      out <<  it->first << " " << it->first.substr(it->first.find("p")+1,it->first.size()-it->first.find("p")-1)  << " " << it->second << endl;
		    out.close();
		    best_comb.swap(compare_limits);
		  }
		  
		  
		  comp_out<<"chi2_tag: "<<chi2_tag<<" chi2_notag: "<<chi2_notag<<" chi2_notop_btag: "<<chi2_const <<" met: "<<met_num<<" Jet_pT: "<<jet_num<<" HT: "<<ht_num<<" discr "<< status<<endl;
		  

		}
		
		start_btag_cut->delete_last(4);
		start_nobtag_cut->delete_last(4);
		start_nobtag2_cut->delete_last(4);
	
	      }
	    }
	  }
	}
      }
      //out.close();
      
      
    }

    comp_out.close();




    vector<double> first_sum;
    double sum_old =0;

    signal_regions.clear();
    signal_regions.push_back(start_btag_cut);
    signal_regions.push_back(start_nobtag_cut);
    //signal_regions.push_back(start_nobtag2_cut);
    cout<<"starting"<<endl;
    out.open("result/starting_cut_test.txt");
    selection_multi_srs("rootfiles/theta_histos_starting_cut.root", samples, obs_sel, signal_regions);
    limits = expected_limits("rootfiles/theta_histos_starting_cut_rebinned.root", signal_processes);
    for(map<string, double>::const_iterator it=limits.begin(); it!=limits.end(); ++it){
      out <<  it->first << " " << it->first.substr(it->first.find("p")+1,it->first.size()-it->first.find("p")-1)  << " " << it->second << endl;
      first_sum.push_back(it->second);
      sum_old+=it->second;
    }
    out.close();
    

    vector<double> num_pass;
    for(unsigned int p = 0; p<pt_bounds_up.size()-1; ++p){
      num_pass.push_back((start_btag_cut->cutinfo())[p]->passed_events());
    }
    
  
    
    ofstream bound("result/bounds.txt");
    out.open("result/vary_cut.txt");
    int loops = 0;

    for(int i = 0; i<loops; ++i){

      vector<double> sum_limit;
      vector<double> num_pass_new;
      double sum_new = 0;
      shared_ptr<CombiCut> vary_btag_cut(new CombiCut());
      shared_ptr<CombiCut> vary_nobtag_cut(new CombiCut());
      shared_ptr<CombiCut> vary_nobtag2_cut(new CombiCut());

      shared_ptr<CombiCut> to_count(new CombiCut());

      vector<double> pt_bounds_upnew( i%2==0 ? vary_bounds(pt_bounds_up,pt_bounds_down,0, 50, 0,2000): pt_bounds_up);  
      vector<double> pt_bounds_downnew( i%2==0 ? pt_bounds_down : vary_bounds(pt_bounds_up,pt_bounds_down,1, 50, 0,2000));
      vector<double> iso_vec_new(iso_vec);

      for(unsigned int i= 0; i<obs_sel.added_names.size()-2;++i){
	
	vector<shared_ptr<double> > threshold;
	vector<multi_cutinfo::e_mode> multi_mode;
	vector<pair<double, double> > range_pairs;
	vector<double> stepsize;
	vector<string> obsnames;

	pair<double, double> pt_range(0,2000);
	pair<double, double> iso_range(0,5);
	shared_ptr<double> iso_threshold (new double(iso_vec_new[i]));
	shared_ptr<double> pt_threshold_down (new double(pt_bounds_downnew[i]));
	shared_ptr<double> pt_threshold_up (new double(pt_bounds_upnew[i+1]));

	threshold.push_back(pt_threshold_down);
	threshold.push_back(pt_threshold_up);
	threshold.push_back(iso_threshold);
	
	multi_mode.push_back(multi_cutinfo::cut_gtr);
	multi_mode.push_back(multi_cutinfo::cut_lt);
	multi_mode.push_back(multi_cutinfo::cut_lt);


	range_pairs.push_back(pt_range);
	range_pairs.push_back(pt_range);
	range_pairs.push_back(iso_range);


	stepsize.push_back(1);
	stepsize.push_back(1);
	stepsize.push_back(0.01);


	obsnames.push_back("pT_mu");
	obsnames.push_back("pT_mu");
	//obsnames.push_back("HT");
	//obsnames.push_back("HT");
	obsnames.push_back(obs_sel.added_names[i]);


	multi_cutinfo test1(obsnames, multi_mode, range_pairs, stepsize, threshold);
	
	vary_btag_cut->add(shared_ptr<OptCut>(new MultiCut(test1, obs_sel)),CombiCut::OrCut);
	vary_nobtag_cut->add(shared_ptr<OptCut>(new MultiCut(test1, obs_sel)),CombiCut::OrCut);
	vary_nobtag_cut->add(shared_ptr<OptCut>(new MultiCut(test1, obs_sel)),CombiCut::OrCut);
      }

      //vary_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(relIso05cutinfo,obs_sel)),CombiCut::OrCut);
      vary_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(relIso04cutinfo,obs_sel)),CombiCut::OrCut);
      vary_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(metcutinfo, obs_sel)),CombiCut::AndCut);
      vary_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(htcutinfo, obs_sel)),CombiCut::AndCut);
      vary_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_210maxcutinfo, obs_sel)),CombiCut::AndCut);
      //vary_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Btag_cutinfo, obs_sel)),CombiCut::AndCut);
      vary_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Toptag_cutinfo, obs_sel)),CombiCut::AndCut);
      vary_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Chi2_45cutinfo, obs_sel)),CombiCut::AndCut);


      vary_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(relIso04cutinfo,obs_sel)),CombiCut::OrCut);
      vary_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(metcutinfo, obs_sel)),CombiCut::AndCut);
      vary_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(htcutinfo, obs_sel)),CombiCut::AndCut);
      vary_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_210maxcutinfo, obs_sel)),CombiCut::AndCut);
      vary_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(NoBtag_cutinfo, obs_sel)),CombiCut::AndCut);
      //vary_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Chi2_cutinfo, obs_sel)),CombiCut::AndCut);
      vary_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(NoToptag_cutinfo, obs_sel)),CombiCut::AndCut);
      vary_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Chi2_45cutinfo, obs_sel)),CombiCut::AndCut);


      vary_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(relIso04cutinfo,obs_sel)),CombiCut::OrCut);
      vary_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(metcutinfo, obs_sel)),CombiCut::AndCut);
      vary_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(htcutinfo, obs_sel)),CombiCut::AndCut);
      vary_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_210maxcutinfo, obs_sel)),CombiCut::AndCut);
      vary_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(NoBtag_cutinfo, obs_sel)),CombiCut::AndCut);
      //vary_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(Chi2_cutinfo, obs_sel)),CombiCut::AndCut);
      vary_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(Toptag_cutinfo, obs_sel)),CombiCut::AndCut);
      vary_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(Chi2_45cutinfo, obs_sel)),CombiCut::AndCut);

      signal_regions.clear();
      signal_regions.push_back(vary_btag_cut);
      signal_regions.push_back(vary_nobtag_cut);
      signal_regions.push_back(vary_nobtag2_cut);
      //signal_regions.push_back(to_count);
      selection_multi_srs("theta_histos_linear_cut.root", samples, obs_sel, signal_regions);
      limits = expected_limits("theta_histos_linear_cut_rebinned.root", signal_processes);
      for(map<string, double>::const_iterator it=limits.begin(); it!=limits.end(); ++it){
	out <<  it->first << " " << it->first.substr(it->first.find("p")+1,it->first.size()-it->first.find("p")-1)  << " " << it->second <<endl;
	sum_limit.push_back(it->second);
	sum_new+=it->second;

      }

      for(unsigned int p = 0; p<pt_bounds_upnew.size()-1; ++p){
	  bound<<pt_bounds_downnew[p]<<" "<<pt_bounds_upnew[p+1]<<" "<< iso_vec_new[p]<<" "<<((vary_btag_cut->cutinfo())[p])->passed_events()<<endl;
	  num_pass_new.push_back(((vary_btag_cut->cutinfo())[p])->passed_events());
      }

      double improv_weight = 1;



      for(int u =4; u<sum_limit.size()-1;++u){
	improv_weight *= first_sum[u]/sum_limit[u];
	//cout<<"sum: "<<first_sum[u]<<" "<<sum_limit[u]<<endl;

      }
      
      //improv_weight = first_sum[3]/sum_limit[3];

      cout<<"step: "<<i <<endl;
      
      cout<<improv_weight<<endl;
      cout<<sum_new<<" ? "<<sum_old<<endl;

      if(improv_weight > 1){
      //if(sum_new<sum_old){
	pt_bounds_up.swap(pt_bounds_upnew);
	pt_bounds_down.swap(pt_bounds_downnew);
	iso_vec.swap(iso_vec_new);
	first_sum.swap(sum_limit);
	num_pass.swap(num_pass_new);
	sum_old = sum_new;

	ofstream best_limits("result/best_limits.txt");
	for(map<string, double>::const_iterator it=limits.begin(); it!=limits.end(); ++it)
	  best_limits<<  it->first << " " << it->first.substr(it->first.find("p")+1,it->first.size()-it->first.find("p")-1)  << " " << it->second <<endl;
	best_limits.close();

	ofstream best_bounds("result/best_bounds.txt");
	for(unsigned int p = 0; p<pt_bounds_up.size()-1; ++p){
	  best_bounds<<pt_bounds_down[p]<<" "<<pt_bounds_up[p+1]<<" "<< iso_vec[p]<<" "<< relIso_vec[p]  <<" "<< num_pass[p] <<endl;
	}
	best_bounds.close();

      }

    

      bound<<"Limits Weight: "<<improv_weight <<endl;
      bound<<"----"<<endl;
      out<<"---"<<endl;
    }
    
    out.close();
    bound.close();
       
    vector<double> relIso_scan;
    relIso_scan.push_back(0.2);
   
    relIso_scan.push_back(0.18);
    relIso_scan.push_back(0.16);
    relIso_scan.push_back(0.14);
    relIso_scan.push_back(0.12);
    
    relIso_scan.push_back(0.1);
    relIso_scan.push_back(0.08);
    relIso_scan.push_back(0.06);
    relIso_scan.push_back(0.04);
    relIso_scan.push_back(0.02);

    relIso_scan.push_back(0.4);

    out.open("result/scan.txt");
    ofstream numbers("result/rel_numbers.txt");
    
    vector<shared_ptr<double> > growing_threshold;
    vector<multi_cutinfo::e_mode> growing_multi_mode;
    vector<pair<double, double> > growing_range_pairs;
    vector<double> growing_stepsize;
    vector<string> growing_obsnames;

    for(int p=0; p<0; ++p){
      
      double radius = 0;
      double sum = 1000;
     
      double total_num=0;
      double pass_num=0;

      shared_ptr<CombiCut> vary_cut(new CombiCut());
      shared_ptr<CombiCut> to_count(new CombiCut());    
	
      pair<double, double> pt_range(0,2000);
      pair<double, double> iso_range(0,5);
      shared_ptr<double> iso_threshold (new double(0.2));

      double pt_down = 200*p;
      double pt_up = 200+200*p;

      shared_ptr<double> pt_threshold_down (new double(pt_down));
      shared_ptr<double> pt_threshold_up (new double(pt_up));

      string obs_result ("relIso04");

      for(unsigned int i= 0; i<obs_sel.added_names.size()-2;++i){


	double sum_new = 0;
	
	vector<shared_ptr<double> > threshold;
	vector<multi_cutinfo::e_mode> multi_mode;
	vector<pair<double, double> > range_pairs;
	vector<double> stepsize;
	vector<string> obsnames;
	
	threshold.push_back(pt_threshold_down);
	threshold.push_back(pt_threshold_up);
	threshold.push_back(iso_threshold);
	
	multi_mode.push_back(multi_cutinfo::cut_gtr);
	multi_mode.push_back(multi_cutinfo::cut_lt);
	multi_mode.push_back(multi_cutinfo::cut_lt);
	
	range_pairs.push_back(pt_range);
	range_pairs.push_back(pt_range);
	range_pairs.push_back(iso_range);
	
	stepsize.push_back(1);
	stepsize.push_back(1);
	stepsize.push_back(0.01);
	
	obsnames.push_back("pT_mu");
	obsnames.push_back("pT_mu");
	obsnames.push_back(obs_sel.added_names[i]);
	
	multi_cutinfo test1(obsnames, multi_mode, range_pairs, stepsize, threshold);
	multi_cutinfo adding_cuts(growing_obsnames, growing_multi_mode, growing_range_pairs, growing_stepsize, growing_threshold);

	vary_cut->add(shared_ptr<OptCut>(new MultiCut(test1, obs_sel)),CombiCut::OrCut);

	if(i==0){
	  threshold.pop_back();
	  stepsize.pop_back();
	  multi_mode.pop_back();
	  range_pairs.pop_back();
	  obsnames.pop_back();
	  
	  multi_cutinfo count_me(obsnames, multi_mode, range_pairs, stepsize, threshold);
	  to_count->add(shared_ptr<OptCut>(new MultiCut(count_me, obs_sel)),CombiCut::OrCut);
	}

	vary_cut->add(shared_ptr<OptCut>(new CutInfoCut(relIso04cutinfo,obs_sel)),CombiCut::OrCut);
	  
	vary_cut->add(shared_ptr<OptCut>(new CutInfoCut(metcutinfo, obs_sel)),CombiCut::AndCut);
	vary_cut->add(shared_ptr<OptCut>(new CutInfoCut(htcutinfo, obs_sel)),CombiCut::AndCut);
	vary_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_maxcutinfo, obs_sel)),CombiCut::AndCut);
   
	signal_regions.clear();
	signal_regions.push_back(vary_cut);
	if(i==0)signal_regions.push_back(to_count);
	selection_multi_srs("theta_histos_linear_cut.root", samples, obs_sel, signal_regions);
	limits = expected_limits("theta_histos_linear_cut_rebinned.root", signal_processes);
        double c = 0;
	for(map<string, double>::const_iterator it=limits.begin(); it!=limits.end(); ++it){
          c++; 
	  sum_new+=it->second*(10-c);
	}
	
	cout<<sum_new<<" "<<sum<<" "<<obs_sel.added_names[i] <<" "<<pt_up <<endl;

	if(sum_new<sum){
	  sum = sum_new;
	  total_num=((to_count->cutinfo())[0])->passed_events();
	  pass_num=((vary_cut->cutinfo())[0])->passed_events();
	  radius = relIso_scan[i];
	  obs_result = obs_sel.added_names[i];

	}
	if(i==0) numbers<<((to_count->cutinfo())[0])->passed_events()<<" "<<((vary_cut->cutinfo())[0])->passed_events();
	numbers<<" "<<((vary_cut->cutinfo())[0])->passed_events();

	vary_cut->delete_last(5);

      }
      
      out<<radius<<" "<<*iso_threshold<<" "<<*pt_threshold_down<<" "<<*pt_threshold_up<<" "<<total_num<<" "<<pass_num<<endl;
      //cout<<radius<<" "<<*iso_threshold<<" "<<*pt_threshold_down<<" "<<*pt_threshold_up<<" "<<total_num<<" "<<pass_num<<endl;
      //numbers<<" "<<*pt_threshold_down<<" "<<*pt_threshold_up<<endl;

      growing_threshold.push_back(pt_threshold_down);
      growing_threshold.push_back(pt_threshold_up);
      growing_threshold.push_back(iso_threshold);
	
      growing_multi_mode.push_back(multi_cutinfo::cut_gtr);
      growing_multi_mode.push_back(multi_cutinfo::cut_lt);
      growing_multi_mode.push_back(multi_cutinfo::cut_lt);
      
      growing_range_pairs.push_back(pt_range);
      growing_range_pairs.push_back(pt_range);
      growing_range_pairs.push_back(iso_range);
	
      growing_stepsize.push_back(1);
      growing_stepsize.push_back(1);
      growing_stepsize.push_back(0.01);
	
      growing_obsnames.push_back("pT_mu");
      growing_obsnames.push_back("pT_mu");
      growing_obsnames.push_back(obs_result);
  
      double print_size = growing_obsnames.size()/3;
   
      for(int n= 0; n<(int)print_size;++n)
 	cout<<growing_obsnames[2+3*n]<<endl; 
    }

    out.close();
    numbers.close(); 

    return 0;
}

int main(){
    try{
      return run();
    }
    catch(const string & s){
      cerr << "error: " << s << endl;
    }
    return 1;
}
