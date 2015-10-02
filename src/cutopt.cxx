// This is the cut optimizer for the Z'->ttbar analysis. It uses ntuples produced by highmass_tuple
// as input.
// Thanks to Jochen Ott for the code.
#include "ProgramWrapper.hpp"
#include "PythonFunction.hpp"

#include "TFile.h"
#include "TTree.h"

#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/scoped_array.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/predicate.hpp>
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
    twodcut_btag->add(shared_ptr<OptCut>(new CutInfoCut(Chi2_10cutinfo, obs_sel)));

    shared_ptr<AndCut> twodcut_nobtag(new AndCut());
    twodcut_nobtag->add(shared_ptr<OptCut>(new CutInfoCut(TwoD_cutinfo, obs_sel)));
    twodcut_nobtag->add(shared_ptr<OptCut>(new CutInfoCut(metcutinfo, obs_sel)));	     
    twodcut_nobtag->add(shared_ptr<OptCut>(new CutInfoCut(htcutinfo, obs_sel)));	     
    twodcut_nobtag->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_maxcutinfo, obs_sel)));
    twodcut_nobtag->add(shared_ptr<OptCut>(new CutInfoCut(NoBtag_cutinfo, obs_sel)));
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

    twoD_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(TwoD_cutinfo, obs_sel)),CombiCut::AndCut);
    twoD_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(metcutinfo, obs_sel)),CombiCut::AndCut);
    twoD_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(htcutinfo, obs_sel)),CombiCut::AndCut);
    twoD_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_maxcutinfo, obs_sel)),CombiCut::AndCut);
    twoD_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Btag_cutinfo, obs_sel)),CombiCut::AndCut);
    twoD_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Toptag_cutinfo, obs_sel)),CombiCut::AndCut);
    twoD_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Chi2_10cutinfo, obs_sel)),CombiCut::AndCut);

    twoD_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(TwoD_cutinfo, obs_sel)),CombiCut::AndCut);
    twoD_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(metcutinfo, obs_sel)),CombiCut::AndCut);
    twoD_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(htcutinfo, obs_sel)),CombiCut::AndCut);
    twoD_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_maxcutinfo, obs_sel)),CombiCut::AndCut);
    twoD_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(NoBtag_cutinfo, obs_sel)),CombiCut::AndCut);
    twoD_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(NoToptag_cutinfo, obs_sel)),CombiCut::AndCut);
    twoD_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Chi2_10cutinfo, obs_sel)),CombiCut::AndCut);

    twoD_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(TwoD_cutinfo, obs_sel)),CombiCut::AndCut);
    twoD_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(metcutinfo, obs_sel)),CombiCut::AndCut);
    twoD_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(htcutinfo, obs_sel)),CombiCut::AndCut);
    twoD_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(Jet_pt_maxcutinfo, obs_sel)),CombiCut::AndCut);
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

    for(unsigned int i=0; i<relIso_vec.size(); ++i){
      double limit =  29.1356/(relIso_vec[i]- 0.023111)-164.383;   
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

    start_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(relIso04cutinfo,obs_sel)),CombiCut::OrCut);
    start_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Btag_cutinfo, obs_sel)),CombiCut::AndCut);
    start_btag_cut->add(shared_ptr<OptCut>(new CutInfoCut(Toptag_cutinfo, obs_sel)),CombiCut::AndCut);
    start_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(relIso04cutinfo,obs_sel)),CombiCut::OrCut);
    start_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(NoBtag_cutinfo, obs_sel)),CombiCut::AndCut);
    start_nobtag_cut->add(shared_ptr<OptCut>(new CutInfoCut(NoToptag_cutinfo, obs_sel)),CombiCut::AndCut);
    start_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(relIso04cutinfo,obs_sel)),CombiCut::OrCut);
    start_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(Btag_cutinfo, obs_sel)),CombiCut::AndCut);
    start_nobtag2_cut->add(shared_ptr<OptCut>(new CutInfoCut(NoToptag_cutinfo, obs_sel)),CombiCut::AndCut);
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
