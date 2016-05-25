//project code
#include "config_parser.h"
#include "sample.h"
#include "histogram_manager.h"
#include "CutManager.h"
#include "Cut.h"
#include "PythonWrapper.h"
#include "BatchManager.h"

//try some multithreading
#include <thread>
#include <future>

//boost libraries in action
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/filesystem.hpp>
//#include <boost/progress.hpp>
#include "boost/program_options.hpp" 
#include <boost/algorithm/string/predicate.hpp>
//#include <boost/timer/timer.hpp>
#include <boost/progress.hpp>
#ifdef  AUTO
#include <boost/timer/timer.hpp>
#endif     //----  AUTO  -----

//std libs I'm using
#include <iostream>
#include <fstream>

using namespace std;

bool run_theta(configInfo & config, vector<sample> & samples, CutManager & ScanCuts, string & result, string & id);
bool print_permutation(CutManager & ScanCuts, string result);
sample thread_function(const std::string & fname_pattern, const std::vector<std::string> & observables_);

int main(int argc, char* argv[]){
#ifdef  AUTO
    boost::timer::auto_cpu_timer t;
#endif     //----  AUTO  -----

  /**s
   ** TO DO: write help command & securit checks
   **
   ** Define and parse the program options
   ** 
   */ 
  
  namespace po = boost::program_options; 
  po::options_description desc("Options"); 
  desc.add_options() 
    ("help,h", "Help message. Usage is /bin/run [options] Config.ini")
    ("cutrange,r",po::value<std::vector<int> >()->multitoken(),"Define which of the i scans are perfomed. One argument is just a specific cut configuration, two give a range and more are treated as vector of cut configurations. Not yet implemented.")
    ("submit,s",po::value<int> (),"Prepare and submit n jobs to the batch.")
    ("result", po::value<string>(),"Name of the result file")
    ("id",po::value<string>()," For splitting jobs. Every job need its utilities. Like analsis.py.")
    ("ConfigFile",po::value<string>()->required(),"Don't forget to specify the Config.ini File.")
    ("testrun",po::bool_switch()->default_value(false),"Just to see how many Scanpoints are beeing performed.");
  

  //po::positional_options_description positionalOptions; 
  //  positionalOptions.add("add", -1); 
  po::positional_options_description p;
  p.add("ConfigFile", -1);
  po::variables_map vm; 
  //po::store(po::parse_command_line(argc, argv, desc),vm);
  po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
  po::notify(vm);
  

  if(vm.count("help")){
    cout<<desc<<endl;
    return 0;
  }

  vector<int> jobs_ids;
  int number_batchJobs=-1;
  string id = "";
  if(vm.count("id"))
    id =  vm["id"].as<string>();
  string result = "result"+id+".txt";
  if(vm.count("cutrange"))
    jobs_ids = vm["cutrange"].as<vector<int>>();
  if(vm.count("result"))
    result = vm["result"].as<string>();
  string configfileName = vm["ConfigFile"].as<string>();
  //for(auto id : jobs_ids)
  //  cout<<"ids "<<id<<endl;
  /*
  if(!boost::algorithm::ends_with(argv[argc-1],".ini")){
  cerr<<"Last argument does not end with ini. Does not seem to be a Config.ini File."<<endl;
  return 1;
  } */   
  auto cwd = boost::filesystem::current_path();
  bool debug = false;
  boost::property_tree::ptree pt;
  //boost::property_tree::ini_parser::read_ini(argv[argc-1],pt);
  boost::property_tree::ini_parser::read_ini(configfileName,pt);
  configInfo currentConfig(pt);

  if(debug)currentConfig.printConfig();
  if(!boost::filesystem::exists(currentConfig.get_Workdir()+id))
    boost::filesystem::create_directories(currentConfig.get_Workdir()+id);
  else
    cout<<"Workdir already exists, I will ignore this and write stuff in there for now"<<endl;
  chdir((currentConfig.get_Workdir()+id).c_str());
 
  //initialise Cuts
  CutManager ScanCuts;
  for(auto & cut : currentConfig.get_Cuts())
    ScanCuts.add_Cut(cut.nick, cut.cutname(),cut.cutInput());
  ScanCuts.set_PermutationVector(jobs_ids);

  cout << "Number of Points that will be scanned "<<ScanCuts.get_NumberScanPoints()<<endl;

  if(vm["testrun"].as<bool>()){
    print_permutation(ScanCuts,"Permuatition_Test.txt");
    return 0;
  }
  if(ScanCuts.get_NumberScanPoints()>1000 && ScanCuts.get_NumberScanPoints()< 10000)
    cout<<"Why do I have to work so much?"<<endl;
  else if(ScanCuts.get_NumberScanPoints()>10000)
    cout<<"You are killing me!"<<endl;
 
  if(vm.count("submit")){
    number_batchJobs = vm["submit"].as<int>();
    cout<<"Number of batch Jobs "<<number_batchJobs<<endl;
    BatchManager batchmanager(number_batchJobs,ScanCuts.get_NumberScanPoints());
    batchmanager.write_ScriptFile(configfileName,cwd.string());
    batchmanager.submit_Jobs();
    return 0;
  }

 
  //load sample trees  
  vector<sample> sample_vec;
  vector<string> observables = currentConfig.get_Variables();

  //serial version
  for(auto & data :  currentConfig.get_Samples()){
    sample mysample(currentConfig.get_Path()+data.dir,observables);
    sample_vec.push_back(mysample);
  } 
  /*
  //load sample trees multithread version
  //this is somehow difficult because of ROOT argh!!!!
  //Error Msg -> Fatal in <TClass::Init>: gInterpreter not initialized
  //You have to use the TThread stuff :'(
  //int number_of_threads = 4;
  std::vector<std::future<sample>> future_results;
  //std::cout <<"Main thread " <<std::this_thread::get_id() <<endl;
  for(auto & data :  currentConfig.get_Samples()){
    //cout<<"creating thread for "<<data.dir<<endl;
    future_results.push_back(async(std::launch::async,thread_function,currentConfig.get_Path()+data.dir,observables));
  }
  for(auto & result :  future_results){
    sample_vec.push_back(result.get());
  } 
  */

  string emptystring = "" ;
  run_theta(currentConfig,sample_vec,ScanCuts,result,emptystring);
  return 0;
}

struct resultInfo{
  string cuts;
  vector<double> limits;
  double compareValue;
};

sample thread_function(const std::string & fname_pattern, const std::vector<std::string> & observables_){
  //std::cout <<"Daughter thread "<< std::this_thread::get_id() <<endl;
  sample mysample(fname_pattern,observables_);
  return mysample;
}

bool run_theta(configInfo & config, vector<sample> & samples, CutManager & ScanCuts, string & result, string & id){
  //Get names of the signal samples
  ofstream out(result);
  out<<"CutName ";
  for(auto & sampleName : config.get_SampleInformation(2)){
    out<<sampleName<<" ";
  }
  out<<endl;
  //Get names of the background process, to rebin them later
  string backgrounds = "";
  bool first = true;
  for(auto & sampleName : config.get_SampleInformation(3)){
    if(first)
      backgrounds = backgrounds+"'"+sampleName+"'";
    else
      backgrounds = backgrounds+",'"+sampleName+"'";
    first = false;
  }
  vector<resultInfo> best_results;
  unsigned int numberOfCuts = ScanCuts.number_Cuts();
  HistoManager ThetaHists;
  unsigned int numberOfObsHist = config.get_PlotInfo().size();
  for(auto & sampleName : config.get_SampleInformation(1)){
    for(auto PlotInfo : config.get_PlotInfo())
      ThetaHists.create_Hists(PlotInfo ,sampleName,config.get_Sys());
  }
  cout<<" number of scanned cuts "<<numberOfCuts<<" total number of points "<<ScanCuts.get_NumberScanPoints()<<endl;
  boost::progress_display show_progress(ScanCuts.get_NumberScanPoints());
  for(int p =0; p<ScanCuts.get_NumberScanPoints(); p++){
    ScanCuts.permute_Cuts(p);
    int i =0;
    for(auto & sample : samples){
      for(int m = 0; m< sample.get_NumberEvents(); m++){
	bool passCuts = true;
	//int frac_samp = int(double(m)/double(sample.get_NumberEvents())*100);
	//if(frac_samp%10)cout<<"\r"<<"Scan Progress [%] "<<int(double(p)/double(ScanCuts.get_NumberScanPoints())*100) <<" "<<ScanCuts.get_OuputName() <<" "<<config.get_Samples()[i].nick<<" [%]: "<<frac_samp<<flush;// " pass "<< pass << " fail "<<fail<<flush;
	//cout<<"==========================================="<<endl;
	for(unsigned int it =0; it <numberOfCuts; it++){
	  passCuts = ScanCuts.check_Cut(it, sample.get_Value(m,it+numberOfObsHist+1));
	  //cout<<passCuts<<" Value "<<sample.get_Value(m,it+numberOfObsHist+1) <<endl;
	  if(!passCuts) break;
	  if(it+1==numberOfCuts){
	    int position = 0; double obsVal;
	    for(unsigned int mp = 0; mp < numberOfObsHist;mp++){
	      obsVal = sample.get_Value(m,mp);
	      position = mp;
	      //cout<<obsVal<<endl;
	      if(obsVal >= config.get_PlotInfo()[mp].min && obsVal <= config.get_PlotInfo()[mp].max)break;
	    }
	    //cout<<"position "<<position<<" Val "<<obsVal<<" weight "<< sample.get_Value(m,numberOfObsHist)<<endl;
	    //ThetaHists.hists[i+position]->Fill(obsVal,sample.get_Value(m,numberOfObsHist));
	    ThetaHists.fill_ithHist(i+position,obsVal,sample.get_Value(m,numberOfObsHist),sample.get_SysWeight(m));
	  }
	}
      }
      i+= numberOfObsHist; 
    }
    resultInfo myInfo;
    myInfo.cuts = ScanCuts.get_OuputName();
    myInfo.compareValue =1;
    ThetaHists.write_HistsToFile(ScanCuts.get_OuputName()+".root");
    //cout<<config.get_ThetaModelDir()<<" "<< config.get_ThetaDir()<<" "<<ScanCuts.get_OuputName()+".root"<<endl;
    map<string, double> limits = expected_limits(config.get_ThetaModelDir(), config.get_ThetaDir(),ScanCuts.get_OuputName()+".root",backgrounds,config.get_Rebinning(),id);
    out<<ScanCuts.get_OuputName()<<" ";
    for(map<string, double>::const_iterator it=limits.begin(); it!=limits.end(); ++it){
      out << it->second<<" ";
      myInfo.limits.push_back(it->second);
      myInfo.compareValue *=it->second;
    }
    out<<endl;
    if(p<5) best_results.push_back(myInfo);
    else{
      int position = -1;
      double distance = 0;
      for(unsigned int ip = 0; ip < 5;ip++){
	if(myInfo.compareValue<best_results[ip].compareValue && distance < fabs(myInfo.compareValue-best_results[ip].compareValue)){
	  position = ip;
	  distance = fabs(myInfo.compareValue-best_results[ip].compareValue);
	}
      }
      if(position != -1)best_results[position]= myInfo;
    }		   
    ThetaHists.reset_Hists();
    ++show_progress;
  }
  out<<"================================================================="<<endl;
  for(auto & result: best_results){
    out<<result.cuts<<" ";
    for(auto & number : result.limits)
      out<<number<<" ";
    out<<endl;
  }
  
  out.close();
  return true;
}


bool print_permutation(CutManager & ScanCuts, string result){
  ofstream out(result);
  for(int p =0; p<ScanCuts.get_NumberScanPoints(); p++){
    ScanCuts.permute_Cuts(p);
    out<<ScanCuts.get_OuputName()<<endl;
  }
  return true;
}






/*
** Ideas:
** -> Parallel reading of samples, don't forget to block mem if writting to it
** -> Split work load into making histograms & rebinning them and running theta
**
** -> Put stuff into a class that manage the histogram production
** -> most of the stuff in run_theta is rather generic and could go in main as part of a class
**
** -> In principale one could also disintangle theta computation and hist production 
**    if something goes wrong, one could restart it easily at any given point. 
**
*/
