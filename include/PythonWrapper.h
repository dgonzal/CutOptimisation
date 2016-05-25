#pragma once 

#include <sys/wait.h>
#include <sys/prctl.h>
#include <map>
#include <iostream>
#include <string>
#include <fstream>

#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/lexical_cast.hpp>

namespace io = boost::iostreams;

class ProgramWrapper{
 public:
   ProgramWrapper(const char* file, char* const argv[]);
   ~ProgramWrapper();
   //read from child until the string appears on a line by itself.
   //This is used to synchronize the output.
   //
   //Will block until the line was there
   void read_until_line(const std::string & line);
   std::string operator()(const std::string & cmd_in);
 private:
   io::file_descriptor_source  * childs_out;
   io::file_descriptor_sink * childs_in;
   io::stream<io::file_descriptor_source> childs_out_buffer;
   io::stream<io::file_descriptor_sink> childs_in_buffer;
   pid_t child_pid;
};

//not yet implimented!
inline bool rebinning(const std::string & theta_fname, const std::vector<std::string> & samples){
  std::cout<<" rebinning not implimented"<<std::endl;
  return true;
}

// return the expected limits by running theta
inline std::map<std::string, double> expected_limits(const std::string & modelDir, const std::string theta_dir, const std::string & theta_rootfile, const std::string & backgrounds ="", const std::string & binningDir ="", const std::string & id=""){
  std::string anaName = "analysis"+id+".py";
  std::map<std::string, double> result;
  std::ofstream ta(anaName);
  if(!binningDir.empty() && !backgrounds.empty()){
    std::size_t position_appendix = theta_rootfile.find(".root");
    ta << "execfile(\""<<binningDir<<"\")"<<std::endl;
    ta << "binFile(0.3,'"+theta_rootfile+"', 'M_{B} [GeV/c^{2}]',["+backgrounds+"])"<<std::endl;
    ta << "execfile(\""<<modelDir <<"\")" << std::endl;
    ta << "model(\"" << theta_rootfile.substr(0,position_appendix)+"_rebinned.root" << "\")" << std::endl;
  }
  else{
    ta << "execfile(\""<<modelDir <<"\")" << std::endl;
    ta << "model(\"" << theta_rootfile << "\")" << std::endl;
  }
  ta.close();
  //run theta-auto:
  char * const arg_char = (char *)anaName.data();
  char * const args[] = {arg_char, 0};
  ProgramWrapper theta_auto((theta_dir + "/utils2/theta-auto.py").c_str(), args);
  theta_auto.read_until_line("expected limit");
  for(size_t i=0; i<100; ++i){
    std::string nextline = theta_auto("");
    if(nextline=="end") break;
    size_t p = nextline.find(" ");
    if(p==std::string::npos){
      throw "unexpected line from theta-auto: " + nextline;
    }
    std::string process = nextline.substr(0, p);
    double limit = boost::lexical_cast<double>(nextline.substr(p+1));
    result[process] = limit;
  }
  return result;
}
