#include "PythonWrapper.h" 

using namespace std;
namespace io = boost::iostreams;


ProgramWrapper::ProgramWrapper(const char * file, char * const argv[]){
  int input[2], output[2];
  pipe(input);
  pipe(output);
  //cout << "creating child ..." << endl;
  if((child_pid = fork())){
    //cout << "created child " << child_pid << endl;
    close(input[0]);
    close(output[1]);
    childs_in = new io::file_descriptor_sink(input[1], io::close_handle);
    childs_in_buffer.open(*childs_in);
    childs_out = new io::file_descriptor_source(output[0], io::close_handle);
    childs_out_buffer.open(*childs_out);
  }
  else{
    close(input[1]);
    close(output[0]);
    dup2(input[0], STDIN_FILENO);
    dup2(output[1], STDOUT_FILENO);
    //commit suicide if parent dies:
    prctl(PR_SET_PDEATHSIG, 9);
    execvp(file, argv);
  }
}

ProgramWrapper::~ProgramWrapper(){
  //cout << "destroying child " << child_pid << endl;
  kill(child_pid, 9);
  waitpid(child_pid, 0, 0);
  delete childs_out;
  delete childs_in;
}

void ProgramWrapper::read_until_line(const string & line){
  string l;
  do{
    std::getline(childs_out_buffer, l);
  } while (l != line);
}

string ProgramWrapper::operator()(const string & cmd_in){
       childs_in_buffer << cmd_in << "\n" << flush;
       string result;
       std::getline(childs_out_buffer, result);
       return result;
   }


