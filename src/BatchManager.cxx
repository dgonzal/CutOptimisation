#include "BatchManager.h"

#include <boost/filesystem.hpp>

/*
** Here the information for submission are stored! 
**
** You could specify the queue here but it is easy to 
** export QUEUE=short.q
** to get any queue you want to have!
**
** myfile << "##Set the queue you want to job to end up"<<endl;
** myfile << "##$ -v QUEUE=short.q"<<endl;
**
**
*/


using namespace std;

int BatchManager::write_ScriptFile(string config, string dir){
  int numberScanPointsPerJob = ceil(float (ScanPoints)/ float (NumberJobs));
  int max = ScanPoints;
  ofstream myfile;
  myfile.open ("cutopt.sh",fstream::out);
  myfile << "#!/bin/bash"<<endl;
  myfile << "##This is a simple example of a SGE batch script"<<endl;
  myfile << "##Use home server with scientific linux 6"<<endl; 
  myfile << "#$ -l os=sld6"<<endl; 
  myfile << "#$ -l site=hh"<<endl; 
  myfile << "#$ -cwd"<<endl;
  myfile << "##Set the queue you want to job to end up"<<endl;
  myfile << "#$ -q "<<queue<<endl;
  myfile << "##You need to set up enviroment"<<endl;
  myfile << "#$ -V"<<endl; 
  myfile << "##email Notification"<<endl;
  myfile << "#$ -m "<<Notification<<endl;
  myfile << "#$ -M "<<Mail<<endl;
  myfile << "##running in local mode with 8-12 cpu slots"<<endl;
  myfile << "##$ -pe local 8-12"<<endl;
  myfile << "##CPU memory"<<endl;
  myfile << "#$ -l h_vmem="<<RAM<<"G"<<endl;
  myfile << "##DISK memory"<<endl;
  myfile << "#$ -l h_fsize="<<DISK<<"G "<<endl;  
  myfile << "##Doing the work"<<endl;
  myfile << "a=$(("<<to_string(numberScanPointsPerJob) <<" * ${SGE_TASK_ID}))"<<endl;
  myfile << "b=$(($a + "<<to_string(numberScanPointsPerJob)<<"))"<<endl;
  myfile << "if [ $b -gt "<<to_string(max)<<" ] ;then "<<endl;
  myfile << "b="<<max<<endl; 
  myfile << "fi"<<endl; 
  myfile << "if [ $a -lt "<<to_string(max)<<" ] ;then "<<endl;
  //myfile << "cd "<<workdir<<endl;
  myfile << dir<<"/bin/run --id ${SGE_TASK_ID} -r $a $b -- "<<dir<<"/"<<config<<endl;
  myfile << "fi"<<endl; 
  myfile.close();
  return 0;
}
 
int BatchManager::submit_Jobs(){
  if(!boost::filesystem::exists(Stream))
    boost::filesystem::create_directories(Stream);
  system(("qsub -t 1-"+to_string(NumberJobs)+" -o "+Stream+"/"+" -e "+Stream+"/"+" cutopt.sh").c_str());
  return 0;
}

