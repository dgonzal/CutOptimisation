#pragma once 

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>       /* ceil */

class BatchManager{
 public:
  BatchManager(int NumberJobs_, int ScanPoints_):NumberJobs(NumberJobs_),ScanPoints(ScanPoints_){}
  //~BatchManager();
  int write_ScriptFile(std::string config, std::string dir);
  int submit_Jobs();
  void resubmit_Jobs();
  void check_Jobs();
  void set_Queue(std::string my_queue){queue = my_queue;}
  void set_Memory(int Disk =2, int Ram=8){DISK=Disk;RAM=Ram;}
  void set_Mail(std::string & Mail_, std::string & Notification_){Mail = Mail_; Notification =Notification_;}
 private:
  int NumberJobs;
  int ScanPoints;
  int DISK = 2;
  int RAM = 8;
  std::string queue = "short.q";
  std::string Stream = "SplitStream";
  std::string Mail = "daniel.gonzalez@dey.de";
  std::string Notification = "as";
};
