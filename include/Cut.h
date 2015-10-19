#pragma once 
/*
class CutInterface {
 public:
  CutInterface(){}
  virtual bool passes(double eventValue)=0;
};
*/

#include <iostream>
#include <stdlib.h>  
#include <string>

class Cut{
 public: 
  Cut(int numberOfSteps_);
  //~Cut();
  ///void next_Cut(){counter++;}
  void set_counter(int counter_){counter=counter_;}
  int get_counter(){return counter;}
  int get_numberOfSteps(){return numberOfSteps;}
  virtual const std::string get_CutValues();
  virtual bool passes(double eventValue);
 private:
  int counter;
  int numberOfSteps;
  //virtual bool operator();//look up how to implement it, see uhh2!!!
};

class CutGreater: public Cut{
 public:
  CutGreater(double step_, double min_, double max_, int numberOfSteps_);
  //~CutGreater();
  bool passes(double eventValue);
  const std::string get_CutValues();
 private:
  double step, min, max;
};

class CutSmaller: public Cut {
 public:
  CutSmaller(double step_, double min_, double max_, int numberOfSteps_);
  //~CutSmaller();
  bool passes(double eventValue);
  const std::string get_CutValues();
 private:
  double step, min, max;
};

