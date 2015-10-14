#pragma once 


/**
 ** Class to translate from string to Cut 
 ** Have to register the classes here since
 ** c++ does not have reflection and root 
 ** is not the best choice
 */

class CutFactory: AndCut{
 public:
  CutFactory();
  void translateCut(cosnt cutInfo & myCutInfo);
 private:

};
