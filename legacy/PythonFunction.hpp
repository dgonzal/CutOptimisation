#ifndef PYTHONFUNCTION_HPP
#define PYTHONFUNCTION_HPP


#include <string>

#include "Python.h"
#include <boost/python.hpp>
#include <boost/filesystem.hpp>


using namespace std;


class PythonFunction{
private:
public:
  ~PythonFunction(){}

  int rebinning(char* format, const char* histo_name){
    PyObject *pName, *pModule, *pDict, *pFunc, *pCall;
    Py_Initialize();
    pName = PyString_FromString("histogram_rebinning"); 
    pModule = PyImport_Import(pName);
    if(pModule == NULL){
      PyErr_Occurred();
      PyErr_Print();
    }
    pDict = PyModule_GetDict(pModule);
    pFunc = PyDict_GetItemString(pDict, "binFile");
    if (PyCallable_Check(pFunc)){
      pCall = PyObject_CallFunction(pFunc,format,0.3,histo_name,"'M_{t#bar{t}} [GeV/c^{2}]'","['ttbar','qcd','wjets','zjets']'");
      if(pCall==NULL){
	PyErr_Occurred();
	PyErr_Print();
      }
    }
    else{
	PyErr_Occurred();
	PyErr_Print();
    }
    Py_DECREF(pFunc);
    Py_DECREF(pModule);  
    Py_DECREF(pName);
    Py_Finalize();
    return 0;
  }
};
#endif
