#include "include/root-utils.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TMethodCall.h"

#include <stdexcept>
#include <cassert>

#include <cxxabi.h>

using namespace std;


std::string demangle(const char * typename_){
    int status;
    char * demangled = abi::__cxa_demangle(typename_, 0, 0, &status);
    if(status != 0 or demangled == 0){
        throw runtime_error("could not demangle type name '" + string(typename_) + "'");
    }
    string result(demangled);
    free(demangled);
    return result;
}

TrivialTree::TrivialTree(TTree & t): tree(t){
}

TBranch * TrivialTree::get_branch(const string & branchname){
    TBranch * b = tree.GetBranch(branchname.c_str());
    if(!b){
        throw runtime_error("branch '" + branchname + "' not found");
    }
    return b;
}

void TrivialTree::SetBranchAddress(const string & bname, void * addr){
    auto it = branches.find(bname);
    if(it == branches.end()){
        branches[bname] = get_branch(bname);
    }
    branches[bname]->SetAddress(addr);
}

size_t TrivialTree::GetEntries(){
    return tree.GetEntries();
}
    
size_t TrivialTree::GetEntry(size_t i){
    size_t res = 0;
    for(auto & bname_tb : branches){
        res += bname_tb.second->GetEntry(i);
    }
    return res;
}

namespace{

  void merge(TDirectory * lhs, TDirectory * rhs){
    map<string, TKey*> lhs_keys;
    map<string, TKey*> rhs_keys;
    {
        TList * l_keys = lhs->GetListOfKeys();
        TIter next(l_keys);
        while(TObject * key_ = next()){
	  TKey * key = dynamic_cast<TKey*>(key_);
	  assert(key);
	  lhs_keys[key->GetName()] = key;
        }
    }
    {
      TList * r_keys = rhs->GetListOfKeys();
      TIter next(r_keys);
      while(TObject * key_ = next()){
	TKey * key = dynamic_cast<TKey*>(key_);
	assert(key);
	rhs_keys[key->GetName()] = key;
      }
    }
    
    if(lhs_keys.size() != rhs_keys.size()) throw runtime_error("merge: directories have not the same number of entries!");
    
    // iterate over lhs keys and merge each one with rhs:
    for(auto & lit : lhs_keys){
      auto rit = rhs_keys.find(lit.first);
      if(rit == rhs_keys.end()) throw runtime_error("merge: object '" + lit.first + "' not found in right list");
      // depending on the type, make different stuff now:
      string l_class = lit.second->GetClassName();
      // but first checki that it is the same class:
      if(l_class != rit->second->GetClassName()) throw runtime_error("merge: object '" + lit.first + "' has different type in left and right hand of merge");
      // now try to call the "Merge" method of the left hand object:
      TObject * left_object = lit.second->ReadObj();
      if(!left_object) throw runtime_error("could not read object '" + lit.first + "'");
      // if it is a TDirectory, merge recursively:
      if(TDirectory * left_tdir = dynamic_cast<TDirectory*>(left_object)){
	TDirectory * right_tdir = dynamic_cast<TDirectory*>(rit->second->ReadObj());
	assert(right_tdir);
	merge(left_tdir, right_tdir);
	continue;
      }
      else{
	// otherwise, call the 'merge' method:
	TMethodCall mergeMethod;
	mergeMethod.InitWithPrototype(left_object->IsA(), "Merge", "TCollection*" );
	if(!mergeMethod.IsValid()){
	  throw runtime_error("object '" + lit.first + "' (class '" + l_class + "') has no 'Merge' method");
	}
	TList l;
	l.Add(rit->second->ReadObj());
	mergeMethod.SetParam((Long_t)&l);
	mergeMethod.Execute(left_object);
	// remove the original TKey in the output file:
	lit.second->Delete();
	delete lit.second;
      }
    }
  }
}

void merge_rootfiles(const std::string & file1, const std::string & file2){
    TFile f1(file1.c_str(), "update");
    TFile f2(file2.c_str(), "read");
    merge(&f1, &f2);
    f1.Write();
    f1.Close();
}

