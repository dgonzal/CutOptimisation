#pragma once

#include "TTree.h"
//#include "log.hpp"

#include <map>
#include <string>
#include <stdexcept>


// a simpler/well-defined Tree interface
class TrivialTree{
public:
    // note: make sure that t lives longer than this wrapper
    explicit TrivialTree(TTree & t);
    
    // t should be a null pointer; the object will be created as required and destroyed if the underlyting TTree's lifetime ends.
    template<typename T>
    void connect_object(const std::string & branchname, T*& t);
    
    // connect a 'plain' variable (float, int, etc.)
    template<typename T>
    void connect_plain(const std::string & branchname, T & t);

    // proxy to TBranch::SetAddress; use this only if connect_object and connect_plain are not working for you,
    // as this method does not perform any type checking.
    void SetBranchAddress(const std::string & bname, void * addr);
    
    size_t GetEntries();
    
    size_t GetEntry(size_t i);
    
    ~TrivialTree(){}

private:
    TBranch * get_branch(const std::string & branchname); // throws if not found.
    TTree & tree;
    std::map<std::string, TBranch*> branches;
    
};

std::string demangle(const char * typename_);


template<typename T>
void TrivialTree::connect_object(const std::string & branchname, T*&t){
    TBranch * b = get_branch(branchname);
    TClass * class_ = 0;
    EDataType type;
    if(0 != b->GetExpectedType(class_, type)){
        throw std::runtime_error("could not determine type of branch '" + branchname + "'");
    }
    if(class_ == 0){
        throw std::runtime_error("branch '" + branchname + "' is of non-class type, but connect_object was called (call connect_variable instead!)");
    }
    const type_info * ti = class_->GetTypeInfo();
    if(ti==0){
        throw std::runtime_error("no type_info for branch '" + branchname + "'");
    }
    if(*ti == typeid(T)){
      std::cout << "Branch '" << branchname << "': type checking successful" << std::endl;
    }
    else{
      std::cout << "Branch '" << branchname << "' has type '" << class_->GetName() << "' but supplied variable has type '" << demangle(typeid(T).name()) << "'" << std::endl;
        throw std::runtime_error("Branch type mismatch (see log for details)");
    }
    if(t==0) t = new T();
    SetBranchAddress(branchname, t);
}

template<typename T>
void TrivialTree::connect_plain(const std::string & branchname, T & t){
    TBranch * b = get_branch(branchname);
    TClass * class_ = 0;
    EDataType type;
    if(0 != b->GetExpectedType(class_, type)){
        throw std::runtime_error("could not determine type of branch '" + branchname + "'");
    }
    if(class_ != 0){
        throw std::runtime_error("branch '" + branchname + "' is of class type, use connect_object instead of connect_plain!");
    }
    if(type == TDataType::GetType(typeid(T))){
        std::cerr << "Branch '" << branchname << "': type checking successful" << std::endl;
    }
    else{
         std::cout << "Branch '" << branchname << "' has type '" << TDataType::GetDataType(type)->GetTypeName() << "' but supplied variable has type '" << demangle(typeid(T).name()) << "'" << std::endl;
         throw std::runtime_error("Branch type mismatch (see log for details)");
    }
    SetBranchAddress(branchname, &t);
}

// merge rootfile file2 into file1. Both files must exist and contain the same objects.
void merge_rootfiles(const std::string & file1, const std::string & file2);

