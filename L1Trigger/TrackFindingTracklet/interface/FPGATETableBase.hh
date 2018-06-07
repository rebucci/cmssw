#ifndef FPGATETABLEBASE_H
#define FPGATETABLEBASE_H

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <vector>


using namespace std;

class FPGATETableBase{

public:

  FPGATETableBase() {
   
  }

  ~FPGATETableBase() {

  }


  void writeVMTable(std::string name) {
    
    ofstream out;
    out.open(name.c_str());
    for(unsigned int i=0;i<table_.size();i++){
      out << i << " "  <<table_[i]<<endl;
    }
    out.close();
  }


protected:

  vector<int> table_;

};



#endif



