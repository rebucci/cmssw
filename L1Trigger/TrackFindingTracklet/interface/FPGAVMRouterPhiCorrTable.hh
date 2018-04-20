#ifndef FPGAVMROUTERPHICORRTABLE_H
#define FPGAVMROUTERPHICORRTABLE_H

#include "FPGATETableBase.hh"
#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <vector>


using namespace std;

class FPGAVMRouterPhiCorrTable: public FPGATETableBase{

public:

  FPGAVMRouterPhiCorrTable() {
   
  }

  ~FPGAVMRouterPhiCorrTable() {

  }


  void init(int layer,
	    int bendbits,
	    int rbits
	    ) {

    assert(bendbits==5);  //FIXME for now this is hardcoded
    
    layer_=layer;
    bendbits_=bendbits;
    rbits_=rbits;

    rbins_=(1<<rbits);
    rmin_=rmean[layer-1]-drmax;
    rmax_=rmean[layer-1]+drmax;
    dr_=2*drmax/rbins_;

    bendbins_=(1<<bendbits);
    bendmin_=-15.0/2.0;
    bendmax_=16.0/2.0;
    dbend_=(bendmax_-bendmin_)/bendbins_;

    rmean_=rmean[layer-1];

    for (int ibend=0;ibend<32;ibend++) {
      for (int irbin=0;irbin<rbins_;irbin++) {
	//int ibin=irbin+izbin*rbins_;
	int value=getphiCorrValue(ibend,irbin);
	table_.push_back(value);
      }
    }

    if (writeVMTables) {
      writeVMTable("VMPhiCorrL"+std::to_string(layer_)+".txt");
    }

    
  }

  int getphiCorrValue(int ibend, int irbin){

    double Delta=(irbin+0.5)*dr_-drmax;
    double dphi=Delta*0.5*(ibend-15)*0.009/0.18/rmean_;

    int idphi=0;
      
    if (layer_<=3) {
      idphi=dphi/kphi;
    } else {
      idphi=dphi/kphi1;
    }

    return idphi;
    
  }
  
  int lookupPhiCorr(int ibend, int rbin) {

    int index=ibend*rbins_+rbin;
    assert(index<(int)table_.size());
    return table_[index];
    
  }

  
private:


  double rmean_;

  double rmin_;
  double rmax_;

  double bendmin_;
  double bendmax_;

  double dr_;
  double dbend_;
  
  int bendbits_;
  int rbits_;

  int bendbins_;
  int rbins_;

  int layer_;
  
  
};



#endif



