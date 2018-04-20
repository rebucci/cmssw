//This class holds the reduced VM stubs
#ifndef FPGAVMSTUBSME_H
#define FPGAVMSTUBSME_H

#include "L1TStub.hh"
#include "FPGAStub.hh"
#include "FPGAMemoryBase.hh"

using namespace std;

class FPGAVMStubsME:public FPGAMemoryBase{

public:

  FPGAVMStubsME(string name, unsigned int iSector, 
	      double phimin, double phimax):
    FPGAMemoryBase(name,iSector){
    phimin_=phimin;
    phimax_=phimax;
  }

  void addStub(std::pair<FPGAStub*,L1TStub*> stub) {
    stubs_.push_back(stub);
    if (stub.first->isBarrel()) { // barrel
      int bin=3-(stub.first->z().value()>>(stub.first->z().nbits()-MEBinsBits));
      //cout << "FPGAVMStubsME::addStub "<<bin<<" "<<stub.first->z().value()<<" "<<stub.first->z().nbits()<<endl;
      assert(bin>=0);
      assert(bin<MEBins);
      if (debug1) {
	cout << getName() << " adding stub to bin "<<bin<<endl;
      }
      binnedstubs_[bin].push_back(stub);
    }
    else { // disk 
      int ir = stub.first->r().value();
      //presumably the VMRouter can use a lookup table for the bin.
      //For now implement a simple cut
      unsigned int bin=3;
      if (stub.first->isPSmodule()){
	if (ir<50.0/kr) bin=2;
	if (ir<35.0/kr) bin=1;
	if (ir<26.0/kr) bin=0;
      }
      if (stub.first->disk().value()<0) bin+=4;
      binnedstubs_[bin].push_back(stub);
    }
  }

  unsigned int nStubs() const {return stubs_.size();}

  FPGAStub* getFPGAStub(unsigned int i) const {return stubs_[i].first;}
  L1TStub* getL1TStub(unsigned int i) const {return stubs_[i].second;}
  std::pair<FPGAStub*,L1TStub*> getStub(unsigned int i) const {return stubs_[i];}

  unsigned int nStubsBin(unsigned int bin) const {
    assert(bin<MEBins);
    return binnedstubs_[bin].size();
  }

  std::pair<FPGAStub*,L1TStub*> getStubBin(unsigned int bin, unsigned int i) const {
    assert(bin<MEBins);
    assert(i<binnedstubs_[bin].size());
    return binnedstubs_[bin][i];
  }
  
  void clean() {
    stubs_.clear();
    for (unsigned int i=0; i<NLONGVMBINS; i++){
      binnedstubs_[i].clear();
    }
  }

  void writeStubs(bool first) {

    std::string fname="MemPrints/VMStubsME/VMStubs_";
    fname+=getName();
    //get rid of duplicates
    int len = fname.size();
    if(fname[len-2]=='n'&& fname[len-1]>'1'&&fname[len-1]<='9') return;
    //
    fname+="_";
    ostringstream oss;
    oss << iSector_+1;
    if (iSector_+1<10) fname+="0";
    fname+=oss.str();
    fname+=".dat";
    if (first) {
      bx_ = 0;
      event_ = 1;
      out_.open(fname.c_str());
    }
    else
      out_.open(fname.c_str(),std::ofstream::app);

    out_ << "BX = "<<(bitset<3>)bx_ << " Event : " << event_ << endl;

    for (unsigned int i=0;i<NLONGVMBINS;i++) {
      for (unsigned int j=0; j<binnedstubs_[i].size();j++) {
        string stub = binnedstubs_[i][j].first->stubindex().str();
        out_ << hex << i << " " << j << dec << " " << stub << endl;
      }
    }
    out_.close();

    bx_++;
    event_++;
    if (bx_>7) bx_=0;

  }

private:

  double phimin_;
  double phimax_;
  std::vector<std::pair<FPGAStub*,L1TStub*> > stubs_;

  std::vector<std::pair<FPGAStub*,L1TStub*> > binnedstubs_[MEBins];

  
  
};

#endif
