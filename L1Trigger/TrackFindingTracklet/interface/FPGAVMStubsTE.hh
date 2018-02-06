//This class holds the reduced VM stubs
#ifndef FPGAVMSTUBSTE_H
#define FPGAVMSTUBSTE_H

#include "L1TStub.hh"
#include "FPGAStub.hh"
#include "FPGAMemoryBase.hh"

using namespace std;

class FPGAVMStubsTE:public FPGAMemoryBase{

public:

  FPGAVMStubsTE(string name, unsigned int iSector, 
	      double phimin, double phimax):
    FPGAMemoryBase(name,iSector){
    phimin_=phimin;
    phimax_=phimax;
    string subname=name.substr(6,2);
    layer_ = 0;
    disk_  = 0;
    phibin_=-1;
    other_=0;
    if (subname=="L1") layer_=1;
    if (subname=="L2") layer_=2;
    if (subname=="L3") layer_=3;
    if (subname=="L4") layer_=4;
    if (subname=="L5") layer_=5;
    if (subname=="L6") layer_=6;
    if (subname=="D1") disk_=1;
    if (subname=="D2") disk_=2;
    if (subname=="D3") disk_=3;
    if (subname=="D4") disk_=4;
    if (subname=="D5") disk_=5;
    if (layer_==0&&disk_==0) {
      cout << name<<" subname = "<<subname<<" "<<layer_<<" "<<disk_<<endl;
    }
    assert((layer_!=0)||(disk_!=0));
    overlap_=false;
    subname=name.substr(11,1);
    if (subname=="X") overlap_=true;
    if (subname=="Y") overlap_=true;
    if (subname=="W") overlap_=true;
    if (subname=="Q") overlap_=true;

    subname=name.substr(12,2);
    phibin_=subname[0]-'0';
    if (subname[1]!='n') {
      phibin_=10*phibin_+(subname[1]-'0');
    }
    
  }
  
  bool addStub(std::pair<FPGAStub*,L1TStub*> stub) {

    int binlookup=stub.first->getVMBits().value();
    assert(binlookup>=0);
    int bin=(binlookup/8);

    if (overlap_) {
	if (disk_==1) {
          bool negdisk=stub.first->disk().value()<0.0;
	  assert(bin<4);
	  if (negdisk) bin+=4;
	  stubsbinned_[bin].push_back(stub);
	  if (debug1) cout << getName()<<" Stub with lookup = "<<binlookup
			   <<" in disk = "<<disk_<<"  in bin = "<<bin<<endl;
	}
    } else {
      if (stub.first->isBarrel()){
	if (layer_%2==0) {
	  stubsbinned_[bin].push_back(stub);
	}
	
      } else {

        bool negdisk=stub.first->disk().value()<0.0;

	if (disk_%2==0) {
	  assert(bin<4);
	  if (negdisk) bin+=4;
	  stubsbinned_[bin].push_back(stub);
	}
        	
      }
    }
    stubs_.push_back(stub);
    return true;
  }

  unsigned int nStubs() const {return stubs_.size();}
  unsigned int nStubsBinned(unsigned int bin) const {return stubsbinned_[bin].size();}

  FPGAStub* getFPGAStub(unsigned int i) const {return stubs_[i].first;}
  L1TStub* getL1TStub(unsigned int i) const {return stubs_[i].second;}
  std::pair<FPGAStub*,L1TStub*> getStub(unsigned int i) const {return stubs_[i];}

  FPGAStub* getFPGAStubBinned(unsigned int bin, unsigned int i) const {return stubsbinned_[bin][i].first;}
  L1TStub* getL1TStubBinned(unsigned int bin, unsigned int i) const {return stubsbinned_[bin][i].second;}
  std::pair<FPGAStub*,L1TStub*> getStubBinned(unsigned int bin,unsigned int i) const {return stubsbinned_[bin][i];}


  
  void clean() {
    stubs_.clear();
    for (unsigned int i=0;i<NLONGVMBINS;i++){
      stubsbinned_[i].clear();
    }
  }

  void writeStubs(bool first) {

    std::string fname="MemPrints/VMStubsTE/VMStubs_";
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

    if (layer_%2==1) {
      for (unsigned int j=0;j<stubs_.size();j++){
	string stub=stubs_[j].first->stubindex().str();
	stub+="|";
    stub+=stubs_[j].first->stubaddressaste().str();
    stub+="|";
	FPGAWord tmp;
	tmp.set(stubs_[j].first->getVMBits().value(),10,true,__LINE__,__FILE__);
	stub+=tmp.str();  
	if (j<16) out_ <<"0";
	out_ << hex << j << dec ;
	out_ <<" "<<stub<< endl;
      }
    } else if (layer_!=0) {
      for (unsigned int i=0;i<NLONGVMBINS;i++) {
	for (unsigned int j=0;j<stubsbinned_[i].size();j++){
	  string stub=stubsbinned_[i][j].first->stubindex().str();
      stub+="|";
      stub+=stubsbinned_[i][j].first->stubaddressaste().str();
	  out_ << hex << i << " " << j << dec << " "<<stub<<endl;
	}
      }
    }

    //Need to be done smarter to handle overlap region
    if (!overlap_) {
      if (disk_%2==1) {
	for (unsigned int j=0;j<stubs_.size();j++){
	  string stub=stubs_[j].first->stubindex().str();
	  stub+="|";
      stub+=stubs_[j].first->stubaddressaste().str();
      stub+="|";
	  FPGAWord tmp;
	  tmp.set(stubs_[j].first->getVMBits().value(),10,true,__LINE__,__FILE__);
	  stub+=tmp.str();  
	  if (j<16) out_ <<"0";
	  out_ << hex << j << dec ;
	  out_ <<" "<<stub<< endl;
	}
      } else if (disk_!=0) {
	for (unsigned int i=0;i<NLONGVMBINS;i++) {
	  for (unsigned int j=0;j<stubsbinned_[i].size();j++){
	    string stub=stubsbinned_[i][j].first->stubindex().str();
        stub+="|";
        stub+=stubsbinned_[i][j].first->stubaddressaste().str();
	    out_ << hex << i << " " << j << dec << " "<<stub<<endl;
	  }
	}
      }
    }
    

    out_.close();

    bx_++;
    event_++;
    if (bx_>7) bx_=0;

  }

  int phibin() const {
    return phibin_;
  }

  void getPhiRange(double &phimin, double &phimax) {

    if (layer_==1 || layer_==3 || layer_==5) {
      int nphibin=24;
      double dphi=dphisector/nphibin;
      phimax=phibin()*dphi;
      phimin=phimax-dphi;
      return;
    }

    if (layer_==2 || layer_==4 || layer_==6) {
      int nphibin=12;
      double dphi=dphisector/nphibin;
      phimax=phibin()*dphi-dphisector/6.0;
      phimin=phimax-dphi;
      return;
    }

    assert(0); //not implemented yet
    
  }
  
  void setother(FPGAVMStubsTE* other){
    other_=other;
  }
  
  FPGAVMStubsTE* other() const {
    return other_;
  }

  

private:

  int layer_;
  int disk_;
  int phibin_;
  FPGAVMStubsTE* other_;
  bool overlap_;
  double phimin_;
  double phimax_;
  std::vector<std::pair<FPGAStub*,L1TStub*> > stubs_;
  std::vector<std::pair<FPGAStub*,L1TStub*> > stubsbinned_[NLONGVMBINS];

};

#endif
