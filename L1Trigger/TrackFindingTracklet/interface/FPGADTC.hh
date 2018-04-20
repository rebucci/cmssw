// This class holds a list of stubs that are in a given layer and DCT region
#ifndef FPGADTC_H
#define FPGADTC_H

#include "L1TStub.hh"
#include "FPGAStub.hh"
#include "FPGADTCLink.hh"

using namespace std;

class FPGADTC{

public:

  FPGADTC(string name=""){
    name_=name;
  }

  void init(string name){
    name_=name;
  }
  
  void addLink(double phimin,double phimax){
    //cout << "FPGADTC addLink "<<name_<<endl;
    FPGADTCLink link(phimin,phimax);
    links_.push_back(link);
  }
  
  int addStub(std::pair<FPGAStub*,L1TStub*> stub) {
    double phi=stub.second->phi()+0.00*two_pi/NSector; //HAck to correct for offset.
    if (phi>0.5*two_pi/NSector) phi-=two_pi;
    bool overlaplayer=((stub.second->layer()+1)%2==0);
    //cout << "layer overlaplayer : "<<stub.second->layer()+1<<" "<<overlaplayer
    //	 <<endl;
    int added=0;
    //cout << "In FPGADTC : "<<name_<<" #links "<<links_.size()<<endl; 
    for (unsigned int i=0;i<links_.size();i++){
      if (links_[i].inRange(phi,overlaplayer)){
	added++;
	//cout << "Added stub in FPGADTC"<<endl;
	links_[i].addStub(stub);
      }
    }
    return added;
  }

  unsigned int nLinks() const {return links_.size();}

  const FPGADTCLink& link(unsigned int i) const {return links_[i];}

  void clean() {
    for (unsigned int i=0;i<links_.size();i++) {
      links_[i].clean();
    }
  }


private:

  string name_;
  std::vector<FPGADTCLink > links_;



};

#endif
