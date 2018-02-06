// This class holds a list of stubs for an input link.
//This modules 'owns' the pointers to the stubs. All
//subsequent modules that handles stubs uses a pointer
//to the original stored here

 
#ifndef FPGAINPUTLINK_H
#define FPGAINPUTLINK_H

#include "L1TStub.hh"
#include "FPGAStub.hh"
#include "FPGAMemoryBase.hh"
#include "FPGAVMRouterPhiCorrTable.hh"
#include <math.h>
#include <sstream>

using namespace std;

class FPGAInputLink:public FPGAMemoryBase{

public:

  FPGAInputLink(string name, unsigned int iSector, 
		double phimin, double phimax):
    FPGAMemoryBase(name,iSector){
    phimin_=phimin;
    phimax_=phimax;
  }

  void addStub(L1TStub& al1stub, FPGAStub& stub) {

    static bool first=true;
    static FPGAVMRouterPhiCorrTable phiCorrLayers[6];

    if (first) {
      for (int l=0;l<6;l++){
        phiCorrLayers[l].init(l+1,5,3);
      }
      first=false;
    }
    

    
    if (stub.layer().value()!=-1) {

      FPGAWord r=stub.r();
      FPGAWord bend=stub.bend();
      int bendbin=bend.value()>>(bend.nbits()-5);
      int rbin=(r.value()+(1<<(r.nbits()-1)))>>(r.nbits()-3);
      
      int iphicorr=phiCorrLayers[stub.layer().value()].getphiCorrValue(bendbin,rbin);

      stub.setPhiCorr(iphicorr);

    }
    
    bool add=false;
    int iphivmRaw=stub.iphivmRaw();
      
    if (stub.layer().value()!=-1) {
      //cout << "FPGAInputLink::addStub "<<stub.layer().value()<<endl;

      
      string subname=getName().substr(5,7);
      string subnamelayer=getName().substr(3,2);
      
      int layer=stub.layer().value()+1;
      
      if (layer==2&&subnamelayer=="L2"&&(subname=="PHIW_ZP"||subname=="PHIQ_ZP")) {
	if (al1stub.z()>25.0) {
	  if (iphivmRaw>=4 && iphivmRaw<=15 && subname=="PHIW_ZP") add=true;     //overlap
	  if (iphivmRaw>=16 && iphivmRaw<=27 && subname=="PHIQ_ZP") add=true;    //overlap
	}
      }
      if (layer==2&&subnamelayer=="L2"&&(subname=="PHIW_ZM"||subname=="PHIQ_ZM")) {
	if (al1stub.z()<-25.0) {
	  if (iphivmRaw>=4 && iphivmRaw<=15 && subname=="PHIW_ZM") add=true;     //overlap
	  if (iphivmRaw>=16 && iphivmRaw<=27 && subname=="PHIQ_ZM") add=true;    //overlap
	}
      }

      //Special case --- L3 is grouped with D1 for overlap seeding
	
      if (layer==3&&subnamelayer=="L3"&&(subname=="PHIW_ZP"||subname=="PHIQ_ZP")) {
	if (al1stub.z()>40.0) {
	    if (iphivmRaw>=0 && iphivmRaw<=19 && subname=="PHIW_ZP") add=true;     //overlap
	    if (iphivmRaw>=12 && iphivmRaw<=31 && subname=="PHIQ_ZP") add=true;    //overlap
	  }
      }

      if (layer==3&&subnamelayer=="L3"&&(subname=="PHIW_ZM"||subname=="PHIQ_ZM")) {
	if (al1stub.z()<-40.0) {
	    if (iphivmRaw>=0 && iphivmRaw<=19 && subname=="PHIW_ZM") add=true;     //overlap
	    if (iphivmRaw>=12 && iphivmRaw<=31 && subname=="PHIQ_ZM") add=true;    //overlap
	  }
      }

	
      if (!((layer==1&&subnamelayer=="L1")||
	    (layer==2&&subnamelayer=="L2")||
	    (layer==3&&subnamelayer=="L3")||
	    (layer==4&&subnamelayer=="L4")||
	    (layer==5&&subnamelayer=="L5")||
	    (layer==6&&subnamelayer=="L6"))){
	return;
      }
      
      //cout << "Stub candidate in "<<getName()<<" "<<subnamelayer<<" "<<subname<<" "<<iphivmRaw<<" "
      //	 <<al1stub.phi()<<" "
      // 	 <<al1stub.z()<<endl;
      
      int ASPHI = stub.phiregion();  // AllStub phi region
      
      if (subnamelayer=="L1"||subnamelayer=="L3"||subnamelayer=="L5"){
	if (stub.z().value()>0) {
      if (ASPHI==1 && subname=="PHI1_ZP") add=true;
      if (ASPHI==2 && subname=="PHI2_ZP") add=true;
      if (ASPHI==3 && subname=="PHI3_ZP") add=true;
      /*
	  if (subnamelayer=="L1") {
	    if (iphivmRaw>=4 && iphivmRaw<=7 && subname=="PHI1_ZP") add=true;
	    if (iphivmRaw>=8 && iphivmRaw<=11 && subname=="PHI2_ZP") add=true;
	    if (iphivmRaw>=12 && iphivmRaw<=15 && subname=="PHI3_ZP") add=true;
	    if (iphivmRaw>=16 && iphivmRaw<=19 && subname=="PHI4_ZP") add=true;
	    if (iphivmRaw>=20 && iphivmRaw<=23 && subname=="PHI5_ZP") add=true;
	    if (iphivmRaw>=24 && iphivmRaw<=27 && subname=="PHI6_ZP") add=true;
	  } else {
	    if (iphivmRaw>=4 && iphivmRaw<=11 && subname=="PHI1_ZP") add=true;
	    if (iphivmRaw>=12 && iphivmRaw<=19 && subname=="PHI2_ZP") add=true;
	    if (iphivmRaw>=20 && iphivmRaw<=27 && subname=="PHI3_ZP") add=true;
	  }
      */
	  //these are for TE
	  if ((subnamelayer=="L1"&&fabs(al1stub.z())<87.0)||subnamelayer=="L3"||subnamelayer=="L5"){
	    if (iphivmRaw>=4 && iphivmRaw<=11 && subname=="PHIA_ZP") add=true;
	    if (iphivmRaw>=12 && iphivmRaw<=19 && subname=="PHIB_ZP") add=true;
	    if (iphivmRaw>=20 && iphivmRaw<=27 && subname=="PHIC_ZP") add=true;
	  }
	  if ((subnamelayer=="L1"&&fabs(al1stub.z())>80.0)||subnamelayer=="L3"||subnamelayer=="L5"){
	    if (iphivmRaw>=4 && iphivmRaw<=15 && subname=="PHIX_ZP") add=true;     //overlap
	    if (iphivmRaw>=16 && iphivmRaw<=27 && subname=="PHIY_ZP") add=true;    //overlap
	  }
	} else {
      if (ASPHI==1 && subname=="PHI1_ZM") add=true;
      if (ASPHI==2 && subname=="PHI2_ZM") add=true;
      if (ASPHI==3 && subname=="PHI3_ZM") add=true;
      /*
	  if (subnamelayer=="L1") {
	    if (iphivmRaw>=4 && iphivmRaw<=7 && subname=="PHI1_ZM") add=true;
	    if (iphivmRaw>=8 && iphivmRaw<=11 && subname=="PHI2_ZM") add=true;
	    if (iphivmRaw>=12 && iphivmRaw<=15 && subname=="PHI3_ZM") add=true;
	    if (iphivmRaw>=16 && iphivmRaw<=19 && subname=="PHI4_ZM") add=true;
	    if (iphivmRaw>=20 && iphivmRaw<=23 && subname=="PHI5_ZM") add=true;
	    if (iphivmRaw>=24 && iphivmRaw<=27 && subname=="PHI6_ZM") add=true;
	  } else {
	    if (iphivmRaw>=4 && iphivmRaw<=11 && subname=="PHI1_ZM") add=true;
	    if (iphivmRaw>=12 && iphivmRaw<=19 && subname=="PHI2_ZM") add=true;
	    if (iphivmRaw>=20 && iphivmRaw<=27 && subname=="PHI3_ZM") add=true;
	  }
      */
	  //these are for TE
	  if ((subnamelayer=="L1"&&fabs(al1stub.z())<87.0)||subnamelayer=="L3"||subnamelayer=="L5"){
	    if (iphivmRaw>=4 && iphivmRaw<=11 && subname=="PHIA_ZM") add=true;
	    if (iphivmRaw>=12 && iphivmRaw<=19 && subname=="PHIB_ZM") add=true;
	    if (iphivmRaw>=20 && iphivmRaw<=27 && subname=="PHIC_ZM") add=true;
	  }
	  if ((subnamelayer=="L1"&&fabs(al1stub.z())>80.0)||subnamelayer=="L3"||subnamelayer=="L5"){
	    if (iphivmRaw>=4 && iphivmRaw<=15 && subname=="PHIX_ZM") add=true;   //overlap
	    if (iphivmRaw>=16 && iphivmRaw<=27 && subname=="PHIY_ZM") add=true;	//overlap  
	  }
	}
      }

      
      if (subnamelayer=="L2"||subnamelayer=="L4"||subnamelayer=="L6"){
	if (stub.z().value()>0) {
	  //remember that these are for ME
      if (ASPHI==1 && subname=="PHI1_ZP") add=true;
      if (ASPHI==2 && subname=="PHI2_ZP") add=true;
      if (ASPHI==3 && subname=="PHI3_ZP") add=true;
      if (ASPHI==4 && subname=="PHI4_ZP") add=true;
      /*
	  if (iphivmRaw>=4 && iphivmRaw<=7 && subname=="PHI1_ZP") add=true;
	  if (iphivmRaw>=8 && iphivmRaw<=15 && subname=="PHI2_ZP") add=true;
	  if (iphivmRaw>=16 && iphivmRaw<=23 && subname=="PHI3_ZP") add=true;
	  if (iphivmRaw>=24 && iphivmRaw<=27 && subname=="PHI4_ZP") add=true;
      */
	  //these are for TE
	  if (iphivmRaw>=0 && iphivmRaw<=7 && subname=="PHIA_ZP") add=true;
	  if (iphivmRaw>=8 && iphivmRaw<=15 && subname=="PHIB_ZP") add=true;
	  if (iphivmRaw>=16 && iphivmRaw<=23 && subname=="PHIC_ZP") add=true;
	  if (iphivmRaw>=24 && iphivmRaw<=31 && subname=="PHID_ZP") add=true;
	} else {
	  //remember that these are for ME
      if (ASPHI==1 && subname=="PHI1_ZM") add=true;
      if (ASPHI==2 && subname=="PHI2_ZM") add=true;
      if (ASPHI==3 && subname=="PHI3_ZM") add=true;
      if (ASPHI==4 && subname=="PHI4_ZM") add=true;
      /*
	  if (iphivmRaw>=4 && iphivmRaw<=7 && subname=="PHI1_ZM") add=true;
	  if (iphivmRaw>=8 && iphivmRaw<=15 && subname=="PHI2_ZM") add=true;
	  if (iphivmRaw>=16 && iphivmRaw<=23 && subname=="PHI3_ZM") add=true;
	  if (iphivmRaw>=24 && iphivmRaw<=27 && subname=="PHI4_ZM") add=true;
      */
	  //these are for TE
	  if (iphivmRaw>=0 && iphivmRaw<=7 && subname=="PHIA_ZM") add=true;
	  if (iphivmRaw>=8 && iphivmRaw<=15 && subname=="PHIB_ZM") add=true;
	  if (iphivmRaw>=16 && iphivmRaw<=23 && subname=="PHIC_ZM") add=true;
	  if (iphivmRaw>=24 && iphivmRaw<=31 && subname=="PHID_ZM") add=true;
	}
      }
    }
    else if (stub.disk().value()!=0) {
      //cout << "FPGAInputLink::addStub in disk "<<stub.disk().value()<<endl;
    
      string subname=getName().substr(5,4);
      string subnamelayer=getName().substr(3,2);
      
      int disk=stub.disk().value();

      if (abs(disk)==1&&(subnamelayer=="F1"||subnamelayer=="B1")&&(subname=="PHIW"||subname=="PHIQ")) {
	if (al1stub.r()>40.0) { 
	  if (al1stub.z()>0.0&&subnamelayer=="F1") {
	    if (iphivmRaw>=0 && iphivmRaw<=19 && subname=="PHIW") add=true;     //overlap
	    if (iphivmRaw>=12 && iphivmRaw<=31 && subname=="PHIQ") add=true;    //overlap
	  }
	  if (al1stub.z()<0.0&&subnamelayer=="B1") {
	    if (iphivmRaw>=0 && iphivmRaw<=19 && subname=="PHIW") add=true;     //overlap
	    if (iphivmRaw>=12 && iphivmRaw<=31 && subname=="PHIQ") add=true;    //overlap
	  }
	}
      }
    
      
      if (!((disk==1&&subnamelayer=="F1")||
	    (disk==2&&subnamelayer=="F2")||
	    (disk==3&&subnamelayer=="F3")||
	    (disk==4&&subnamelayer=="F4")||
	    (disk==5&&subnamelayer=="F5")||
	    (disk==-1&&subnamelayer=="B1")||
	    (disk==-2&&subnamelayer=="B2")||
	    (disk==-3&&subnamelayer=="B3")||
	    (disk==-4&&subnamelayer=="B4")||
	    (disk==-5&&subnamelayer=="B5"))){
	return;
      }
      
      //cout << "Stub candidate in "<<getName()<<" "<<subnamelayer<<" "<<subname<<" "<<iphivmRaw<<" "
      //	 <<al1stub.phi()<<" "
      // 	 <<al1stub.z()<<endl;

      int ASPHI = stub.phiregion();  // AllStub phi region
      
      if (subnamelayer=="F1"||subnamelayer=="F3"||subnamelayer=="F5"||
	  subnamelayer=="B1"||subnamelayer=="B3"||subnamelayer=="B5"){
        if (ASPHI==1 && subname=="PHI1") add=true;
        if (ASPHI==2 && subname=="PHI2") add=true;
        if (ASPHI==3 && subname=="PHI3") add=true;
        /*
	if (iphivmRaw>=4 && iphivmRaw<=11 && subname=="PHI1") add=true;
	if (iphivmRaw>=12 && iphivmRaw<=19 && subname=="PHI2") add=true;
	if (iphivmRaw>=20 && iphivmRaw<=27 && subname=="PHI3") add=true;
        */
	//these are for TE
	if ((subnamelayer=="F1"||subnamelayer=="B1")&&(subname=="PHIA"||subname=="PHIB"||subname=="PHIC"||subname=="PHID")){
	  if (iphivmRaw>=4 && iphivmRaw<=17 && subname=="PHIA") add=true;
	  if (iphivmRaw>=4 && iphivmRaw<=17 && subname=="PHIB") add=true;
	  if (iphivmRaw>=14 && iphivmRaw<=27 && subname=="PHIC") add=true;
	  if (iphivmRaw>=14 && iphivmRaw<=27 && subname=="PHID") add=true;
	}else{
	  if (iphivmRaw>=4 && iphivmRaw<=15 && subname=="PHIA") add=true;
	  if (iphivmRaw>=16 && iphivmRaw<=27 && subname=="PHIB") add=true;
	  if (iphivmRaw>=0 && iphivmRaw<=19 && subname=="PHIX") add=true;
	  if (iphivmRaw>=12 && iphivmRaw<=31 && subname=="PHIY") add=true;
	}
      }

      
      if (subnamelayer=="F2"||subnamelayer=="F4"||
	  subnamelayer=="B2"||subnamelayer=="B4"){
	//remember that these are for ME
        if (ASPHI==1 && subname=="PHI1") add=true;
        if (ASPHI==2 && subname=="PHI2") add=true;
        if (ASPHI==3 && subname=="PHI3") add=true;
        /*
	if (iphivmRaw>=4 && iphivmRaw<=11 && subname=="PHI1") add=true;
	if (iphivmRaw>=12 && iphivmRaw<=19 && subname=="PHI2") add=true;
	if (iphivmRaw>=20 && iphivmRaw<=27 && subname=="PHI3") add=true;
        */
	//these are for TE
	if ((subnamelayer=="F2"||subnamelayer=="B2")&&(subname=="PHIA"||subname=="PHIB"||subname=="PHIC"||subname=="PHID")){
	  if (iphivmRaw>=0 && iphivmRaw<=13 && subname=="PHIA") add=true;
	  if (iphivmRaw>=4 && iphivmRaw<=19 && subname=="PHIB") add=true;
	  if (iphivmRaw>=12 && iphivmRaw<=27 && subname=="PHIC") add=true;
	  if (iphivmRaw>=18 && iphivmRaw<=31 && subname=="PHID") add=true;
	}else{
	  if (iphivmRaw>=0 && iphivmRaw<=19 && subname=="PHIA") add=true;
	  if (iphivmRaw>=12 && iphivmRaw<=31 && subname=="PHIB") add=true;
	}
      }
    }
  
      
    if (!add) {
      //cout << "Will not add stub" << endl;
      return;
    }
    if (debug1) {
      cout << "Will add stub in "<<getName()<<" iphiwmRaw = "<<iphivmRaw<<" phi="<<al1stub.phi()<<" z="<<al1stub.z()<<" r="<<al1stub.r()<<endl;
    }
    if (stubs_.size()<MAXSTUBSLINK) {
      L1TStub* l1stub=new L1TStub(al1stub);
      //FPGAStub* stub=new FPGAStub(*l1stub,phimin_,phimax_);
      FPGAStub* stubptr=new FPGAStub(stub);
      std::pair<FPGAStub*,L1TStub*> tmp(stubptr,l1stub);
      stubs_.push_back(tmp);
    }
  }

  unsigned int nStubs() const {return stubs_.size();}

  FPGAStub* getFPGAStub(unsigned int i) const {return stubs_[i].first;}
  L1TStub* getL1TStub(unsigned int i) const {return stubs_[i].second;}
  std::pair<FPGAStub*,L1TStub*> getStub(unsigned int i) const {return stubs_[i];}

  void writeStubs(bool first, bool w2, bool padded) {

    //Barrel
    std::string fname="MemPrints/InputStubs/InputStubs_";
    fname+=getName();
    fname+="_";
    ostringstream oss;
    oss << iSector_+1;
    if (iSector_+1<10) fname+="0";
    fname+=oss.str();
    if(w2)
      fname+="_in2.dat";
    else
      fname+=".dat";

    if (first) {
      bx_=0;
      event_=1;
      out_.open(fname.c_str());
    }
    else 
      out_.open(fname.c_str(),std::ofstream::app);

    int nlay[6];
    for (unsigned int i=0;i<6;i++) {
      nlay[i]=0;
    }

    unsigned int maxstubslink=0;
    if (TMUX==4) maxstubslink=21;
    else if (TMUX==6) maxstubslink=33;
    else if (TMUX==8) maxstubslink=45;
    else {
      cout << "ERROR! Only TMUX=4/6/8 are supported! Exiting..." << endl;
      return;
    }

    for (unsigned int i=0;i<stubs_.size();i++){
      if(stubs_[i].first->isBarrel()){
	int lay=stubs_[i].first->layer().value();
	assert(lay>=0);
	assert(lay<6);
	nlay[lay]++;
      }
    }
    unsigned long int nlay2[6];
    nlay2[0] = nlay[0];
    for(unsigned int i=1; i<6; ++i) {
      nlay2[i] = nlay2[i-1]+nlay[i];
      if (nlay2[i]>maxstubslink) {
	nlay2[i]=maxstubslink;
      }
      assert(nlay2[i]<=maxstubslink);
    }

    string marker="";

    if(!w2) {
      //Write out the header
      out_ << "11111111111"<<marker;  
      out_ << (bitset<3>)bx_; //Dummy BX
      for (unsigned int i=0;i<3;i++) {
	out_ <<marker<<(bitset<6>)nlay[i];
      }
      out_ << endl;
      out_ << "00000000000"<<marker<<(bitset<3>)bx_;  
      for (unsigned int i=3;i<6;i++) {
	out_ << marker<<(bitset<6>)nlay[i];
      }
      out_ << endl;
      
      //now write the stubs in layer order
      string leftover="";

      
      for (unsigned int i=0;i<6;i++) {
	for (unsigned int j=0;j<stubs_.size();j++){
	  if(stubs_[j].first->isBarrel()){
	    unsigned int lay=stubs_[j].first->layer().value();
	    if (lay==i) {
	      string tmp=stubs_[j].first->strbare();
	      string stub=leftover+tmp.substr(0,32-leftover.size());
	      leftover=tmp.substr(32-leftover.size(),leftover.size()+4);
	      out_ << stub << endl;
	      if (leftover.size()==32){
		out_ << leftover << endl;
		leftover="";
	      }
	    }
	  }
	}
      }
      if (leftover.size()!=0) {
	string tmp="00000000000000000000000000000000";
	out_ << leftover << tmp.substr(0,32-leftover.size())<<endl;
      }
    }
    else{ 
      //format in2 for IPBUS
      unsigned int xline = 10;
      unsigned long int xword;
      //Header
      xword = ((unsigned long int)7<<33)|(((unsigned long int) bx_)<<25)|(0x1ffffff);
      out_<< std::hex << xline<<"000";
      out_.width(4);
      out_<<(xword>>20);
      out_<< " "<< xline<<"00";
      out_.width(5);
      out_<<(xword & 0xfffff);
      out_<<" 51000003 51000007\n";
      //Stub counts
      xline = (xline+1)&15;
      xword = (nlay2[5])|(nlay2[4]<<6)|(nlay2[3]<<12)|(nlay2[2]<<18)
	|(nlay2[1]<<24)|(nlay2[0]<<30); 
      out_<< std::hex << xline<<"000";
      out_.fill('0');
      out_.width(4);
      out_<<(xword>>20);
      out_<< " "<< xline<<"00";
      out_.width(5);
      out_<<(xword & 0xfffff);
      out_<<" 51000003 51000007\n";
      //the stubs      
      int scount = 0;
      for (unsigned int i=0;i<6;i++) {
	for (unsigned int j=0;j<stubs_.size();j++){
	  if(stubs_[j].first->isBarrel()){
	    unsigned int lay=stubs_[j].first->layer().value();
	    if (lay==i) {
	      scount++;
	      if (scount>(int)maxstubslink) continue;
	      string tmp=stubs_[j].first->strbare();
	      bitset<36> btmp(tmp);
	      xword = btmp.to_ulong();
	      xline = (xline+1)&15;
	      out_<< std::hex << xline<<"000";
	      out_.width(4);
	      out_<<(xword>>20);
	      out_<< " "<< xline<<"00";
	      out_.width(5);
	      out_<<(xword & 0xfffff);
	      out_<<" 51000003 51000007\n";
	    }
	  }
	}
      }
      if (padded) {
    for (int k=0;k<(int)maxstubslink-scount;k++) {
	  xline = (xline+1)&15;
	  out_<< std::hex << xline<<"0000000 "<< std::hex << xline<<"0000001";
	  out_<<" 51000003 51000007\n";
	}
      }
      //Trailer
      xline = (xline+1)&15;
      xword = ((unsigned long int)7<<33)|(((unsigned long int) bx_)<<25);
      out_<< std::hex << xline<<"000";
      out_.width(4);
      out_<<(xword>>20);
      out_<< " "<< xline<<"00";
      out_.width(5);
      out_<<(xword & 0xfffff);
      out_<<" 51000003 51000007\n";
    }    
  
    out_.close();
      
    //Now the disks

    //forward

    fname="MemPrints/InputStubs/InputStubsDiskF_";
    fname+=getName();
    fname+="_";
    ostringstream ossF;
    ossF << iSector_+1;
    if (iSector_+1<10) fname+="0";
    fname+=ossF.str();
    if(w2)
      fname+="_in2.dat";
    else
      fname+=".dat";

    if (first) {
      bx_=0;
      event_=1;
      out_.open(fname.c_str());
    }
    else 
      out_.open(fname.c_str(),std::ofstream::app);

    for (unsigned int i=0;i<5;i++) {
      nlay[i]=0;
      nlay2[i]=0;
    }

    for (unsigned int i=0;i<stubs_.size();i++){
      if(stubs_[i].first->isDisk()&&stubs_[i].first->disk().value()>0){
	int lay=stubs_[i].first->disk().value()-1;
	assert(lay>=0);
	assert(lay<5);
	nlay[lay]++;
      }
    }
    nlay[5] = 0;
    nlay2[5] = 0;
    nlay2[0] = nlay[0];
    for(unsigned int i=1; i<5; ++i)
      nlay2[i] = nlay2[i-1]+nlay[i];

    if(!w2) {
      //Write out the header
      out_ << "11111111111"<<marker;  
      out_ << (bitset<3>)bx_; //Dummy BX
      for (unsigned int i=0;i<3;i++) {
	out_ <<marker<<(bitset<6>)nlay[i];
      }
      out_ << endl;
      out_ << "00000000000"<<marker<<(bitset<3>)bx_;  
      for (unsigned int i=3;i<6;i++) {
	out_ << marker<<(bitset<6>)nlay[i];
      }
      out_ << endl;
      
      //now write the stubs in disk order
      string leftover="";
      
      for (unsigned int i=0;i<5;i++) {
	for (unsigned int j=0;j<stubs_.size();j++){
	  if(stubs_[j].first->isDisk()&&stubs_[j].first->disk().value()>0){
	    unsigned int lay=stubs_[j].first->disk().value()-1;
	    if (lay==i) {
	      string tmp=stubs_[j].first->strbare();
	      string stub=leftover+tmp.substr(0,32-leftover.size());
	      leftover=tmp.substr(32-leftover.size(),leftover.size()+4);
	      out_ << stub << endl;
	      if (leftover.size()==32){
		out_ << leftover << endl;
		leftover="";
	      }
	    }
	  }
	}
      }
      if (leftover.size()!=0) {
	string tmp="00000000000000000000000000000000";
	out_ << leftover << tmp.substr(0,32-leftover.size())<<endl;
      }
    }
    else{ 
      //format in2 for IPBUS
      unsigned int xline = 10;
      unsigned long int xword;
      //Header
      xword = ((unsigned long int)7<<33)|(((unsigned long int) bx_)<<25)|(0x1ffffff);
      out_<< std::hex << xline<<"000";
      out_.width(4);
      out_<<(xword>>20);
      out_<< " "<< xline<<"00";
      out_.width(5);
      out_<<(xword & 0xfffff);
      out_<<" 51000003 51000007\n";
      //Stub counts
      xline = (xline+1)&15;
      xword = (nlay2[5])|(nlay2[4]<<6)|(nlay2[3]<<12)|(nlay2[2]<<18)
	|(nlay2[1]<<24)|(nlay2[0]<<30); 
      out_<< std::hex << xline<<"000";
      out_.fill('0');
      out_.width(4);
      out_<<(xword>>20);
      out_<< " "<< xline<<"00";
      out_.width(5);
      out_<<(xword & 0xfffff);
      out_<<" 51000003 51000007\n";
      //the stubs      
      for (unsigned int i=0;i<5;i++) {
	for (unsigned int j=0;j<stubs_.size();j++){
	  if(stubs_[j].first->isDisk()&&stubs_[j].first->disk().value()>0){
	    unsigned int lay=stubs_[j].first->disk().value()-1;
	    if (lay==i) {
	      string tmp=stubs_[j].first->strbare();
	      bitset<36> btmp(tmp);
	      xword = btmp.to_ulong();
	      xline = (xline+1)&15;
	      out_<< std::hex << xline<<"000";
	      out_.width(4);
	      out_<<(xword>>20);
	      out_<< " "<< xline<<"00";
	      out_.width(5);
	      out_<<(xword & 0xfffff);
	      out_<<" 51000003 51000007\n";
	    }
	  }
	}
      }

      //Trailer
      xline = (xline+1)&15;
      xword = ((unsigned long int)7<<33)|(((unsigned long int) bx_)<<25);
      out_<< std::hex << xline<<"000";
      out_.width(4);
      out_<<(xword>>20);
      out_<< " "<< xline<<"00";
      out_.width(5);
      out_<<(xword & 0xfffff);
      out_<<" 51000003 51000007\n";

    }    
  
    out_.close();

    //back

    fname="MemPrints/InputStubs/InputStubsDiskB_";
    fname+=getName();
    fname+="_";
    ostringstream ossB;
    ossB << iSector_+1;
    if (iSector_+1<10) fname+="0";
    fname+=ossB.str();
    if(w2)
      fname+="_in2.dat";
    else
      fname+=".dat";

    if (first) {
      bx_=0;
      event_=1;
      out_.open(fname.c_str());
    }
    else 
      out_.open(fname.c_str(),std::ofstream::app);

    for (unsigned int i=0;i<5;i++) {
      nlay[i]=0;
      nlay2[i]=0;
    }

    for (unsigned int i=0;i<stubs_.size();i++){
      if(stubs_[i].first->isDisk()&&stubs_[i].first->disk().value()<0){
	int lay=-stubs_[i].first->disk().value()-1;
	assert(lay>=0);
	assert(lay<5);
	nlay[lay]++;
      }
    }
    nlay[5] = 0;
    nlay2[5] = 0;
    nlay2[0] = nlay[0];
    for(unsigned int i=1; i<5; ++i)
      nlay2[i] = nlay2[i-1]+nlay[i];

    if(!w2) {
      //Write out the header
      out_ << "11111111111"<<marker;  
      out_ << (bitset<3>)bx_; //Dummy BX
      for (unsigned int i=0;i<3;i++) {
	out_ <<marker<<(bitset<6>)nlay[i];
      }
      out_ << endl;
      out_ << "00000000000"<<marker<<(bitset<3>)bx_;  
      for (unsigned int i=3;i<6;i++) {
	out_ << marker<<(bitset<6>)nlay[i];
      }
      out_ << endl;
      
      //now write the stubs in disk order
      string leftover="";
      
      for (unsigned int i=0;i<5;i++) {
	for (unsigned int j=0;j<stubs_.size();j++){
	  if(stubs_[j].first->isDisk()&&stubs_[j].first->disk().value()<0){
	    unsigned int lay=-stubs_[j].first->disk().value()-1;
	    if (lay==i) {
	      string tmp=stubs_[j].first->strbare();
	      string stub=leftover+tmp.substr(0,32-leftover.size());
	      leftover=tmp.substr(32-leftover.size(),leftover.size()+4);
	      out_ << stub << endl;
	      if (leftover.size()==32){
		out_ << leftover << endl;
		leftover="";
	      }
	    }
	  }
	}
      }
      if (leftover.size()!=0) {
	string tmp="00000000000000000000000000000000";
	out_ << leftover << tmp.substr(0,32-leftover.size())<<endl;
      }
    }
    else{ 
      //format in2 for IPBUS
      unsigned int xline = 10;
      unsigned long int xword;
      //Header
      xword = ((unsigned long int)7<<33)|(((unsigned long int) bx_)<<25)|(0x1ffffff);
      out_<< std::hex << xline<<"000";
      out_.width(4);
      out_<<(xword>>20);
      out_<< " "<< xline<<"00";
      out_.width(5);
      out_<<(xword & 0xfffff);
      out_<<" 51000003 51000007\n";
      //Stub counts
      xline = (xline+1)&15;
      xword = (nlay2[5])|(nlay2[4]<<6)|(nlay2[3]<<12)|(nlay2[2]<<18)
	|(nlay2[1]<<24)|(nlay2[0]<<30); 
      out_<< std::hex << xline<<"000";
      out_.fill('0');
      out_.width(4);
      out_<<(xword>>20);
      out_<< " "<< xline<<"00";
      out_.width(5);
      out_<<(xword & 0xfffff);
      out_<<" 51000003 51000007\n";
      //the stubs      
      for (unsigned int i=0;i<5;i++) {
	for (unsigned int j=0;j<stubs_.size();j++){
	  if(stubs_[j].first->isDisk()&&stubs_[j].first->disk().value()<0){
	    unsigned int lay=-stubs_[j].first->disk().value()-1;
	    if (lay==i) {
	      string tmp=stubs_[j].first->strbare();
	      bitset<36> btmp(tmp);
	      xword = btmp.to_ulong();
	      xline = (xline+1)&15;
	      out_<< std::hex << xline<<"000";
	      out_.width(4);
	      out_<<(xword>>20);
	      out_<< " "<< xline<<"00";
	      out_.width(5);
	      out_<<(xword & 0xfffff);
	      out_<<" 51000003 51000007\n";
	    }
	  }
	}
      }

      //Trailer
      xline = (xline+1)&15;
      xword = ((unsigned long int)7<<33)|(((unsigned long int) bx_)<<25);
      out_<< std::hex << xline<<"000";
      out_.width(4);
      out_<<(xword>>20);
      out_<< " "<< xline<<"00";
      out_.width(5);
      out_<<(xword & 0xfffff);
      out_<<" 51000003 51000007\n";

    }    
  
    out_.close();

    //increment BX
    bx_++;
    if (bx_>7) bx_=0;

  }
  
  

  void clean() {
    for(unsigned int i=0;i<stubs_.size();i++){
      delete stubs_[i].first;
      delete stubs_[i].second;
    }
    stubs_.clear();
  }


private:

  double phimin_;
  double phimax_;
  std::vector<std::pair<FPGAStub*,L1TStub*> > stubs_;

};

#endif
