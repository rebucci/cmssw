//Base class for processing modules
#ifndef FPGAMLAYERPROJECTION_H
#define FPGAMLAYERPROJECTION_H

using namespace std;

class FPGALayerProjection{

public:

  FPGALayerProjection(){
    valid_=false;
  }

  void init(int projlayer,
	    double rproj,
	    int iphiproj,
	    int izproj,
	    int iphider,
	    int izder,
	    bool minusNeighbor,
	    bool plusNeighbor,
	    double phiproj,
	    double zproj,
	    double phiprojder,
	    double zprojder,
	    double phiprojapprox,
	    double zprojapprox,
	    double phiprojderapprox,
	    double zprojderapprox) {

    assert(projlayer>=1);
    assert(projlayer<=6);

    //cout << "Initiating projection to layer = "<<projlayer<< " at radius = "<<rproj<<endl;
    
    valid_=true;

    rproj_=rproj;
    
    projlayer_=projlayer;

    if (iphiproj<0) iphiproj=0; //FIXME should be assert?

    if (rproj<60.0) {
      fpgaphiproj_.set(iphiproj,nbitsphiprojL123,true,__LINE__,__FILE__);
      int iphivm=(iphiproj>>(nbitsphiprojL123-5))&0x7;
      if ((projlayer_%2)==1) {
	iphivm^=4;
      }
      fpgaphiprojvm_.set(iphivm,3,true,__LINE__,__FILE__);
      fpgazproj_.set(izproj,nbitszprojL123,false,__LINE__,__FILE__);
      //if (fpgazproj_.atExtreme()) {
      //	cout << "FPGALayerProjection fpgazproj at extreme (case1)" << endl;
      //}
      int izvm=izproj>>(12-7)&0xf; 
      fpgazprojvm_.set(izvm,4,true,__LINE__,__FILE__);
      fpgaphiprojder_.set(iphider,nbitsphiprojderL123,false,__LINE__,__FILE__);
      fpgazprojder_.set(izder,nbitszprojderL123,false,__LINE__,__FILE__);
    } else {
      fpgaphiproj_.set(iphiproj,nbitsphiprojL456,true,__LINE__,__FILE__);
      int iphivm=(iphiproj>>(nbitsphiprojL456-5))&0x7;
      if ((projlayer_%2)==1) {
	iphivm^=4;
      }
      fpgaphiprojvm_.set(iphivm,3,true,__LINE__,__FILE__);
      fpgazproj_.set(izproj,nbitszprojL456,false,__LINE__,__FILE__);
      int izvm=izproj>>(8-7)&0xf;
      fpgazprojvm_.set(izvm,4,true,__LINE__,__FILE__);
      fpgaphiprojder_.set(iphider,nbitsphiprojderL456,false,__LINE__,__FILE__);
      fpgazprojder_.set(izder,nbitszprojderL456,false,__LINE__,__FILE__); 
    }

    ////Separate the vm projections into zbins
    ////This determines the central bin:
    ////int zbin=3-(zproj.value()>>(zproj.nbits()-3));
    ////But we need some range:
    unsigned int zbin1=3-(((fpgazproj_.value()>>(fpgazproj_.nbits()-5))+1)>>2);
    unsigned int zbin2=3-(((fpgazproj_.value()>>(fpgazproj_.nbits()-5))-1)>>2);
    if (zbin1>7) zbin1=0; //note that zbin1 is unsigned
    if (zbin2>7) zbin2=7;
    assert(zbin1<=zbin2);
    assert(zbin2-zbin1<=1); 
    fpgazbin1projvm_.set(zbin1,3,true,__LINE__,__FILE__); // first z bin
    if (zbin1==zbin2) fpgazbin2projvm_.set(0,1,true,__LINE__,__FILE__); // don't need to check adjacent z bin
    else              fpgazbin2projvm_.set(1,1,true,__LINE__,__FILE__); // do need to check next z bin

    
    minusNeighbor_=minusNeighbor;
    plusNeighbor_=plusNeighbor;

    phiproj_=phiproj;
    zproj_=zproj;
    phiprojder_=phiprojder;
    zprojder_=zprojder;

    phiprojapprox_=phiprojapprox;
    zprojapprox_=zprojapprox;
    phiprojderapprox_=phiprojderapprox;
    zprojderapprox_=zprojderapprox;

  }

  virtual ~FPGALayerProjection(){}

  bool valid() const {
    return valid_;
  }
  
  int projlayer() const {
    assert(valid_);
    return projlayer_;
  };

  double rproj() const {
    assert(valid_);
    return rproj_;
  };

  
  FPGAWord fpgaphiproj() const {
    assert(valid_);
    return fpgaphiproj_;
  };
  
  FPGAWord fpgazproj() const {
    assert(valid_);
    return fpgazproj_;
  };

  FPGAWord fpgaphiprojder() const {
    assert(valid_);
    return fpgaphiprojder_;
  };

  FPGAWord fpgazprojder() const {
    assert(valid_);
    return fpgazprojder_;
  };

  bool minusNeighbor() const {
    assert(valid_);
    return minusNeighbor_;
  };
  
  bool plusNeighbor() const {
    assert(valid_);
    return plusNeighbor_;
  };

  FPGAWord fpgaphiprojvm() const {
    assert(valid_);
    return fpgaphiprojvm_;
  };

  FPGAWord fpgazbin1projvm() const {
    assert(valid_);
    return fpgazbin1projvm_;
  };
  
  FPGAWord fpgazbin2projvm() const {
    assert(valid_);
    return fpgazbin2projvm_;
  };

  FPGAWord fpgazprojvm() const {
    assert(valid_);
    return fpgazprojvm_;
  };

  double phiproj() const {
    assert(valid_);
    return phiproj_;
  };

  double zproj() const {
    assert(valid_);
    return zproj_;
  };

  double phiprojder() const {
    assert(valid_);
    return phiprojder_;
  };

  double zprojder() const {
    assert(valid_);
    return zprojder_;
  };

  double phiprojapprox() const {
    assert(valid_);
    return phiprojapprox_;
  };

  double zprojapprox() const {
    assert(valid_);
    return zprojapprox_;
  };

  double phiprojderapprox() const {
    assert(valid_);
    return phiprojderapprox_;
  };

  double zprojderapprox() const {
    assert(valid_);
    return zprojderapprox_;
  };

  

protected:

  bool valid_;

  int projlayer_;

  double rproj_;
  
  FPGAWord fpgaphiproj_;
  FPGAWord fpgazproj_;
  FPGAWord fpgaphiprojder_;
  FPGAWord fpgazprojder_;

  bool minusNeighbor_;
  bool plusNeighbor_;


  FPGAWord fpgaphiprojvm_;
  FPGAWord fpgazprojvm_;
  
  FPGAWord fpgazbin1projvm_;
  FPGAWord fpgazbin2projvm_;
 
  double phiproj_;
  double zproj_;
  double phiprojder_;
  double zprojder_;

  double zbin1_;
  double zbin2_;

  double phiprojapprox_;
  double zprojapprox_;
  double phiprojderapprox_;
  double zprojderapprox_;
    
};

#endif
