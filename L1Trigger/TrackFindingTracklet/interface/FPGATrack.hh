#ifndef FPGATRACK_HH
#define FPGATRACK_HH

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <vector>
#include <map>

using namespace std;

class FPGATrack{

public:

  FPGATrack(int irinv, int iphi0, int it, int iz0, int ichisq,
            double chisq,
            std::map<int, int> stubID, std::vector<L1TStub*> l1stub,
            int seed){

    irinv_=irinv;
    iphi0_=iphi0;
    iz0_=iz0;
    it_=it;
    ichisq_=ichisq;

    chisq_=chisq;

    stubID_=stubID;
    l1stub_=l1stub;

    seed_=seed;
    duplicate_=false;
    sector_=NSector;

  }


  ~FPGATrack() {

  }
  
  void setDuplicate(bool flag) { duplicate_=flag; }
  void setSector(int nsec) { sector_=nsec; }

  int irinv() const { return irinv_; }
  int iphi0() const { return iphi0_; }
  int iz0()   const { return iz0_; }
  int it()    const { return it_; }
  int ichisq() const {return ichisq_;}

  std::map<int, int> stubID() const { return stubID_; }
  std::vector<L1TStub*> stubs() const { return l1stub_; }

  int seed() const { return seed_; }
  int duplicate() const { return duplicate_; }
  int sector() const { return sector_; }

  double pt(double bfield=3.811202) const {
    if(abs(seed_) < 10) return (0.3*bfield/100.0)/(irinv_*krinvpars);
    if(abs(seed_) > 10 && abs(seed_) < 20) return (0.3*bfield/100.0)/(irinv_*krinvparsdisk);
    if(abs(seed_) > 20) return (0.3*bfield/100.0)/(irinv_*krinvparsdisk);
    return 0; }
  double phi0() const {
    if(abs(seed_) < 10) return iphi0_*kphi0pars+( (two_pi/NSector)*sector_ - (two_pi/NSector)/6. );
    if(abs(seed_) > 10 && abs(seed_) < 20) return iphi0_*kphi0parsdisk+( (two_pi/NSector)*sector_ - (two_pi/NSector)/6. );
    if(abs(seed_) > 20) return iphi0_*kphi0parsdisk+( (two_pi/NSector)*sector_ - (two_pi/NSector)/6. );
    return 0; }
  double eta() const {
    if(abs(seed_) < 10) return -log(tan(0.5*(0.25*two_pi-atan(it_*ktpars))));
    if(abs(seed_) > 10 && abs(seed_) < 20) return -log(tan(0.5*(0.25*two_pi-atan(it_*ktparsdisk))));
    if(abs(seed_) > 20) return -log(tan(0.5*(0.25*two_pi-atan(it_*ktparsdisk))));
    return 0; }
  double z0() const {
    if(abs(seed_) < 10) return iz0_*kz;
    if(abs(seed_) > 10 && abs(seed_) < 20) return iz0_*kzdisk;
    if(abs(seed_) > 20) return iz0_*kzdisk;
    return 0; }
  double rinv() const {
    if(abs(seed_) < 10) return irinv_*krinvpars;
    if(abs(seed_) > 10 && abs(seed_) < 20) return irinv_*krinvparsdisk;
    if(abs(seed_) > 20) return irinv_*krinvparsdisk;
    return 0; }
  double d0() const {return 0.0;} //Fix when fit for 5 pars
  double chisq() const {return chisq_;}



private:
  
  int irinv_;
  int iphi0_;
  int iz0_;
  int it_;
  int ichisq_;

  double rinv_;
  double phi0_;
  double z0_;
  double t_;
  double chisq_;

  std::map<int, int> stubID_;
  std::vector<L1TStub*> l1stub_;

  int seed_;
  bool duplicate_;
  int sector_;

};

#endif



