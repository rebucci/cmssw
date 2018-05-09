#ifndef FPGASTUB_H
#define FPGASTUB_H

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <vector>
#include "L1TStub.hh"

#include "FPGAWord.hh"
#include "FPGAConstants.hh"

using namespace std;

class FPGAStub{

public:

  FPGAStub() {
  
  }
  

  FPGAStub(const L1TStub& stub,double phiminsec, double phimaxsec) {

    //cout << "FPGASTub : making stub"<<endl;
    
    double r=stub.r();
    double z=stub.z();
    double ptinv=1.0/stub.pt();
    double sbend = stub.bend();

    int ibend=sbend*2+15;
    if (ibend<0) ibend=0;
    if (ibend>31) ibend=31;

    bend_.set(ibend,5,true,__LINE__,__FILE__);
    
    isPSmodule_ = false;
    if (stub.isPSmodule()) isPSmodule_=true;
    
    //HACK!!! seems like stubs in negative disk has wrong sign!
    if (z<-120.0) ptinv=-ptinv;
    //cout << "z stub.pt() : "<<z<<" "<<stub.pt()<<endl;
    int ipt=0;
    int layer = stub.layer()+1; 

    if(!enstubbend){
      if (fabs(ptinv)<0.4) ipt=0;
      if (fabs(ptinv)<0.3) ipt=1;
      if (fabs(ptinv)<0.2) ipt=2;
      if (fabs(ptinv)<0.1) ipt=3;
      if (ptinv<0.0) ipt=7-ipt;
    }
    //Stub pt encoding based on bend
    if(enstubbend){
      if(layer == 1 || layer == 2){                      
	if(fabs(sbend)  >=  2   ) ipt = 7;               
	if(fabs(sbend)  ==  1.5 ) ipt = 6;               
	if(fabs(sbend)  ==  1   ) ipt = 5;               
	if(fabs(sbend)  <=  0.5 ) ipt = 4;               
	if( sbend < -0.5 ) ipt = 8 - ipt;                                                       
      }                                                  
      if(layer == 3){                                    
	if(fabs(sbend)  >=  2.5 || fabs(sbend)  ==  3  )	ipt = 7;                                              
	if(fabs(sbend)  ==  2 || fabs(sbend)  ==  1.5) ipt = 6;                                            
	if(fabs(sbend)  ==  1 ) ipt = 5;                 
	if(fabs(sbend)  <=  0.5 ) ipt = 4;               
	if( sbend < -0.5 ) ipt = 8 - ipt;                
      }                                                  
      
      if(layer == 4){                                    
	if(fabs(sbend)  >=  4   || fabs(sbend) == 3.5 || fabs(sbend) == 3 ) ipt = 7;                           
	if(fabs(sbend)  ==  2.5 || fabs(sbend) == 2 ) ipt = 6;                                              
	if(fabs(sbend)  ==  1.5 ) ipt = 5;               
	if(fabs(sbend)  <=  1   ) ipt = 4;               
	if( sbend < -1 ) ipt = 8 - ipt;                  
      }                                                  
      
      
      if(layer == 5){                                    
	if(fabs(sbend)  >=  5.5 || fabs(sbend)  ==  5  || fabs(sbend)  ==  4.5 ) ipt = 7;                     
	if(fabs(sbend)  ==  4   || fabs(sbend)  ==  3.5|| fabs(sbend)  ==  3   ) ipt = 6;                     
	if(fabs(sbend)  ==  2.5 || fabs(sbend)  ==  2  || fabs(sbend)  ==  1.5 ) ipt = 5;                     
	if(fabs(sbend)  <=  1 ) ipt = 4;                 
	if( sbend < -1 ) ipt = 8 - ipt;                  
      }                                  
      
      if(layer == 6){                                                                          
	if(fabs(sbend)  >=  6.5 || fabs(sbend)  == 6   || fabs(sbend)  == 5.5 ) ipt= 7;        
	if(fabs(sbend)  ==  5   || fabs(sbend)  == 4.5 || fabs(sbend)  == 4   ) ipt= 6;        
	if(fabs(sbend)  ==  3.5 || fabs(sbend)  == 3   || fabs(sbend)  == 2.5 || fabs(sbend) == 2 ) ipt = 5;           
	if(fabs(sbend)  <=  1.5 ) ipt = 4;                                                     
	if( sbend < -1.5 ) ipt = 8 - ipt;                                                      
      }        
    }

    // hold the real values from L1Stub	
    stubphi_=stub.phi();
    stubr_  =stub.r();
    stubz_  =stub.z();
    stubrpt_=stub.pt();

    stubphimaxsec_ = phimaxsec;
    stubphiminsec_ = phiminsec;
   

    isbarrel_=false;

    if (layer<999) {

      isbarrel_=true;

      disk_.set(0,4,false,__LINE__,__FILE__);
	  
      double rmin=-1.0;
      double rmax=-1.0;

      if (layer==1) {rmin=rminL1; rmax=rmaxL1;}
      if (layer==2) {rmin=rminL2; rmax=rmaxL2;}
      if (layer==3) {rmin=rminL3; rmax=rmaxL3;}
      if (layer==4) {rmin=rminL4; rmax=rmaxL4;}
      if (layer==5) {rmin=rminL5; rmax=rmaxL5;}
      if (layer==6) {rmin=rminL6; rmax=rmaxL6;}


      assert(rmin>0.0);
      assert(rmax>0.0);
      if (r<rmin||r>rmax) cout << "Error r, rmin, rmax :"<<r
			       <<" "<<rmin<<" "<<rmax<<endl;

      int irbits=nbitsrL123;
      if (layer>=4) irbits=nbitsrL456;
      
      int ir=round_int((1<<irbits)*((r-rmean[layer-1])/(rmax-rmin)));

      //cout << "r rmean ir "<<r<<" "<<rmean[layer-1]<<" "<<ir<<endl;

      double zmin=-zlength;
      double zmax=zlength;
    
      if (z<zmin||z>zmax) cout << "Error z, zmin, zmax :"<<z
			     <<" "<<zmin<<" "<<zmax<<endl;
    
      int izbits=nbitszL123;
      if (layer>=4) izbits=nbitszL456;
      
      int iz=round_int((1<<izbits)*z/(zmax-zmin));
      
      if (z<zmin||z>zmax) cout << "Error phi, phimin, phimax :"<<stubphi_
			       <<" "<<phiminsec<<" "<<phimaxsec<<endl;
      
      assert(phimaxsec-phiminsec>0.0);
      //cout << "stubphi_ phiminsec phiminsec-(phimaxsec-phiminsec)/6.0 : "
      //	   << stubphi_<<" "<<phiminsec<<" "
      //	   <<phiminsec-(phimaxsec-phiminsec)/6.0<<endl;
      if (stubphi_<phiminsec-(phimaxsec-phiminsec)/6.0) {
	stubphi_+=two_pi;
      }
      assert((phimaxsec-phiminsec)>0.0);

      //assert(stubphi_-phimin>0.0);  //These two are not correct when
      //assert(stubphi_<phimax);      //we allow for duplications!

      int iphibits=nbitsphistubL123;
      if (layer>=4) iphibits=nbitsphistubL456;


      //cout << "phimax-phimin : "<<phimax-phimin<<" "<<two_pi/NSector<<endl;

      //cout << "phi phimin phimax : "<<stubphi_<<" "<<phiminsec
      //	   <<" "<<phimaxsec<<endl;

      double deltaphi=stubphi_-phiminsec;
      if (deltaphi>0.5*two_pi) deltaphi-=two_pi;
      
      int iphi=(1<<iphibits)*(0.125+0.75*(deltaphi)/(phimaxsec-phiminsec));
      
      //double dphi=stubphi_-phiminsec+(phimaxsec-phiminsec)/6.0-iphi*kphi;

      //cout << "FPGAStub: "<<deltaphi<<" "<<phiminsec<<" "<<phimaxsec<<" "<<kphi<<endl;

      phitmp_=stubphi_-phiminsec+(phimaxsec-phiminsec)/6.0;


      phimin_=phiminsec;


      //cout << "iphi second :"<<iphi<<" "<<(iphi&0xffffc)<<endl;

      layer_.set(layer-1,3,true,__LINE__,__FILE__);
      stubpt_.set(ipt,3,true,__LINE__,__FILE__);
      r_.set(ir,irbits,false,__LINE__,__FILE__);
      z_.set(iz,izbits,false,__LINE__,__FILE__);
      phi_.set(iphi,iphibits,true,__LINE__,__FILE__);

      /*
      double Delta=stub.r()-rmean[layer-1];
      double dphi=Delta*stub.bend()*0.009/0.18/rmean[layer-1];

      int idphi=0;
  
      if (layer<=3) {
        idphi=dphi/kphi;
      } else {
        idphi=dphi/kphi1;
      }

     
      //cout << "iphi idphi "<<phi_.value()<<" "<<idphi<<endl;
    
      int iphicorr=iphi+idphi;

      if (iphicorr<0) iphicorr=0;
      if (iphicorr>=(1<<phi_.nbits())) iphicorr=(1<<phi_.nbits())-1;

        
      */

      phicorr_.set(iphi,iphibits,true,__LINE__,__FILE__);

      
      int iphivm=0;
      
      iphivm=(iphi>>(iphibits-(Nphibits+VMphibits)))&((1<<VMphibits)-1);
      
      if (layer==1||layer==3||layer==5) {
	iphivm^=(1<<(VMphibits-1));
      }

      phivm_.set(iphivm,VMphibits,true,__LINE__,__FILE__);

    } else {
      
      // Here we handle the hits on disks.

      int disk=stub.module();
      assert(disk>0);
      if (z<0.0) disk=-disk;
      int sign=1;
      if (disk<0) sign=-1;

      double zmin=0.0;
      double zmax=0.0;
      if(enstubbend){
      //pt encoding based on position                                                        
                                                                                             
      if(r  > 245.559 && r  < 407.415 ){//rings 1-4                                          
        if(fabs(sbend)  ==  2  ) ipt = 7;                                                    
        if(fabs(sbend)  ==  1.5) ipt = 6;                                                    
        if(fabs(sbend)  ==  1  ) ipt = 5;                                                    
        if(fabs(sbend)  <=  0.5) ipt = 4;                                                    
        if( sbend < -0.5 ) ipt = 8 - ipt;                                                    
      }                                                                                      
      if(r  > 407.415 && r  < 516.901 ){//rings 5-7                                          
        if(fabs(sbend)  ==  2.5) ipt = 7;                                                    
        if(fabs(sbend)  ==  2  ) ipt = 6;                                                    
        if(fabs(sbend)  ==  1.5) ipt = 5;                                                    
        if(fabs(sbend)  <=  1  ) ipt = 4;                                                    
        if( sbend < -1  ) ipt = 8 -  ipt;                                                     
      }                                                                                      
      if(r  > 519.273 && r  < 558.049 ){//ring 8                                             
        if(fabs(sbend)  ==  3  ) ipt = 7;                                                    
        if(fabs(sbend)  ==  2.5) ipt = 6;                                                    
        if(fabs(sbend)  ==  2 || fabs(sbend)  == 1.5 ) ipt = 5;                              
        if(fabs(sbend)  <=  1  ) ipt = 4;                                                    
        if( sbend < -1  ) ipt = 8 - ipt;                                                     
      }                                                                                      
      if(r  > 558.049 && r  < 596.825 ){//ring 9                                             
        if(fabs(sbend)  ==  3.5) ipt = 7;                                                    
        if(fabs(sbend)  ==  3 || fabs(sbend) == 2.5 ) ipt = 6;                               
        if(fabs(sbend)  ==  2 || fabs(sbend) == 1.5 ) ipt = 5;                               
        if(fabs(sbend)  <=  1) ipt = 4;                                                      
        if( sbend < -1 ) ipt = 8 - ipt;                                                      
      }                                                                                      
      if(r  > 600.567 && r  < 1100.000 ){//ring 10 to 15                                     
        if(fabs(sbend)  ==  5.5 || fabs(sbend)  ==  5) ipt = 7;                              
        if(fabs(sbend)  ==  4.5 || fabs(sbend)  ==  4 || fabs(sbend)  ==  3.5 || fabs(sbend) ==  3) ipt = 6;
        if(fabs(sbend)  ==  2.5 || fabs(sbend)  ==  2 || fabs(sbend)  ==  1.5  ) ipt = 5;    
        if(fabs(sbend)  <=  1) ipt = 4;                                                      
        if( sbend < -1 ) ipt = 8 - ipt ;                                                      
      } 
      }

      if (disk==1) {zmin=zminD1; zmax=zmaxD1;}
      if (disk==2) {zmin=zminD2; zmax=zmaxD2;}
      if (disk==3) {zmin=zminD3; zmax=zmaxD3;}
      if (disk==4) {zmin=zminD4; zmax=zmaxD4;}
      if (disk==5) {zmin=zminD5; zmax=zmaxD5;}

      if (disk==-1) {zmax=-zminD1; zmin=-zmaxD1;}
      if (disk==-2) {zmax=-zminD2; zmin=-zmaxD2;}
      if (disk==-3) {zmax=-zminD3; zmin=-zmaxD3;}
      if (disk==-4) {zmax=-zminD4; zmin=-zmaxD4;}
      if (disk==-5) {zmax=-zminD5; zmin=-zmaxD5;}

      if ((z>zmax)||(z<zmin)) {
	cout << "Error disk z, zmax, zmin: "<<z<<" "<<zmax<<" "<<zmin<<endl;
      }

      int iz=(1<<nzbitsdisk)*((z-sign*zmean[abs(disk)-1])/fabs(zmax-zmin));

      assert(phimaxsec-phiminsec>0.0);
      if (stubphi_<phiminsec-(phimaxsec-phiminsec)/6.0) {
	stubphi_+=two_pi;
      }

      //Generates errors for overlap stubs
      //if (stubphi_<phiminsec||stubphi_>phimaxsec) {
      //	cout << "Error disk phi, phimin, phimax :"
      //	     <<stubphi_
      //	     <<" "<<phiminsec<<" "<<phimaxsec<<endl;
      //}
      
      assert(phimaxsec-phiminsec>0.0);
      if (stubphi_<phiminsec-(phimaxsec-phiminsec)/6.0) {
	stubphi_+=two_pi;
      }

      int iphibits=nbitsphistubL123;
      //if (layer>=4) iphibits=nbitsphistubL456; //Need to figure out this...
      //cout << "phimax-phimin : "<<phimax-phimin<<" "<<two_pi/NSector<<endl;

      double deltaphi=stubphi_-phiminsec;
      if (deltaphi>0.5*two_pi) deltaphi-=two_pi;
      
      int iphi=(1<<iphibits)*(0.125+0.75*(deltaphi)/(phimaxsec-phiminsec));

      double rmin=rmindisk;
      double rmax=rmaxdisk;
    
      if (r<rmin||r>rmax) cout << "Error disk r, rmin, rmax :"<<r
			     <<" "<<rmin<<" "<<rmax<<endl;
    
      int ir=(1<<nrbitsdisk)*(r-rmin)/(rmax-rmin);

      int irSS = -1;
      if (!isPSmodule_) {
	for (int i=0; i<20; ++i){
	  if (fabs(r-rDSS[i])<.2){
	    irSS = i;
	    break;
	  }
	}
	if (irSS<0) {
	  cout << "WARNING! didn't find rDSS value! r = " << r << endl;
	  assert(0);
	}
      }
      if(irSS < 0){
	//PS modules
      	//cout << "ir irbits : "<<ir<<" "<<irbits<<endl;
	r_.set(ir,nrbitsdisk,true,__LINE__,__FILE__);
      }
      else {
	//SS modules
	r_.set(irSS,5,true,__LINE__,__FILE__);  // in case of SS modules, store index, not r itself
      }

      //cout << "iz izbits : "<<iz<<" "<<izbits<<" "<<disk<<endl;
      z_.set(iz,nzbitsdisk,false,__LINE__,__FILE__);
      phi_.set(iphi,iphibits,true,__LINE__,__FILE__);
      stubpt_.set(ipt,3,true,__LINE__,__FILE__);

      phicorr_.set(iphi,iphibits,true,__LINE__,__FILE__);    //FIXME For now not corrected  
      
      int iphivm=0;

      iphivm=(iphi>>(iphibits-(Nphibits+VMphibits)))&((1<<VMphibits)-1);
      
      if ((abs(disk)%2)==0) {
        iphivm^=(1<<(VMphibits-1));
      }

      //iphivm=(iphi>>(iphibits-5))&0x7;
      //if ((abs(disk)%2)==1) {
      //  iphivm^=4;
      //}

      //cout << "iphivm :"<<iphivm<<endl;

      disk_.set(disk,4,false,__LINE__,__FILE__);    
      phivm_.set(iphivm,3,true,__LINE__,__FILE__);

      double alpha=stub.alpha();
      assert(fabs(alpha)<alphamax);
      int ialpha=round_int(alpha/kalpha);

      //if (stub.r()>100) cout << "r alpha ialpha : "<<stub.r()<<" "<<alpha<<" "<<ialpha<<endl;
      
      alpha_.set(ialpha,nbitsalpha,false,__LINE__,__FILE__);

      double alphanew=stub.alphanew();
      assert(fabs(alphanew)<1.0);
      int ialphanew=alphanew*(1<<(nbitsalpha-1));
      assert(ialphanew<(1<<(nbitsalpha-1)));
      assert(ialphanew>=-(1<<(nbitsalpha-1)));
      alphanew_.set(ialphanew,nbitsalpha,false,__LINE__,__FILE__);

      //cout << "alphanew ialphanew "<<alphanew<<" "<<ialphanew<<" "<<alphanew_.value()<<endl;
	    

      
    }

  }

 
  ~FPGAStub() {

  }


  //Returns a number from 0 to 31
  int iphivmRaw() const {

    int iphivm=(phicorr_.value()>>(phicorr_.nbits()-5));
    assert(iphivm>=0);
    assert(iphivm<32);
    return iphivm;
    
  }

  //VMbits is the number of bits for the fine bins. E.g. 32 bins would use VMbits=5
  //finebits is the number of bits within the VM 
  
  int iphivmFineBins(int VMbits, int finebits) const {

    return (phicorr_.value()>>(phicorr_.nbits()-VMbits-finebits))&((1<<finebits)-1);

  }


    //Returns a number from 0 to 31
  int iphivmRawPlus() const {

    //cout << layer_.value()<<" "<<disk_.value()<<endl;

    //cout << "bits : "<<phi_.value()<<" "<<(1<<(phi_.nbits()-8))<<endl;
    
    int iphivm=((phi_.value()+(1<<(phi_.nbits()-7)))>>(phi_.nbits()-5));
    if (iphivm<0) iphivm=0;
    if (iphivm>31) iphivm=0;
    return iphivm;
    
  }

  int iphivmRawMinus() const {

    //cout << layer_.value()<<" "<<disk_.value()<<endl;
    
    int iphivm=((phi_.value()-(1<<(phi_.nbits()-7)))>>(phi_.nbits()-5));
    if (iphivm<0) iphivm=0;
    if (iphivm>31) iphivm=0;
    return iphivm;
    
  }
  
  std::string str() const {
    
    std::ostringstream oss;
    //oss << stubpt_.str()<<"|"<<r_.str()<<"|"
    //    << z_.str()<<"|"<< phi_.str();
    oss << bend_.str()<<"|"<<r_.str()<<"|"
        << z_.str()<<"|"<< phi_.str();

    return oss.str();

  }
  std::string strdisk() const {
   
    std::ostringstream oss;
    //if (isPSmodule())
    //  oss << stubpt_.str()<<"||"<<r_.str()<<"|"
    //      << z_.str()<<"|"<< phi_.str();
    //else
    //  oss << stubpt_.str()<<"|"<<alpha_.str()<<"|0"<<r_.str()<<"|"
    //      << z_.str()<<"|"<< phi_.str();
    if (isPSmodule())
      oss << bend_.str()<<"|"<<r_.str()<<"|"
          << z_.str()<<"|"<< phi_.str();
    else
      oss << bend_.str()<<"|"<<alpha_.str()<<"|0"<<r_.str()<<"|"
          << z_.str()<<"|"<< phi_.str();

    return oss.str();

  }

  std::string str_phys() const {

    std::ostringstream oss;
   
    int ilz_   = 1;
    float nbsr = 128.0;
    float nphibs = 16384;
    float rmean_ = rmean[layer_.value() ];
    int Layer = layer_.value() + 1;

    double rmin=1.0;
    double rmax=10.0;

    if (Layer==1) {rmin=rminL1; rmax=rmaxL1;}
    if (Layer==2) {rmin=rminL2; rmax=rmaxL2;}
    if (Layer==3) {rmin=rminL3; rmax=rmaxL3;}
    if (Layer==4) {rmin=rminL4; rmax=rmaxL4;}
    if (Layer==5) {rmin=rminL5; rmax=rmaxL5;}
    if (Layer==6) {rmin=rminL6; rmax=rmaxL6;}
    if(Layer > 3){
      ilz_   = 16;
      nbsr   = 256.0;
      nphibs = 131072;
    }

    double pt_r = -99;
    if(stubpt_.value() == 0 ) pt_r = 1/0.4;
    if(stubpt_.value() == 1 ) pt_r = 1/0.3;
    if(stubpt_.value() == 2 ) pt_r = 1/0.2;
    if(stubpt_.value() == 3 ) pt_r = 1/0.1;
    if(stubpt_.value() == 4 ) pt_r = -1/0.1;
    if(stubpt_.value() == 5 ) pt_r = -1/0.2;
    if(stubpt_.value() == 6 ) pt_r = -1/0.3;
    if(stubpt_.value() == 7 ) pt_r = -1/0.4;

    double rreal = r_.value()/nbsr*(rmax - rmin) + rmean_;
    double phireal = ((phi_.value()*1.0)/nphibs  - 0.125 )*(1.0/0.75)*(stubphimaxsec_ - stubphiminsec_) + stubphiminsec_;




    oss << pt_r<<" "
	<<rreal<<" "
        << z_.value()*kz*ilz_<<" "
        <<phireal;


    /*  //For comparison with the real quantities                                                                 
	oss << pt_r<<"| fl :"<<stubrpt_ <<"|"                                                                         
        <<rreal<<"| fl :"<<stubr_ <<"|"                                                                           
        << z_.value()*kz*ilz_<<"|fl :"<<stubz_ <<"|"                                                              
        <<phireal<<"| fl :"<<stubphi_ ;                                                                           
    */

    return oss.str();

  }





  std::string strbare() const {
    
    std::ostringstream oss;
    //oss << stubpt_.str()<<r_.str()
	//<< z_.str()<< phi_.str();
    oss << bend_.str()<<r_.str()
	<< z_.str()<< phi_.str();

    return oss.str();

  }

  std::string strbareUNFLIPPED() const {
    
    std::ostringstream oss;
    oss << r_.str()
	<< z_.str()<< phi_.str()<<stubpt_.str();

    return oss.str();

  }

  std::string inputstr() const {
    
    std::ostringstream oss;
    oss << r_.str()<< z_.str()
	<< phi_.str()<<stubpt_.str();

    return oss.str();

  }


  std::string rawstr() const {
    
    std::ostringstream oss;
    oss << layer_.str()<<"|"<<stubpt_.str()<<"|"<< r_.str()<<"|" 
  	<< z_.str() <<"|"<< phi_.str();
    
    return oss.str();
    
  }

  std::string vmstr() const {
    
    std::ostringstream oss;
    oss << stubpt_.str() <<"|"<<stubindex_.str()<<"|"<< phivm_.str();

    return oss.str();

  }


  std::string phiregionaddressstr() {

    std::ostringstream oss;
	assert(phiregion().value()>-1);
	oss << phiregion().str() << stubindex_.str();

	return oss.str();
	
  }

  int ilink() const {

    //changed pow(2,phi_.nbits()) to (1<<phi_.nbits()), etc
    if (phi_.value()<0.33*(1<<phi_.nbits()) ) return 1;
    if (phi_.value()<0.66*(1<<phi_.nbits()) ) return 2;
    return 3;

  }
	
  FPGAWord phiregion() const {
	// 3 bits
	int iphiregion = 7;
	
	if (layer_.value()==0 or layer_.value()==2 or layer_.value()==4) { // L1, L3, L5
      if (iphivmRaw()>=4 and iphivmRaw()<=11) iphiregion = 0;
	  else if (iphivmRaw()>=12 and iphivmRaw()<=19) iphiregion = 1;
	  else if (iphivmRaw()>=20 and iphivmRaw()<=27) iphiregion = 2;
    }
    else if (layer_.value()==1 or layer_.value()==3 or layer_.value()==5) { // L2, L4, L6
      if (iphivmRaw()>=4 and iphivmRaw()<=7) iphiregion = 0;
	  else if (iphivmRaw()>=8 and iphivmRaw()<=15) iphiregion = 1;
	  else if (iphivmRaw()>=16 and iphivmRaw()<=23) iphiregion = 2;
	  else if (iphivmRaw()>=24 and iphivmRaw()<=27) iphiregion = 3;
    }
    else if (abs(disk_.value())>=1 and abs(disk_.value())<=5) { // Disk
      if (iphivmRaw()>=4 and iphivmRaw()<=11) iphiregion = 0;
	  else if (iphivmRaw()>=12 and iphivmRaw()<=19) iphiregion = 1;
	  else if (iphivmRaw()>=20 and iphivmRaw()<=27) iphiregion = 2;
    }
    //else
    //  return -1;

	FPGAWord phi;
	phi.set(iphiregion,3);

	return phi;
  }

 
  void setAllStubIndex(int nstub){
    if (nstub>=(1<<6)){
      if (debug1) cout << "Warning too large stubindex!"<<endl;
      nstub=(1<<6)-1;
    }

    stubindex_.set(nstub,6);
  }

  void setAllStubAddressTE(int nstub){
    if (nstub>=(1<<6)){
      if (debug1) cout << "Warning too large stubindex!"<<endl;
      nstub=(1<<6)-1;
    }

    stubaddressaste_.set(nstub,6);
  }


  void setPhiCorr(int phiCorr){

    /*    
    int layer=layer_.value()+1;
    
    double Delta=stubr_-rmean[layer-1];
    double dphi=Delta*0.5*(bend_.value()-15.0)*0.009/0.18/rmean[layer-1];

    int idphi=0;

    if (layer<=3) {
      idphi=dphi/kphi;
    } else {
        idphi=dphi/kphi1;
    }
    */
     
    //cout << "iphi idphi "<<phi_.value()<<" "<<idphi<<endl;
    
    //int iphicorr=phi_.value()+idphi;
    int iphicorr=phi_.value()+phiCorr;

    //cout << "phiCorr: layer bend old, new : "<<layer<<" "<<bend_.value()<<" "<<idphi<<" "<<phiCorr<<endl;
    
    if (iphicorr<0) iphicorr=0;
    if (iphicorr>=(1<<phi_.nbits())) iphicorr=(1<<phi_.nbits())-1;

    phicorr_.set(iphicorr,phi_.nbits(),true,__LINE__,__FILE__);

  }
  
  void setVMBits(int bits){
    int nbits=-1;
    assert(vmbits_.value()==-1); //Should never change the value; -1 means uninitialized
    if (layer_.value()==0 or layer_.value()==2 or layer_.value()==4) { // L1, L3, L5
      nbits=2*NLONGVMBITS+1+3;
    }
    if (layer_.value()==1 or layer_.value()==3 or layer_.value()==5) { // L2, L4, L6
      nbits=2*NLONGVMBITS+1+3;
    }
    int disk=abs(disk_.value());
    if (disk==1 or disk==3) { // D1, D3
      nbits=2*NLONGVMBITS+3;
    }
    if (disk==2 or disk==4) { // D2, D4
      nbits=2*NLONGVMBITS;
    }
    //cout << "layer, disk : "<<layer_.value()<<" "<<disk_.value()<<endl;
    assert(nbits>=0);
    vmbits_.set(bits,nbits,true,__LINE__,__FILE__);
  }

  FPGAWord getVMBits() const { return vmbits_; }
  
  FPGAWord stubpt() const { return stubpt_; }

  FPGAWord phivm() const {return phivm_; }

  FPGAWord bend() const {return bend_; }

  FPGAWord r() const { return r_; }
  FPGAWord z() const { return z_; }
  FPGAWord phi() const { return phi_; }
  FPGAWord alpha() const { return alpha_; }
  FPGAWord alphanew() const { return alphanew_; }


  int ir() const { return r_.value(); }
  int iz() const { return z_.value(); }
  int iphi() const { return phi_.value(); }

  double phitmp() const {return phitmp_;}
  double phimin() const {return phimin_;}

  FPGAWord stubindex() const {return stubindex_;}
  FPGAWord stubaddressaste() const {return stubaddressaste_;} 

  FPGAWord layer() const {return layer_;}

  FPGAWord disk() const {return disk_;}

  double stubr() const { return stubr_;}
  double stubphi() const { return stubphi_;}
  double stubz() const { return stubz_;}
  double stubrpt() const { return stubrpt_;}

  bool isBarrel() const {return isbarrel_;}
  bool isDisk() const {return !isbarrel_;}

  bool isPSmodule() const {return isPSmodule_;}

  int round_int( double r ) {
    return (r > 0.0) ? (r + 0.5) : (r - 0.5); 
  }

  double rapprox(){
    if (disk_.value()==0){
      int lr=2;
      if (layer_.value()>=3) {
	lr=1;
      }
      return r_.value()*kr*lr+rmean[layer_.value()];
    }
    return r_.value()*kr+rmindisk;
  }

  double zapprox() {
    if (disk_.value()==0){
      int lz=1;
      if (layer_.value()>=3) {
	lz=16;
      }
      return z_.value()*kz*lz;
    }
    int sign=1;
    if (disk_.value()<0) sign=-1;
    return z_.value()*kz+sign*zmean[abs(disk_.value())-1];
  }

  double phiapprox(double phimin, double phimax){
    int lphi=1;
    if (layer_.value()>=3) {
      lphi=8;
    }
    double phi=phimin-(phimax-phimin)/6.0+phi_.value()*kphi/lphi;
    if (phi>0.5*two_pi) phi-=two_pi;
    if (phi<-0.5*two_pi) phi+=two_pi;
    return phi;
  }

  
private:

  bool isbarrel_;
  bool isPSmodule_;
  FPGAWord layer_;  
  FPGAWord disk_;  
  FPGAWord stubpt_;
  FPGAWord r_;
  FPGAWord z_;
  FPGAWord phi_;
  FPGAWord alpha_;
  FPGAWord alphanew_;

  FPGAWord bend_;
  
  FPGAWord phicorr_;  //Corrected for bend to nominal radius
  
  FPGAWord phivm_;
  FPGAWord stubindex_;
  FPGAWord stubaddressaste_;

  FPGAWord vmbits_;
  
  double stubphi_;
  double stubr_;
  double stubz_;
  double stubrpt_;

  double stubphiminsec_;
  double stubphimaxsec_;


  double phitmp_;
  double phimin_;

};



#endif



