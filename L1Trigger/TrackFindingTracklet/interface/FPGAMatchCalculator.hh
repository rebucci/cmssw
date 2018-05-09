//This class implementes the match calculator
#ifndef FPGAMATCHCALCULATOR_H
#define FPGAMATCHCALCULATOR_H

#include "FPGAProcessBase.hh"

using namespace std;

class FPGAMatchCalculator:public FPGAProcessBase{

public:

  FPGAMatchCalculator(string name, unsigned int iSector):
    FPGAProcessBase(name,iSector){
    
    fullmatchesToPlus_=0;
    fullmatchesToMinus_=0;
    double dphi=two_pi/NSector;
    phimin_=iSector*dphi;
    phimax_=phimin_+dphi;
    if (phimin_>0.5*two_pi) phimin_-=two_pi;
    if (phimax_>0.5*two_pi) phimax_-=two_pi;
    phioffset_=phimin_-(dphi)/6.0;
    string subname=name.substr(8,2);
    fullmatchesToPlus_=0;
    fullmatchesToMinus_=0;
    layer_=0;
    disk_=0;
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
    if (layer_==0 && disk_==0) {
      cout << "name subname "<<name<<" "<<subname<<endl;
      assert(0);
    }
    icorrshift_=2+idrinvbits+phi0bitshift-rinvbitshift-phiderbitshift;
    icorzshift_=idrinvbits-zderbitshift-tbitshift;
    phi0shift_=3;
    fact_=1;
    if (layer_>=4) {
      fact_=(1<<(nbitszprojL123-nbitszprojL456));
      icorrshift_-=2;
      icorzshift_+=(nbitszprojL123-nbitszprojL456+1);
      phi0shift_=0;
    }
    
    if (layer_==1){
      phimatchcut_[0]=-1;
      zmatchcut_[0]=-1;
      phimatchcut_[1]=0.08/(kphi1*rmean[layer_-1]);
      zmatchcut_[1]=7.0/kz;
      phimatchcut_[2]=0.08/(kphi1*rmean[layer_-1]);
      zmatchcut_[2]=14.0/kz;
    }
    if (layer_==2){
      phimatchcut_[0]=-1;
      zmatchcut_[0]=-1;
      phimatchcut_[1]=0.08/(kphi1*rmean[layer_-1]);
      zmatchcut_[1]=4.0/kz;
      phimatchcut_[2]=0.08/(kphi1*rmean[layer_-1]);
      zmatchcut_[2]=13.0/kz;
    }
    if (layer_==3){
      phimatchcut_[0]=0.1/(kphi1*rmean[layer_-1]);
      zmatchcut_[0]=0.7/kz;
      phimatchcut_[1]=-1;
      zmatchcut_[1]=-1;
      phimatchcut_[2]=0.08/(kphi1*rmean[layer_-1]);
      zmatchcut_[2]=10.0/kz;
    }
    if (layer_==4){
      phimatchcut_[0]=0.19/(kphi1*rmean[layer_-1]);
      zmatchcut_[0]=4.0/kz;
      phimatchcut_[1]=-1;
      zmatchcut_[1]=-1;
      phimatchcut_[2]=0.05/(kphi1*rmean[layer_-1]);
      zmatchcut_[2]=8.0/kz;
    }
    if (layer_==5){
      phimatchcut_[0]=0.4/(kphi1*rmean[layer_-1]);
      zmatchcut_[0]=4.0/kz;
      phimatchcut_[1]=0.08/(kphi1*rmean[layer_-1]);
      zmatchcut_[1]=9.0/kz;
      phimatchcut_[2]=-1;
      zmatchcut_[2]=-1;
    }
    if (layer_==6){
      phimatchcut_[0]=0.5/(kphi1*rmean[layer_-1]);
      zmatchcut_[0]=4.0/kz;
      phimatchcut_[1]=0.19/(kphi1*rmean[layer_-1]);
      zmatchcut_[1]=14.0/kz;
      phimatchcut_[2]=-1;
      zmatchcut_[2]=-1;
    }

    for(int iseedindex=0;iseedindex<7;iseedindex++){
      rphicutPS_[iseedindex]=-1.0;
      rphicut2S_[iseedindex]=-1.0;
      rcutPS_[iseedindex]=-1.0;
      rcut2S_[iseedindex]=-1.0;
    }

    if (abs(disk_)==1){
      rphicutPS_[0]=0.3;
      rcutPS_[0]=0.5;
      rphicut2S_[0]=0.5;
      rcut2S_[0]=3.8;      

      rphicut2S_[1]=0.8;
      rcut2S_[1]=3.8;      

      rphicutPS_[4]=0.15;
      rcutPS_[4]=0.5;

    }

    if (abs(disk_)==2){
      rphicutPS_[0]=0.2;
      rcutPS_[0]=0.5;
      rphicut2S_[0]=0.5;
      rcut2S_[0]=3.8;      

      rphicut2S_[1]=0.8;
      rcut2S_[1]=3.8;      

      rphicutPS_[4]=0.15;
      rcutPS_[4]=0.5;

      rphicutPS_[5]=0.15;
      rcutPS_[5]=0.5;
      rphicut2S_[5]=0.5;
      rcut2S_[5]=3.8;      

      rphicutPS_[6]=0.15;
      rcutPS_[6]=0.5;
      rphicut2S_[6]=0.15;
      rcut2S_[6]=3.4;      
    }

    if (abs(disk_)==3){
      rphicutPS_[0]=0.25;
      rcutPS_[0]=0.5;
      rphicut2S_[0]=0.5;
      rcut2S_[0]=3.6;      


      rphicutPS_[3]=0.15;
      rcutPS_[3]=0.5;
      rphicut2S_[3]=0.15;
      rcut2S_[3]=3.6;      

      rphicutPS_[5]=0.2;
      rcutPS_[5]=0.6;
      rphicut2S_[5]=0.2;
      rcut2S_[5]=3.6;

      rphicutPS_[6]=0.2;
      rcutPS_[6]=0.8;
      rphicut2S_[6]=0.5;
      rcut2S_[6]=3.8;      
    }


    if (abs(disk_)==4){
      rphicutPS_[0]=0.5;
      rcutPS_[0]=0.5;      
      rphicut2S_[0]=0.5;
      rcut2S_[0]=3.6;      


      rphicutPS_[3]=0.2;
      rcutPS_[3]=0.5;
      rphicut2S_[3]=0.2;
      rcut2S_[3]=3.6;      

      rphicutPS_[5]=0.4;
      rcutPS_[5]=1.0;
      rphicut2S_[5]=0.4;
      rcut2S_[5]=3.5;      

      rphicut2S_[6]=0.5;
      rcut2S_[6]=3.8;      
    }



    if (abs(disk_)==5){
      rphicutPS_[3]=0.25;
      rcutPS_[3]=0.75;
      rphicut2S_[3]=0.4;
      rcut2S_[3]=3.6;      

      rphicutPS_[4]=0.15;
      rcutPS_[4]=0.5;
      rphicut2S_[4]=0.2;
      rcut2S_[4]=3.4;      

      rphicutPS_[5]=0.5;
      rcutPS_[5]=2.0;
      rphicut2S_[5]=0.4;
      rcut2S_[5]=3.7;      


    }

    
  }

  void addOutput(FPGAMemoryBase* memory,string output){
    if (writetrace) {
      cout << "In "<<name_<<" adding output to "<<memory->getName()
	   << " to output "<<output<<endl;
    }
    if (output=="matchout1"||
	output=="matchout2"||
	output=="matchout3"||
	output=="matchout4"||
	output=="matchout5"||
	output=="matchout6"||
	output=="matchout7"||
	output=="matchout8"
	){
      FPGAFullMatch* tmp=dynamic_cast<FPGAFullMatch*>(memory);
      assert(tmp!=0);
      fullmatches_.push_back(tmp);
      return;
    }
    if (output=="matchoutplus"){
      FPGAFullMatch* tmp=dynamic_cast<FPGAFullMatch*>(memory);
      assert(tmp!=0);
      assert(fullmatchesToPlus_==0);
      fullmatchesToPlus_=tmp;
      return;
    }
    if (output=="matchoutminus"){
      FPGAFullMatch* tmp=dynamic_cast<FPGAFullMatch*>(memory);
      assert(tmp!=0);
      assert(fullmatchesToMinus_==0);
      fullmatchesToMinus_=tmp;
      return;
    }
    cout << "Count not fined output = "<<output<<endl;
    assert(0);
  }

  void addInput(FPGAMemoryBase* memory,string input){
    if (writetrace) {
      cout << "In "<<name_<<" adding input from "<<memory->getName()
	   << " to input "<<input<<endl;
    }
    if (input=="allstubin"){
      FPGAAllStubs* tmp=dynamic_cast<FPGAAllStubs*>(memory);
      assert(tmp!=0);
      allstubs_=tmp;
      return;
    }
    if (input=="allprojin"){
      FPGAAllProjections* tmp=dynamic_cast<FPGAAllProjections*>(memory);
      assert(tmp!=0);
      allprojs_=tmp;
      return;
    }
    if (input=="match1in"||
	input=="match2in"||
	input=="match3in"||
	input=="match4in"||
	input=="match5in"||
	input=="match6in"||
	input=="match7in"||
	input=="match8in"){
      FPGACandidateMatch* tmp=dynamic_cast<FPGACandidateMatch*>(memory);
      assert(tmp!=0);
      matches_.push_back(tmp);
      return;
    }
    assert(0);
  }

  void execute() {

    
    //Check that order is OK
    checkOrder();
    
    assert(fullmatches_.size()!=0);

    unsigned int countall=0;
    unsigned int countsel=0;

    FPGATracklet* oldTracklet=0;

    std::vector<std::pair<FPGATracklet*,std::pair<FPGAStub*,L1TStub*> > > mergedMatches=mergeMatches(matches_);
    
    for(unsigned int j=0;j<mergedMatches.size();j++){
	
      if (debug1&&j==0) {
        cout << getName() <<" has "<<mergedMatches.size()<<" candidate matches"<<endl;
      }
      
      countall++;
	
      L1TStub* stub=mergedMatches[j].second.second;
      FPGAStub* fpgastub=mergedMatches[j].second.first;
      FPGATracklet* tracklet=mergedMatches[j].first;

      if (oldTracklet!=0) {
	//allow equal here since we can have more than one cadidate match per tracklet projection
	if (iSector_==tracklet->homeSector()&&
	    iSector_==oldTracklet->homeSector()) {
	  assert(oldTracklet->TCID()<=tracklet->TCID());
	}
      }
      oldTracklet=tracklet;
      
      if (layer_!=0) {
	  
	int seedlayer=tracklet->layer();
	
	double phi=stub->phi();
	if (phi<0) phi+=two_pi;
	phi-=phioffset_;
	
	double dr=stub->r()-tracklet->rproj(layer_);
      	assert(fabs(dr)<drmax);

	int ir=fpgastub->r().value();
       	int iphi=tracklet->fpgaphiproj(layer_).value();

	int icorr=(ir*tracklet->fpgaphiprojder(layer_).value())>>icorrshift_;

	iphi+=icorr;
	
	int iz=tracklet->fpgazproj(layer_).value();

	
	int izcor=(ir*tracklet->fpgazprojder(layer_).value())>>icorzshift_;
	
	//if (seedlayer==1&&layer_==6) {
	//  cout << "layer iphi iphicor : "<<layer_<<" "<<tracklet->fpgaphiproj(layer_).value()*kphi1<<" "<<icorr*kphi1
	//       <<"   FP : "<<tracklet->phiproj(layer_)<<" "<<dr*tracklet->phiprojder(layer_)<<endl;
	//}

	iz+=izcor;	
	
	int ideltaz=fpgastub->z().value()-iz;

	int ideltaphi=(fpgastub->phi().value()<<phi0shift_)-(iphi<<(phi0bitshift-1+phi0shift_)); 
	
	if (dumpmatch) {
	  int l=1;
	  if (layer_>3) {
	    l=8; 
	  }
	  
	  cout << "DUMPMATCH : layer = "<<layer_<<endl;
	  cout << "DUMPMATCH : phiproj  = "<<tracklet->phiproj(layer_)<<" "
	       << tracklet->fpgaphiproj(layer_).value()*kphiproj123/l<<endl;
	  cout << "DUMPMATCH : phi corr = "<< dr*tracklet->phiprojder(layer_)
	       << " "<<icorr*kphiproj123/l<<endl;
	  cout << "DUMPMATCH : dr       = "<<dr<<" "<<ir*kr<<endl;
	  cout << "DUMPMATCH : zproj    = "<<tracklet->zproj(layer_)<<" "
	       <<  tracklet->fpgazproj(layer_).value()*kz*fact_<<" "
	       <<  tracklet->zproj(layer_)-tracklet->fpgazproj(layer_).value()*kz*fact_
	       << endl;
	  cout << "DUMPMATCH : z der   = "<<layer_<<" "<<tracklet->zprojder(layer_)
	       << " "<<tracklet->fpgazprojder(layer_).value()*kz/(kr*(1<<(idrinvbits-zderbitshift-tbitshift)))
	       <<" ("<<tracklet->fpgazprojder(layer_).value()<<")"<<endl;  
	  cout << "DUMPMATCH : z icorr corr  = "<<layer_<< " "<<izcor*kz*fact_<<" "<<dr*tracklet->zprojder(layer_)<<endl;
	  double phistub=stub->phi();
	  double dphi=two_pi/NSector;
	  phistub+=dphi/6.0;
	  for(;phistub>dphi;phistub-=dphi);
	  for(;phistub<0;phistub+=dphi);
	  cout << "DUMPMATCH : phi stub project   = "<<phistub<<" "<<tracklet->phiproj(layer_)<<endl;
	  cout << "DUMPMATCH : deltaphi   = "<<ideltaphi*kphi1*rmean[layer_-1]*10000.0<<" "
	       << " "<<(phistub-(tracklet->phiproj(layer_)+dr*tracklet->phiprojder(layer_)))*rmean[layer_-1]*10000.0<<endl;
	  cout << "DUMPMATCH : ideltaz deltaz   = "<<layer_<<" "<<ideltaz*kz<<" "<<stub->z()-(tracklet->zproj(layer_)+
											      dr*tracklet->zprojder(layer_))<<endl;
	  
	}
	
	assert(fabs(dr)<drmax);

	double dphi=phi-(tracklet->phiproj(layer_)+
			 dr*tracklet->phiprojder(layer_));
	  
	if (dphi>0.5*two_pi) dphi-=two_pi;
	if (dphi<-0.5*two_pi) dphi+=two_pi;
	
	
	double dz=stub->z()-(tracklet->zproj(layer_)+
			     dr*tracklet->zprojder(layer_));
	
	double dphiapprox=phi-(tracklet->phiprojapprox(layer_)+
			       dr*tracklet->phiprojderapprox(layer_));
	
	if (dphiapprox>0.5*two_pi) dphiapprox-=two_pi;
	if (dphiapprox<-0.5*two_pi) dphiapprox+=two_pi;
	
	
	double dzapprox=stub->z()-(tracklet->zprojapprox(layer_)+
				   dr*tracklet->zprojderapprox(layer_));
	
	int seedindex=(seedlayer-1)/2;
	if (seedlayer==0) seedindex=1;  //FIXME not applying the correct cuts
	if (seedlayer==2) seedindex=1;  //FIXME not applying the correct cuts
	assert(seedindex>=0);
	assert(seedindex<=2);
	assert(phimatchcut_[seedindex]>0);
	assert(zmatchcut_[seedindex]>0);

	
	if (writeResiduals) {
	  static ofstream out("layerresiduals.txt");
	  
	  double pt=0.003*3.8/fabs(tracklet->rinv());
	  
	  out << layer_<<" "<<seedlayer<<" "<<pt<<" "
	      <<ideltaphi*kphi1*rmean[layer_-1]<<" "<<dphi*rmean[layer_-1]<<" "<<phimatchcut_[seedindex]*kphi1*rmean[layer_-1]
	      <<" "<<ideltaz*fact_*kz<<" "<<dz<<" "<<zmatchcut_[seedindex]*kz<<endl;	  
	}

	
	bool imatch=(fabs(ideltaphi)<=phimatchcut_[seedindex])&&(fabs(ideltaz*fact_)<=zmatchcut_[seedindex]);

	//cout << getName() << "barrel ideltaphi deltaphi "<<ideltaphi*kphi1*rmean[layer_-1]<<" "<<dphi*rmean[layer_-1]<<endl;

	
	if (imatch) {
	  
	  std::pair<FPGAStub*,L1TStub*> tmp(fpgastub,stub);
	  
	  countsel++;

	  tracklet->addMatch(layer_,ideltaphi,ideltaz,
			     dphi,dz,dphiapprox,dzapprox,
				 fpgastub->stubindex().value()+(1<<6)*(fpgastub->phiregion().value()),
			     stub->r(),tmp);
	  

	  if (debug1) {
	    cout << "Accepted full match in layer " <<getName()
		 << " "<<tracklet
		 << " "<<iSector_<<endl;	   
	  }
	      
	  if (tracklet->plusNeighbor(layer_)){
	    
	    assert(fullmatchesToMinus_!=0);
	    
	    int nmatch=fullmatchesToMinus_->nMatches();
	    if (nmatch>1) {
	      assert(fullmatchesToMinus_->getFPGATracklet(nmatch-2)->TCID()<
		     fullmatchesToMinus_->getFPGATracklet(nmatch-1)->TCID());
	    }
	    
	    fullmatchesToMinus_->addMatch(tracklet,tmp);
	  } else if (tracklet->minusNeighbor(layer_)) {
	    assert(fullmatchesToPlus_!=0);
	    fullmatchesToPlus_->addMatch(tracklet,tmp);
	  } else {
	    for (unsigned int l=0;l<fullmatches_.size();l++){
	      //cout << "Adding match to: "<<fullmatches_[l]->getName()<<endl;
	      if ((tracklet->layer()==1&&fullmatches_[l]->getName().substr(3,2)=="L1")||
		  (tracklet->layer()==3&&fullmatches_[l]->getName().substr(3,2)=="L3")||
		  (tracklet->layer()==5&&fullmatches_[l]->getName().substr(3,2)=="L5")){
		assert(tracklet->homeSector()==iSector_);
		fullmatches_[l]->addMatch(tracklet,tmp);
	      }
	    }
	  }
	  
	}
	
	} else {
	
	
	//check that stubs and projections in same half of detector
	assert(stub->z()*tracklet->t()>0.0);
	
	int disk=disk_;
	if (tracklet->t()<0) disk=-disk_;
	
	assert(disk!=0);
	  
	double phi=stub->phi();
	if (phi<0) phi+=two_pi;
	phi-=phioffset_;
	
	double dz=stub->z()-tracklet->zprojdisk(disk);
	
	assert(fabs(dz)<dzmax);
	
	int iz=fpgastub->z().value();
	  
	int iphi=tracklet->fpgaphiprojdisk(disk).value();

	int shifttmp=t2bits+tbitshift+phi0bitshift+2-rinvbitshiftdisk-phiderdiskbitshift;

	assert(shifttmp>=0);
	int iphicorr=(iz*tracklet->fpgaphiprojderdisk(disk).value())>>shifttmp;

	iphi+=iphicorr;
	  
	double phicorr=dz*tracklet->phiprojderdisk(disk);
	
	assert(fabs(tracklet->phiprojderdisk(disk))<0.1);
	assert(fabs(phicorr)<0.1);
	
	double phiproj=tracklet->phiprojdisk(disk)+phicorr;
	
	int ir=tracklet->fpgarprojdisk(disk).value()+rmindisk/krprojshiftdisk;
	
	  
	int shifttmp2=rprojdiskbitshift+t3shift-rderdiskbitshift;
	
	assert(shifttmp2>=0);
	int ircorr=(iz*tracklet->fpgarprojderdisk(disk).value())>>shifttmp2;
	  
	int dorcorr=1;
	
	ir+=ircorr*dorcorr;

	double rcorr=dz*tracklet->rprojderdisk(disk);
	
	double rproj=tracklet->rprojdisk(disk)+rcorr*dorcorr;
	
	int ideltaphi=fpgastub->phi().value()*kphi/kphiproj123-iphi; 
	
	double deltar=stub->r()-rproj;
	  
	int irstub = fpgastub->r().value();
	if(!stub->isPSmodule()){
	  double rstub = rDSS[irstub];
	  assert (rstub < rmaxdisk);
	  irstub = (1<<nrbitsdisk)*(rstub)/(rmaxdisk-rmindisk);
	  } else {
	  irstub+=rmindisk/kr;
	}
	
	int ideltar=(irstub*krdisk)/krprojshiftdisk-ir;
	
	double dr=stub->r()-(tracklet->rprojdisk(disk)+
			     dz*tracklet->rprojderdisk(disk)*dorcorr);
	
	double dphi=phi-(tracklet->phiprojdisk(disk)+
			 dz*tracklet->phiprojderdisk(disk));
	
	if (dphi>0.5*two_pi) dphi-=two_pi;
	if (dphi<-0.5*two_pi) dphi+=two_pi;

	if (dphi>0.5*two_pi/NSector) dphi-=two_pi/NSector;
	if (dphi<-0.5*two_pi/NSector) dphi+=two_pi/NSector;
	
	double dphiapprox=phi-(tracklet->phiprojapproxdisk(disk)+
			       dz*tracklet->phiprojderapproxdisk(disk));
	  
	if (dphiapprox>0.5*two_pi) dphiapprox-=two_pi;
	if (dphiapprox<-0.5*two_pi) dphiapprox+=two_pi;

	if (dphiapprox>0.5*two_pi/NSector) dphiapprox-=two_pi/NSector;
	if (dphiapprox<-0.5*two_pi/NSector) dphiapprox+=two_pi/NSector;
	
	  
	double drapprox=stub->r()-(tracklet->rprojapproxdisk(disk)+
				   dz*tracklet->rprojderapproxdisk(disk)*dorcorr);
	
	double alpha=0.0;

	double drphi=dphi*stub->r();
	double drphiapprox=dphiapprox*stub->r();
	
	if (!stub->isPSmodule()) {
	  alpha=stub->alpha(); 	
	  dphi+=dr*alpha;
	  dphiapprox+=drapprox*alpha;
	  
	  double alphanew=stub->alphanew();
	  drphi+=dr*alphanew*4.57/stub->r();
	  drphiapprox+=dr*alphanew*4.57/stub->r();

	  int ialphanew=fpgastub->alphanew().value();

	  int alphashift=12;
	  double fact=(1<<alphashift)*krprojshiftdisk*4.57/(1<<(nbitsalpha-1))/stub->r2()/kphiproj123;
	  int ifact=fact;
	  
	  //cout << "fact ifact "<<fact<<" "<<ifact<<endl;
	  //cout << "alphanew  "<<alpha<<" "<<alphanew<<" "<<alphanew*4.57/stub->r2()<<" "<<ialphanew*4.57/(1<<(nbitsalpha-1))/stub->r2()<<endl;
	  //cout << "rdphi corr  "<<dr*alpha*stub->r()<<" "<<dr*alphanew*4.57/stub->r()<<" "<<ideltar*krprojshiftdisk*ialphanew*4.57/(1<<(nbitsalpha-1))/stub->r()<<endl;
	  //cout << "dphi drphi "<<dphi*stub->r()<<" "<<drphi<<endl;
	  //cout << "dr ideltar : "<<dr<<" "<<ideltar*krprojshiftdisk<<endl;

	  int iphialphacor=((ideltar*ialphanew*ifact)>>alphashift);
	  
	  //cout << "ialphacor "<<iphialphacor<<" "<<ideltar*fpgastub->alpha().value()*krprojshiftdisk*kalpha/kphiproj123<<endl;
	  
	  //ideltaphi+=ideltar*fpgastub->alpha().value()*krprojshiftdisk*kalpha/kphiproj123;
	  ideltaphi+=iphialphacor;
	}



	int seedlayer=tracklet->layer();
	int seeddisk=tracklet->disk();

	int seedindex=-1;

	if (seedlayer==1&&seeddisk==0) seedindex=0;  //L1L2
	if (seedlayer==3&&seeddisk==0) seedindex=1;  //L3L4
	if (seedlayer==5&&seeddisk==0) seedindex=2;  //L5L6
	if (seedlayer==0&&abs(seeddisk)==1) seedindex=3;  //D1D2
	if (seedlayer==0&&abs(seeddisk)==3) seedindex=4;  //D3D4
	if (seedlayer==1&&abs(seeddisk)==1) seedindex=5;  //L1D1
	if (seedlayer==2&&abs(seeddisk)==1) seedindex=6;  //L2D1

	if (seedindex<0) {
	  cout << "seedlayer abs(seeddisk) : "<<seedlayer<<" "<<abs(seeddisk)<<endl;
	  assert(0);
	}
	
	double drphicut=rphicutPS_[seedindex];
	double drcut=rcutPS_[seedindex]; 
	if (!stub->isPSmodule()) {
	  drphicut=rphicut2S_[seedindex];
	  drcut=rcut2S_[seedindex]; 
	}
	

	if (drphicut<0.0 || drcut<0.0) {
	  cout << "drphicut drcut : "<<drphicut<<" "<<drcut<<endl;
	  cout << "disk_ isPS seedindex : "<<disk_<<" "<<stub->isPSmodule()<<" "<<seedindex<<endl;
	  assert(0);
	}
	
	if (writeResiduals) {
	  static ofstream out("diskresiduals.txt");
	  
	  double pt=0.003*3.8/fabs(tracklet->rinv());
	  
	  out << disk_<<" "<<stub->isPSmodule()<<" "<<seedlayer<<" "<<abs(seeddisk)<<" "<<pt<<" "
	      <<ideltaphi*kphiproj123*stub->r()<<" "<<drphi<<" "
	    //<<phimatchcut_[seedindex]*kphi1*rmean[layer_-1]<<" "
	      <<drphicut<<" "
	      <<ideltar*krprojshiftdisk<<" "<<deltar<<" "
	    //<<zmatchcut_[seedindex]*kz
	      <<drcut<<" "
	      <<endl;	  
	}

	
	
	bool match=(fabs(drphi)<drphicut)&&(fabs(deltar)<drcut);
	
	bool imatch=(fabs(ideltaphi)<drphicut/(kphiproj123*stub->r()))&&(fabs(ideltar)<drcut/krprojshiftdisk);
	
	//cout << getName() << " ideltaphi deltaphi "<<ideltaphi*kphiproj123*stub->r()<<" "<<dphi*stub->r()<<endl;
	
	if (debug1) {
	  cout << "imatch disk: "<<imatch<<" "<<fabs(ideltaphi)<<" "<<drphicut/(kphiproj123*stub->r())<<" "
	       <<fabs(ideltar)<<" "<<drcut/krprojshiftdisk<<" r = "<<stub->r()<<endl;
	}
	
	if (writeDiskMatch1) {

	  static ofstream out1("diskmatch1.txt");
	  
	  out1 << disk<<" "
	       << phiproj<<" "
	       << rproj<<" "
	       << dphi<<" "
	       << deltar<<"    "
	       << iphi*kphiprojdisk<<" "
	       << ir*krprojshiftdisk<<"  "
	       << ideltaphi*kphiprojdisk<<" "
	       << ideltar*krprojshiftdisk<<" "
	       << endl;
	  
	}
	
	if (dumpmatch) {
	  
	  cout << "DUMPMATCHDISK "<<tracklet->layer()<<" "<<tracklet->disk()<<" "<<stub->r()
	       <<" dr (float) = "<<deltar<<" "
	       <<" dr (int) = "<<ideltar*krprojshiftdisk
	       <<" dr res = "<<krprojshiftdisk
	       <<" rcorr = "<<rcorr
	       <<" ircorr = "<<ircorr*krprojshiftdisk
	       <<endl;
	  
	  cout << "DUMPMATCHDISK rproj(float) = "<<tracklet->rprojdisk(disk)
	       <<" rproj(int) = "<<(tracklet->fpgarprojdisk(disk).value()+rmindisk/krprojshiftdisk)*krprojshiftdisk
	       <<endl;
	  
	  cout << "DUMPMATCHDISK1 : "<<disk<<" "
	       <<tracklet->phiprojdisk(disk)<<" "
	       <<dz<<" "
	       <<phicorr<<" "
	       <<phiproj<<" | "
	       <<tracklet->rprojdisk(disk)<<" "
	       <<tracklet->rprojderdisk(disk)<<" "
	       <<rcorr<<" "
	       <<rproj<<" dphi= "
	       <<dphi<<" dr= "
	       <<deltar<<" "
	       <<match
	       <<endl;
	  cout << "DUMPMATCHDISK2 : "<<disk<<" "
	       <<tracklet->fpgaphiprojdisk(disk).value()*kphiproj123<<" "
	       <<iz*kz<<" "
	       <<iphicorr*kphiproj123<<" ("<<iphicorr<<") "
	       <<iphi*kphiproj123<<" | "
	       <<tracklet->fpgarprojdisk(disk).value()*krprojshiftdisk<<" "
	       <<tracklet->fpgarprojderdisk(disk).value()*krprojderdiskshift<<" "
	       <<ircorr*krprojshiftdisk<<" "
	       <<ir*krprojshiftdisk<<" dphi= "
	       <<ideltaphi*kphiproj123<<" dr= "
	       <<ideltar*krprojshiftdisk<<" "
	       <<imatch
	       <<endl;
	}
	
	  
	if (imatch) {
	  
	  std::pair<FPGAStub*,L1TStub*> tmp(fpgastub,stub);
	    
	  countsel++;
	  
	  if (debug1) {
	    cout << "FPGAMatchCalculator found match in disk "<<getName()<<endl;
	  }
	  
	  assert(fabs(dphi)<0.2);
	  assert(fabs(dphiapprox)<0.2);
	  
	  tracklet->addMatchDisk(disk,ideltaphi,ideltar,
				 drphi/stub->r(),dr,drphiapprox/stub->r(),drapprox,
				 stub->alpha(),
				 fpgastub->stubindex().value()+((fpgastub->phiregion().value())<<6),
				 stub->z(),tmp);
	  if (debug1) {
	    cout << "Accepted full match in disk " <<getName()
		 << " "<<tracklet
		 << " "<<iSector_<<endl;	   
	  }
	  
	  if (tracklet->plusNeighborDisk(disk)){
	    fullmatchesToMinus_->addMatch(tracklet,tmp);
	    if (debug1) {
	      cout << "Accepted full match to minus in disk " <<getName()<<" "<<tracklet
		   <<" "<<fullmatchesToMinus_->getName()<<endl;
	    }
	    int nmatch=fullmatchesToMinus_->nMatches();
	    if (nmatch>1) {
	      assert(fullmatchesToMinus_->getFPGATracklet(nmatch-2)->TCID()<
		     fullmatchesToMinus_->getFPGATracklet(nmatch-1)->TCID());
	    }
	  } else if (tracklet->minusNeighborDisk(disk)) {
	    fullmatchesToPlus_->addMatch(tracklet,tmp);
	    if (debug1) {
	      cout << "Accepted full match to plus in disk " <<getName()<<" "<<tracklet
		   <<" "<<fullmatchesToPlus_->getName()<<endl;
	    }
	  } else {
	    for (unsigned int l=0;l<fullmatches_.size();l++){
	      if (((abs(tracklet->disk())==1&&tracklet->layer()==1)&&fullmatches_[l]->getName().substr(3,4)=="D1L1")||
		  (tracklet->layer()==2&&fullmatches_[l]->getName().substr(3,4)=="D1L2")||    //dangerous to check only layer!!!
		  ((abs(tracklet->disk())==1&&tracklet->layer()==0)&&fullmatches_[l]->getName().substr(3,4)=="D1D2")||
		  ((tracklet->disk()==0&&tracklet->layer()==1)&&fullmatches_[l]->getName().substr(3,4)=="L1L2")||
		  ((tracklet->disk()==0&&tracklet->layer()==3)&&fullmatches_[l]->getName().substr(3,4)=="L3L4")||
		  ((abs(tracklet->disk())==3&&tracklet->layer()==0)&&fullmatches_[l]->getName().substr(3,4)=="D3D4")){
		fullmatches_[l]->addMatch(tracklet,tmp);
		if (debug1) {
		  cout << "In "<<getName()<<" added match to "<<fullmatches_[l]->getName()<<endl;
		}
	      }
	    }
	  }
	}
      }
      if (countall>=MAXMC) break;
    }


    if (writeMatchCalculator) {
      static ofstream out("matchcalculator.txt");
      out << getName()<<" "<<countall<<" "<<countsel<<endl;
    }

    
  }


    std::vector<std::pair<FPGATracklet*,std::pair<FPGAStub*,L1TStub*> > > mergeMatches(vector<FPGACandidateMatch*>& candmatch) {

    std::vector<std::pair<FPGATracklet*,std::pair<FPGAStub*,L1TStub*> > >  tmp;

    std::vector<unsigned int> indexArray;
    indexArray.reserve(candmatch.size());
    for (unsigned int i=0;i<candmatch.size();i++) {
      indexArray.push_back(0);
    }

    int bestIndex=-1;
    do {
      int bestSector=100;
      int bestTCID=(1<<16);
      bestIndex=-1;
      for (unsigned int i=0;i<candmatch.size();i++) {
	if (indexArray[i]>=candmatch[i]->nMatches()) {
	  // skip as we were at the end
	  continue;
	}
	int TCID=candmatch[i]->getFPGATracklet(indexArray[i])->TCID();
	int dSector=candmatch[i]->getFPGATracklet(indexArray[i])->homeSector()-iSector_;
	if (dSector>2) dSector-=NSector;
	if (dSector<-2) dSector+=NSector;
	assert(abs(dSector)<2);
	if (dSector==-1) dSector=2;
	if (dSector<bestSector) {
	  bestSector=dSector;
	  bestTCID=TCID;
	  bestIndex=i;
	}
	if (dSector==bestSector) {
	  if (TCID<bestTCID) {
	    bestTCID=TCID;
	    bestIndex=i;
	  }
	}
      }
      if (bestIndex!=-1) {
	tmp.push_back(candmatch[bestIndex]->getMatch(indexArray[bestIndex]));
	indexArray[bestIndex]++;
      }
    } while (bestIndex!=-1);
    
    if (layer_>0) {
      
      int lastTCID=-1;
      int lastTCIDplus=-1;
      int lastTCIDminus=-1;
      bool error=false;
      
      //Allow equal TCIDs since we can have multiple candidate matches
      for(unsigned int i=1;i<tmp.size();i++){
	if (tmp[i].first->minusNeighbor(layer_)) {
	  //cout << "For minus: tracklet TCID "<<tracklet<<" "<<tracklet->TCID()<<" "<<inputproj_[j]->getName()<<endl;
	  if (lastTCIDminus>tmp[i].first->TCID()) {
	    cout << "Wrong TCID ordering for Minus projections in "<<getName()<<" last "<<lastTCIDminus<<" "<<tmp[i].first->TCID()<<endl;
	    error=true;
	  } else {
	    lastTCIDminus=tmp[i].first->TCID();
	  }
	} else if (tmp[i].first->plusNeighbor(layer_)) {
	  if (lastTCIDplus>tmp[i].first->TCID()) {
	    cout << "Wrong TCID ordering for Plus projections in "<<getName()<<" last "<<lastTCIDplus<<" "<<tmp[i].first->TCID()<<endl;
	    error=true;
	  } else {
	    lastTCIDplus=tmp[i].first->TCID();
	      }
	} else {
	  if (lastTCID>tmp[i].first->TCID()) {
	    cout << "Wrong TCID ordering for projections in "<<getName()<<" last "<<lastTCID<<" "<<tmp[i].first->TCID()<<endl;
	    error=true;
	  } else {
	    lastTCID=tmp[i].first->TCID();
	  }
	}
      }
      
      if (error) {
	for(unsigned int i=1;i<tmp.size();i++){
	  cout << "Wrong order for in "<<getName()<<" "<<i<<" "<<tmp[i].first<<" "<<tmp[i].first->TCID()<<" "
	       <<tmp[i].first->plusNeighbor(layer_)<<" "<<tmp[i].first->minusNeighbor(layer_)<<endl;
	}
      }
      
    }
    
    for (unsigned int i=0;i<tmp.size();i++) {
      if (i>0) {
	//This allows for equal TCIDs. This means that we can e.g. have a track seeded
	//in L1L2 that projects to both L3 and D4. The algorithm will pick up the first hit and
	//drop the second

	if (iSector_==tmp[i-1].first->homeSector()&&
	    iSector_==tmp[i].first->homeSector()) {
	  assert(tmp[i-1].first->TCID()<=tmp[i].first->TCID());
	}
      }
    }
    
    return tmp;
  }

  void checkOrder(){
  
    if (layer_>0) {
      for (unsigned int l=0;l<matches_.size();l++){
	for(unsigned int j=1;j<matches_[l]->nMatches();j++) {
	  bool firstMinus=matches_[l]->getMatch(j-1).first->minusNeighbor(layer_);
	  bool firstPlus=matches_[l]->getMatch(j-1).first->plusNeighbor(layer_);
	  bool firstCenter=!(firstMinus||firstPlus);
	  bool secondMinus=matches_[l]->getMatch(j).first->minusNeighbor(layer_);
	  bool secondPlus=matches_[l]->getMatch(j).first->plusNeighbor(layer_);
	  bool secondCenter=!(secondMinus||secondPlus);
	  if (((!firstCenter)&&secondCenter)||
	      (firstPlus&&secondMinus)){
	    cout << "Wrong order in  "<<matches_[l]->getName()<<" "<<matches_[l]->getMatch(j-1).first->plusNeighbor(layer_)<<" "<<matches_[l]->getMatch(j-1).first->minusNeighbor(layer_)<<"  -  "<<
	      matches_[l]->getMatch(j).first->plusNeighbor(layer_)<<" "<<matches_[l]->getMatch(j).first->minusNeighbor(layer_)<<endl;
	  }
	}
      } 
    }
  }

    
private:

  int layer_;
  int disk_;
  int fact_;
  int icorrshift_;
  int icorzshift_;
  int phi0shift_;
  int phimatchcut_[3];
  int zmatchcut_[3];
  double phimin_;
  double phimax_;
  double phioffset_;

  double rphicutPS_[7];
  double rphicut2S_[7];
  double rcutPS_[7];
  double rcut2S_[7];

  FPGAAllStubs* allstubs_;
  FPGAAllProjections* allprojs_;

  vector<FPGACandidateMatch*> matches_;

  vector<FPGAFullMatch*> fullmatches_;
  FPGAFullMatch* fullmatchesToPlus_;
  FPGAFullMatch* fullmatchesToMinus_;

};

#endif
