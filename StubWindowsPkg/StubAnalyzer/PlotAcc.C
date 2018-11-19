#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TError.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <iomanip>
#include <cmath>
#include <fstream>

void SetPlotStyle();

void PlotAcc(TString sourcefile, TString outputname, int pdid=11, int nhits=0, int nevt=99999, float ptmin=2.0, float ptmax=100.0, float etamax=2.4)
{
  // pdid     PDG ID of particle
  //            13  = single muon+/muon-
  //            11  = single electron/positron
  //            211 = single pion+/pion-
  // nhits    how many layers/disks are hit by the particle?
  // nevt     change if you want to look at a # of events less than the total

  
  // -----------------------------------------------------------------------------------------------------------
  // Preliminaries
  // -----------------------------------------------------------------------------------------------------------

  gROOT->SetBatch();
  gErrorIgnoreLevel = kError;
  SetPlotStyle();
  
  TString srcDIR  = "RootFiles/20181102/"; 
  TString subDIR  = "RootFiles/20181102/";
  TString plotDIR = "Plots/20181102/"; 
  TString plotPRE = "detAcc_";
  
  // Read in File
  TChain *FullI = new TChain("FullInfo");
  FullI->Add(srcDIR+sourcefile+".root"); 

  // Output file
  TFile* fout = new TFile(subDIR+"plots_"+plotPRE+outputname+".root","recreate");
  

  // -----------------------------------------------------------------------------------------------------------
  // Definitions
  // -----------------------------------------------------------------------------------------------------------

  // -----------------------------------------------------------------------------------------------------------
  // Variables
  int n_part;                           // The total number of particles inducing at least one stub in the event

  std::vector<int>     *part_pdg;       // PDG id of the particles
  std::vector<int>     *part_nstubs;    // How many stubs are induced by the particle in the tracker?
  std::vector<int>     *part_nhits;     // How many different layers/disks are hit by the particle?
  std::vector<int>     *part_nhits_fe;  // How many different layers/disks are hit by the particle and pass the FE cuts?

  std::vector<float>   *part_pt;     // pt of the particles
  std::vector<float>   *part_rho;    // rho0 of the particles
  std::vector<float>   *part_d0;     // d0 of the particles
  std::vector<float>   *part_z0;     // z0 of the particles
  std::vector<float>   *part_eta;    // eta of the particles
  std::vector<float>   *part_phi;    // phi of the particles
  std::vector<float>   *part_dist;

  // -----------------------------------------------------------------------------------------------------------
  // Initialize
  part_pdg=0;
  part_nstubs=0;
  part_nhits=0;
  part_nhits_fe=0;

  part_pt=0;
  part_rho=0;
  part_d0=0;
  part_z0=0;
  part_eta=0;
  part_phi=0;
  part_dist=0;

  // -----------------------------------------------------------------------------------------------------------
  // Constants
  float min_z0  = -15.;
  float max_z0  = 15.;
  float min_d0  = -5.;
  float max_d0  = 5.;
  
  float rhomax  = 100.0;
  float domax   = 5.0;

  float min_pt  = ptmin;
  float min_eta = -etamax;

  float max_pt  = ptmax;
  float max_eta = etamax;


  // -----------------------------------------------------------------------------------------------------------
  // Branch Addresses 
  TEfficiency *myEff = new TEfficiency();

  FullI->SetBranchAddress("n_part",       &n_part);
  FullI->SetBranchAddress("part_pdg",     &part_pdg);
  FullI->SetBranchAddress("part_nstubs",  &part_nstubs);
  FullI->SetBranchAddress("part_nhits",   &part_nhits);
  FullI->SetBranchAddress("part_nhits_fe",&part_nhits_fe);

  FullI->SetBranchAddress("part_pt",      &part_pt);
  FullI->SetBranchAddress("part_rho",     &part_rho);
  FullI->SetBranchAddress("part_d0",      &part_d0);
  FullI->SetBranchAddress("part_z0",      &part_z0);
  FullI->SetBranchAddress("part_eta",     &part_eta);
  FullI->SetBranchAddress("part_phi",     &part_phi);
  FullI->SetBranchAddress("part_dist",    &part_dist);

  // Initialize
  int n_entries = FullI->GetEntries();

  


  // ---------------------------------------------------------------------------------------------------------
  //        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
  // ---------------------------------------------------------------------------------------------------------


  // -----------------------------------------------------------------------------------------------------------
  // Detector Acceptance
  // -----------------------------------------------------------------------------------------------------------
    const int nbins_det = 25;
    
    float bin_pt_det  = (max_pt-min_pt)/nbins_det;
    float bin_eta_det = (max_eta-min_eta)/nbins_det;
    float bin_z0_det  = (max_z0-min_z0)/nbins_det;
    float bin_d0_det  = (max_d0-min_d0)/nbins_det;

    // Histograms
    TH2F *h_eff_pt_det  = new TH2F("eff_pt_det","eff_pt",nbins_det,min_pt,max_pt,200,0.1,1.05);
    TH2F *h_eff_eta_det = new TH2F("eff_eta_det","eff_eta",nbins_det,min_eta,max_eta,200,0.1,1.05);
    TH2F *h_eff_z_det   = new TH2F("eff_z_det","eff_z",nbins_det,min_z0,max_z0,200,0.1,1.05);
    TH2F *h_eff_d_det   = new TH2F("eff_d_det","eff_d",nbins_det,min_d0,max_d0,200,0.1,1.05);

    float eff_val_det[5][20][nbins_det];
    float sig_val_det[5][20][nbins_det];
    float count_det[5][20][nbins_det];
    
    float eff_vals_det[5][20];
    
    float abscissa_det[5][2][nbins_det];
    
    for (int k=0;k<5;++k) {
      for (int l=0;l<20;++l) {
        eff_vals_det[k][l]  = 0.;
      } // close l loop
    } // close k loop
    
    for (int j=0;j<nbins_det;++j) {
      abscissa_det[0][0][j] = min_pt+(j+0.5)*bin_pt_det;
      abscissa_det[1][0][j] = min_eta+(j+0.5)*bin_eta_det;
      abscissa_det[3][0][j] = min_z0+(j+0.5)*bin_z0_det;
      abscissa_det[4][0][j] = min_d0+(j+0.5)*bin_d0_det;
        
      for (int k=0;k<5;++k) {
        abscissa_det[k][1][j]  = 0.;
          
        for (int l=0;l<20;++l) {
          eff_val_det[k][l][j]  = 0.;
          sig_val_det[k][l][j]  = 0.;
          count_det[k][l][j]    = 0.;
        } // close l loop    
      } // close k loop  
    } // close j loop


    // Definition of eff_val values
    // eff_val_det[][0] : all the particles in the acceptance with at least 1 stub
    // eff_val_det[][1] : all the particles in [0] with at least nhits stubs
    // eff_val_det[][2] : all the particles in [1] with stubs in at least nhits layer/disks
    // eff_val_det[][3] : all the particles in [2] with stubs passing FE in at least nhits layer/disks

    
    // -----------------------------------------------------------------------------------------------------------
    // Loop over Events and Particles
    int bin_det[5];
    
    float npart_det_0=0;
    float npart_det_1=0;
    float npart_det_2=0;
    float npart_det_3=0;
    float npart_det_4=0;
    
    // Event Loop
    for (int j=0;j<std::min(n_entries,nevt);++j) {
      if (j%1000==0) cout << j << endl;
      
      int n_comb_tot_det = 0;
      
      FullI->GetEntry(j); // Get the tree info
      
      // First,count the number of patterns in the sector
      int npatt_det=0;
      int npatt_full_det=0;
      int ptbin_det;
      int inacc_det;
      int sectype_det;
      
      // Loop over all the particles
      for (int k=0;k<n_part;++k) {
        if (abs(part_pdg->at(k))!=pdid && pdid!=-1) continue;
        if (part_nstubs->at(k)<1) continue;
        if (part_pt->at(k)<ptmin) continue;
        if (part_pt->at(k)>ptmax) continue;
        if (part_rho->at(k)>rhomax) continue;
        if (fabs(part_d0->at(k))>domax) continue;

        if (part_z0->at(k)>max_z0) continue;
        if (part_z0->at(k)<min_z0) continue;
        if (part_d0->at(k)>max_d0) continue;
        if (part_d0->at(k)<min_d0) continue;
        if (part_eta->at(k)>max_eta) continue;
        if (part_eta->at(k)<min_eta) continue;

        bin_det[0] = std::min(int((part_pt->at(k)-min_pt)/bin_pt_det),nbins_det-1);
        bin_det[1] = std::min(int((part_eta->at(k)-min_eta)/bin_eta_det),nbins_det-1);
        bin_det[2] = 0;
        bin_det[3] = std::min(int((part_z0->at(k)-min_z0)/bin_z0_det),nbins_det-1);
        bin_det[4] = std::min(int((part_d0->at(k)-min_d0)/bin_d0_det),nbins_det-1);
        
        for (int l=0;l<5;++l) ++eff_val_det[l][0][bin_det[l]];
        
        if (part_nstubs->at(k)<nhits) continue;
        
        for (int l=0;l<5;++l) ++eff_val_det[l][1][bin_det[l]];

        if (part_nhits->at(k)<nhits) continue;            

        for (int l=0;l<5;++l) ++eff_val_det[l][2][bin_det[l]];

        if (part_nhits_fe->at(k)<nhits) continue;            

        for (int l=0;l<5;++l) ++eff_val_det[l][3][bin_det[l]];
      } // End of loop over particles
    } // End of loop on events


    for (int j=0;j<nbins_det;++j) {
      
      for (int k=0;k<5;++k) {
          
        for (int l=0;l<20;++l) eff_vals_det[k][l] += eff_val_det[k][l][j];
        
        if (eff_val_det[k][0][j]>=1) eff_val_det[k][6][j] = eff_val_det[k][1][j]/eff_val_det[k][0][j]; // SR
        if (eff_val_det[k][1][j]>=1) eff_val_det[k][7][j] = eff_val_det[k][2][j]/eff_val_det[k][1][j]; // SRNHITS
        if (eff_val_det[k][2][j]>=1) eff_val_det[k][8][j] = eff_val_det[k][3][j]/eff_val_det[k][2][j]; // SRNHITS_FE
        if (eff_val_det[k][0][j]>=1) eff_val_det[k][9][j] = eff_val_det[k][3][j]/eff_val_det[k][0][j]; // ALL
        
        if (eff_val_det[k][0][j]>=1) sig_val_det[k][1][j] = sqrt(eff_val_det[k][6][j]*(1-eff_val_det[k][6][j])/eff_val_det[k][0][j]);
        if (eff_val_det[k][1][j]>=1) sig_val_det[k][2][j] = sqrt(eff_val_det[k][7][j]*(1-eff_val_det[k][7][j])/eff_val_det[k][1][j]);
        if (eff_val_det[k][2][j]>=1) sig_val_det[k][3][j] = sqrt(eff_val_det[k][8][j]*(1-eff_val_det[k][8][j])/eff_val_det[k][2][j]);            
        if (eff_val_det[k][0][j]>=1) sig_val_det[k][4][j] = sqrt(eff_val_det[k][9][j]*(1-eff_val_det[k][9][j])/eff_val_det[k][0][j]);            
      }
    }
    
    cout << eff_vals_det[0][0] << " particles pass the kinematic cuts and have at least one stub" << endl;
    cout << eff_vals_det[0][1] << " of them have at least " << nhits << " stubs" << endl;
    cout << eff_vals_det[0][2] << " of them have stubs in at least " << nhits << " layer disks" << endl;
    cout << eff_vals_det[0][3] << " of them have stubs passing FE cuts in at least " << nhits << " layer disks" << endl;
  
    // -----------------------------------------------------------------------------------------------------------
    // Plots Errors
    TGraphErrors  *pt_stubs_eff = new TGraphErrors(nbins_det,abscissa_det[0][0],eff_val_det[0][6],abscissa_det[0][1],sig_val_det[0][1]);
    TGraphErrors  *pt_tower_eff = new TGraphErrors(nbins_det,abscissa_det[0][0],eff_val_det[0][7],abscissa_det[0][1],sig_val_det[0][2]);
    TGraphErrors  *pt_AM_eff    = new TGraphErrors(nbins_det,abscissa_det[0][0],eff_val_det[0][8],abscissa_det[0][1],sig_val_det[0][3]);
    
    TGraphErrors  *eta_stubs_eff = new TGraphErrors(nbins_det,abscissa_det[1][0],eff_val_det[1][6],abscissa_det[1][1],sig_val_det[1][1]);
    TGraphErrors  *eta_tower_eff = new TGraphErrors(nbins_det,abscissa_det[1][0],eff_val_det[1][7],abscissa_det[1][1],sig_val_det[1][2]);
    TGraphErrors  *eta_AM_eff    = new TGraphErrors(nbins_det,abscissa_det[1][0],eff_val_det[1][8],abscissa_det[1][1],sig_val_det[1][3]);

    TGraphErrors  *z_stubs_eff  = new TGraphErrors(nbins_det,abscissa_det[3][0],eff_val_det[3][6],abscissa_det[3][1],sig_val_det[3][1]);
    TGraphErrors  *z_tower_eff  = new TGraphErrors(nbins_det,abscissa_det[3][0],eff_val_det[3][7],abscissa_det[3][1],sig_val_det[3][2]);
    TGraphErrors  *z_AM_eff     = new TGraphErrors(nbins_det,abscissa_det[3][0],eff_val_det[3][8],abscissa_det[3][1],sig_val_det[3][3]);

    TGraphErrors  *d_stubs_eff  = new TGraphErrors(nbins_det,abscissa_det[4][0],eff_val_det[4][6],abscissa_det[4][1],sig_val_det[4][1]);
    TGraphErrors  *d_tower_eff  = new TGraphErrors(nbins_det,abscissa_det[4][0],eff_val_det[4][7],abscissa_det[4][1],sig_val_det[4][2]);
    TGraphErrors  *d_AM_eff     = new TGraphErrors(nbins_det,abscissa_det[4][0],eff_val_det[4][8],abscissa_det[4][1],sig_val_det[4][3]);

    TGraphErrors  *pt_inclusive_eff  = new TGraphErrors(nbins_det,abscissa_det[0][0],eff_val_det[0][9],abscissa_det[0][1],sig_val_det[0][4]);    
    TGraphErrors  *eta_inclusive_eff = new TGraphErrors(nbins_det,abscissa_det[1][0],eff_val_det[1][9],abscissa_det[1][1],sig_val_det[1][4]);    
    TGraphErrors  *z_inclusive_eff   = new TGraphErrors(nbins_det,abscissa_det[3][0],eff_val_det[3][9],abscissa_det[3][1],sig_val_det[3][4]);    
    TGraphErrors  *d_inclusive_eff   = new TGraphErrors(nbins_det,abscissa_det[4][0],eff_val_det[4][9],abscissa_det[4][1],sig_val_det[4][4]);    

    // -----------------------------------------------------------------------------------------------------------
    // Canvases
    char buffer_det[80];
    
    float low_det,up_det;
  
    // -----------------------------------------------------------------------------------------------------------
    // Canvas 1: c_pteff
    TCanvas *c_pteff = new TCanvas("c_pteff","pT efficiencies",201,77,1170,608);
      c_pteff->Divide(3,1);
      
      for (int i=0;i<3;++i) {
        c_pteff->cd(i+1);
        c_pteff->cd(i+1)->Range(-5.887851,-1.930603,70.65421,17.37543);
        c_pteff->cd(i+1)->SetFillColor(0);
        c_pteff->cd(i+1)->SetBorderMode(0);
        c_pteff->cd(i+1)->SetBorderSize(2);
        c_pteff->cd(i+1)->SetGridx();
        c_pteff->cd(i+1)->SetGridy();
        c_pteff->cd(i+1)->SetLeftMargin(0.13); //0.07692308
        c_pteff->cd(i+1)->SetTopMargin(0.065); //0.07124352
        c_pteff->cd(i+1)->SetFrameBorderMode(0);
        h_eff_pt_det->GetXaxis()->SetTitle("Tracking particle p_{T} (in GeV/c)");
        h_eff_pt_det->GetYaxis()->SetTitle("Efficiency");
        h_eff_pt_det->Draw();
        
        if (i==0) {
          pt_stubs_eff->SetMarkerStyle(20);
          pt_stubs_eff->Draw("Psame");
          
          low_det = 100*(myEff->ClopperPearson(eff_vals_det[0][0],eff_vals_det[0][1],0.68,0)-eff_vals_det[0][1]/eff_vals_det[0][0]);
          up_det  = 100*(myEff->ClopperPearson(eff_vals_det[0][0],eff_vals_det[0][1],0.68,1)-eff_vals_det[0][1]/eff_vals_det[0][0]);
          
          TPaveText *txt_pteff_0 = new TPaveText(min_pt+0.1*(max_pt-min_pt),0.26,min_pt+0.7*(max_pt-min_pt),0.34,"br");
          sprintf(buffer_det,"#epsilon_{stubs}=%.2f_{%.2f}^{+%.2f}%%",100*eff_vals_det[0][1]/eff_vals_det[0][0],low_det,up_det);
          txt_pteff_0->AddText(buffer_det);
          txt_pteff_0->SetTextFont(102);
          txt_pteff_0->SetTextSize(0.04);
          txt_pteff_0->Draw();                                    
        }
        
        if (i==1) {
          pt_tower_eff->SetMarkerStyle(20);
          pt_tower_eff->Draw("Psame");
          
          low_det = 100*(myEff->ClopperPearson(eff_vals_det[0][1],eff_vals_det[0][2],0.68,0)-eff_vals_det[0][2]/eff_vals_det[0][1]);
          up_det  = 100*(myEff->ClopperPearson(eff_vals_det[0][1],eff_vals_det[0][2],0.68,1)-eff_vals_det[0][2]/eff_vals_det[0][1]);
          
          TPaveText *txt_pteff_1 = new TPaveText(min_pt+0.1*(max_pt-min_pt),0.26,min_pt+0.7*(max_pt-min_pt),0.34,"br");
          sprintf(buffer_det,"#epsilon_{layers}=%.2f_{%.2f}^{+%.2f}%%",100*eff_vals_det[0][2]/eff_vals_det[0][1],low_det,up_det);
          txt_pteff_1->AddText(buffer_det);
          txt_pteff_1->SetTextFont(102);
          txt_pteff_1->SetTextSize(0.04);
          txt_pteff_1->Draw();
        }
        
        if (i==2) {
          pt_AM_eff->SetMarkerStyle(20);
          pt_AM_eff->Draw("Psame");
          pt_inclusive_eff->SetMarkerStyle(24);
          pt_inclusive_eff->SetMarkerColor(4);
          pt_inclusive_eff->Draw("Psame");
          
          low_det = 100*(myEff->ClopperPearson(eff_vals_det[0][2],eff_vals_det[0][3],0.68,0)-eff_vals_det[0][3]/eff_vals_det[0][2]);
          up_det  = 100*(myEff->ClopperPearson(eff_vals_det[0][2],eff_vals_det[0][3],0.68,1)-eff_vals_det[0][3]/eff_vals_det[0][2]);
          
          TPaveText *txt_pteff_2 = new TPaveText(min_pt+0.1*(max_pt-min_pt),0.26,min_pt+0.7*(max_pt-min_pt),0.34,"br");
          sprintf(buffer_det,"#epsilon_{layers_fe}=%.2f_{%.2f}^{+%.2f}%%",100*eff_vals_det[0][3]/eff_vals_det[0][2],low_det,up_det);
          txt_pteff_2->AddText(buffer_det);
          txt_pteff_2->SetTextFont(102);
          txt_pteff_2->SetTextSize(0.04);
          txt_pteff_2->Draw();
     
          low_det = 100*(myEff->ClopperPearson(eff_vals_det[0][0],eff_vals_det[0][4],0.68,0)-eff_vals_det[0][4]/eff_vals_det[0][0]);
          up_det  = 100*(myEff->ClopperPearson(eff_vals_det[0][0],eff_vals_det[0][4],0.68,1)-eff_vals_det[0][4]/eff_vals_det[0][0]);
          
          TPaveText *txt_pteff = new TPaveText(min_pt+0.1*(max_pt-min_pt),0.15,min_pt+0.9*(max_pt-min_pt),0.23,"br");
          sprintf(buffer_det,"#epsilon_{incl}^{N#geq%d,p_{T}#geq%.1f}=%.2f_{%.2f}^{+%.2f}%%",nhits,ptmin,100*eff_vals_det[0][3]/eff_vals_det[0][0],low_det,up_det);
          txt_pteff->AddText(buffer_det);
          txt_pteff->SetFillColor(1);
          txt_pteff->SetTextColor(0);
          txt_pteff->SetTextFont(112);
          txt_pteff->SetTextSize(0.05);
          txt_pteff->Draw();
        }    
      }
      c_pteff->Update();
      c_pteff->Write();
      c_pteff->SaveAs(plotDIR+plotPRE+outputname+"_eff_pt.png");
    

    // -----------------------------------------------------------------------------------------------------------
    // Canvas 2: c_etaeff
    TCanvas *c_etaeff = new TCanvas("c_etaeff","eta efficiencies",201,77,1170,608);
      c_etaeff->Divide(3,1);
      
      for (int i=0;i<3;++i) {
        c_etaeff->cd(i+1);
        c_etaeff->cd(i+1)->Range(-5.887851,-1.930603,70.65421,17.37543);
        c_etaeff->cd(i+1)->SetFillColor(0);
        c_etaeff->cd(i+1)->SetBorderMode(0);
        c_etaeff->cd(i+1)->SetBorderSize(2);
        c_etaeff->cd(i+1)->SetGridx();
        c_etaeff->cd(i+1)->SetGridy();
        c_etaeff->cd(i+1)->SetLeftMargin(0.13); //0.07692308
        c_etaeff->cd(i+1)->SetTopMargin(0.065); //0.07124352
        c_etaeff->cd(i+1)->SetFrameBorderMode(0);
        h_eff_eta_det->GetXaxis()->SetTitle("Tracking particle #eta");
        h_eff_eta_det->GetYaxis()->SetTitle("Efficiency");
        h_eff_eta_det->Draw();
        
        if (i==0) {
          eta_stubs_eff->SetMarkerStyle(20);
          eta_stubs_eff->Draw("Psame");
          
          low_det = 100*(myEff->ClopperPearson(eff_vals_det[0][0],eff_vals_det[0][1],0.68,0)-eff_vals_det[0][1]/eff_vals_det[0][0]);
          up_det  = 100*(myEff->ClopperPearson(eff_vals_det[0][0],eff_vals_det[0][1],0.68,1)-eff_vals_det[0][1]/eff_vals_det[0][0]);
          
          TPaveText *txt_etaeff_0 = new TPaveText(min_eta+0.1*(max_eta-min_eta),0.26,min_eta+0.7*(max_eta-min_eta),0.34,"br");
          sprintf(buffer_det,"#epsilon_{stubs}=%.2f_{%.2f}^{+%.2f}%%",100*eff_vals_det[0][1]/eff_vals_det[0][0],low_det,up_det);
          txt_etaeff_0->AddText(buffer_det);
          txt_etaeff_0->SetTextFont(102);
          txt_etaeff_0->SetTextSize(0.04);
          txt_etaeff_0->Draw(); 
        }
        
        if (i==1) {
          eta_tower_eff->SetMarkerStyle(20);
          eta_tower_eff->Draw("Psame");

          low_det = 100*(myEff->ClopperPearson(eff_vals_det[0][1],eff_vals_det[0][2],0.68,0)-eff_vals_det[0][2]/eff_vals_det[0][1]);
          up_det  = 100*(myEff->ClopperPearson(eff_vals_det[0][1],eff_vals_det[0][2],0.68,1)-eff_vals_det[0][2]/eff_vals_det[0][1]);
          
          TPaveText *txt_etaeff_1 = new TPaveText(min_eta+0.1*(max_eta-min_eta),0.26,min_eta+0.7*(max_eta-min_eta),0.34,"br");
          sprintf(buffer_det,"#epsilon_{layers}=%.2f_{%.2f}^{+%.2f}%%",100*eff_vals_det[0][2]/eff_vals_det[0][1],low_det,up_det);
          txt_etaeff_1->AddText(buffer_det);
          txt_etaeff_1->SetTextFont(102);
          txt_etaeff_1->SetTextSize(0.04);
          txt_etaeff_1->Draw();
        }
        
        if (i==2) {
          eta_AM_eff->SetMarkerStyle(20);
          eta_AM_eff->Draw("Psame");
          eta_inclusive_eff->SetMarkerStyle(24);
          eta_inclusive_eff->SetMarkerColor(4);
          eta_inclusive_eff->Draw("Psame");

          low_det = 100*(myEff->ClopperPearson(eff_vals_det[0][2],eff_vals_det[0][3],0.68,0)-eff_vals_det[0][3]/eff_vals_det[0][2]);
          up_det  = 100*(myEff->ClopperPearson(eff_vals_det[0][2],eff_vals_det[0][3],0.68,1)-eff_vals_det[0][3]/eff_vals_det[0][2]);
          
          TPaveText *txt_etaeff_2 = new TPaveText(min_eta+0.1*(max_eta-min_eta),0.26,min_eta+0.7*(max_eta-min_eta),0.34,"br");
          sprintf(buffer_det,"#epsilon_{layers_fe}=%.2f_{%.2f}^{+%.2f}%%",100*eff_vals_det[0][3]/eff_vals_det[0][2],low_det,up_det);
          txt_etaeff_2->AddText(buffer_det);
          txt_etaeff_2->SetTextFont(102);
          txt_etaeff_2->SetTextSize(0.04);
          txt_etaeff_2->Draw();
     
          low_det = 100*(myEff->ClopperPearson(eff_vals_det[0][0],eff_vals_det[0][4],0.68,0)-eff_vals_det[0][4]/eff_vals_det[0][0]);
          up_det  = 100*(myEff->ClopperPearson(eff_vals_det[0][0],eff_vals_det[0][4],0.68,1)-eff_vals_det[0][4]/eff_vals_det[0][0]);
          
          TPaveText *txt_etaeff = new TPaveText(min_eta+0.1*(max_eta-min_eta),0.15,min_eta+0.9*(max_eta-min_eta),0.23,"br");
          sprintf(buffer_det,"#epsilon_{incl}^{N#geq%d,p_{T}#geq%.1f}=%.2f_{%.2f}^{+%.2f}%%",nhits,ptmin,100*eff_vals_det[0][3]/eff_vals_det[0][0],low_det,up_det);
          txt_etaeff->AddText(buffer_det);
          txt_etaeff->SetFillColor(1);
          txt_etaeff->SetTextColor(0);
          txt_etaeff->SetTextFont(112);
          txt_etaeff->SetTextSize(0.05);
          txt_etaeff->Draw();
        }
      }
      c_etaeff->Update();
      c_etaeff->Write();
      c_etaeff->SaveAs(plotDIR+plotPRE+outputname+"_eff_eta.png");
    
    
    // -----------------------------------------------------------------------------------------------------------
    // Canvas 3: c_zeff
    TCanvas *c_zeff = new TCanvas("c_zeff","z efficiencies",201,77,1170,608);
      c_zeff->Divide(3,1);
      
      for (int i=0;i<3;++i) {
        c_zeff->cd(i+1);
        c_zeff->cd(i+1)->Range(-5.887851,-1.930603,70.65421,17.37543);
        c_zeff->cd(i+1)->SetFillColor(0);
        c_zeff->cd(i+1)->SetBorderMode(0);
        c_zeff->cd(i+1)->SetBorderSize(2);
        c_zeff->cd(i+1)->SetGridx();
        c_zeff->cd(i+1)->SetGridy();
        c_zeff->cd(i+1)->SetLeftMargin(0.13); //0.07692308
        c_zeff->cd(i+1)->SetTopMargin(0.065); //0.07124352
        c_zeff->cd(i+1)->SetFrameBorderMode(0);
        h_eff_z_det->GetXaxis()->SetTitle("Tracking particle #z");
        h_eff_z_det->GetYaxis()->SetTitle("Efficiency");
        h_eff_z_det->Draw();
        
        if (i==0) {
          z_stubs_eff->SetMarkerStyle(20);
          z_stubs_eff->Draw("Psame");
        }
        
        if (i==1) {
          z_tower_eff->SetMarkerStyle(20);
          z_tower_eff->Draw("Psame");
        }
        
        if (i==2) {
          z_AM_eff->SetMarkerStyle(20);
          z_AM_eff->Draw("Psame");
          z_inclusive_eff->SetMarkerStyle(24);
          z_inclusive_eff->SetMarkerColor(4);
          z_inclusive_eff->Draw("Psame");
        }
      }
      c_zeff->Update(); 
      c_zeff->Write();
      c_zeff->SaveAs(plotDIR+plotPRE+outputname+"_eff_z.png");
    
    
    // -----------------------------------------------------------------------------------------------------------
    // Canvas 4: c_d0eff
    TCanvas *c_d0eff = new TCanvas("c4","d0 efficiencies",201,77,1170,608);
      c_d0eff->Divide(3,1);
      
      for (int i=0;i<3;++i) {
        c_d0eff->cd(i+1);
        c_d0eff->cd(i+1)->Range(-5.887851,-1.930603,70.65421,17.37543);
        c_d0eff->cd(i+1)->SetFillColor(0);
        c_d0eff->cd(i+1)->SetBorderMode(0);
        c_d0eff->cd(i+1)->SetBorderSize(2);
        c_d0eff->cd(i+1)->SetGridx();
        c_d0eff->cd(i+1)->SetGridy();
        c_d0eff->cd(i+1)->SetLeftMargin(0.13); //0.07692308
        c_d0eff->cd(i+1)->SetTopMargin(0.065); //0.07124352
        c_d0eff->cd(i+1)->SetFrameBorderMode(0);
        h_eff_d_det->GetXaxis()->SetTitle("Tracking particle d_{0}(in cm)");
        h_eff_d_det->GetYaxis()->SetTitle("Efficiency");
        h_eff_d_det->Draw();
        
        if (i==0) {
          d_stubs_eff->SetMarkerStyle(20);
          d_stubs_eff->Draw("Psame");
        }
        
        if (i==1) {
          d_tower_eff->SetMarkerStyle(20);
          d_tower_eff->Draw("Psame");
        }
        
        if (i==2) {
          d_AM_eff->SetMarkerStyle(20);
          d_AM_eff->Draw("Psame");
          d_inclusive_eff->SetMarkerStyle(24);
          d_inclusive_eff->SetMarkerColor(4);
          d_inclusive_eff->Draw("Psame");
        }
      }
      c_d0eff->Update();
      c_d0eff->Write();
      c_d0eff->SaveAs(plotDIR+plotPRE+outputname+"_eff_d0.png");
    


  // -----------------------------------------------------------------------------------------------------------
  // Acceptance Study
  // -----------------------------------------------------------------------------------------------------------
    const int nbins_std = 25;
      
    float min_pt_std  = 0;
    float max_pt_std  = 10;

    float min_z0_std  = -15.;
    float max_z0_std  = 15.;

    float max_eta_std = 2.4;

    float bin_pt_std  = (max_pt_std-min_pt_std)/nbins_std;
    float bin_eta_std = max_eta_std/nbins_std;
      
    float sig_val_std[5][20][nbins_std];
    float count_std[5][20][nbins_std];  
    float abscissa_std[5][2][nbins_std];  

    float eff_val_std[nbins_std][nbins_std][4][3][4];
    float eff_vals_std[4][3][4];
    
    for (int i=0;i<nbins_std;++i) {
      for (int j=0;j<nbins_std;++j) {
        for (int k=0;k<4;++k) {
          for (int l=0;l<3;++l) {
            for (int m=0;m<4;++m) {
              eff_vals_std[k][l][m]      = 0.;
              eff_val_std[i][j][k][l][m] = 0.;
            } // close m loop
          } // close l loop
        } // close k loop
      } // close j loop
    } // close i loop


    // -----------------------------------------------------------------------------------------------------------
    // Definition of eff_val_std values
    
      // eff_val[ptbin_std][etabin_std][type][nhits][j]
      // type = 0 all primaries (rho<1cm)
      // type = 1 only mu (rho<1cm)
      // type = 2 only e  (rho<1cm)
      // type = 3 only pi (rho<1cm)
      //
      // nhits = 0 <-> 4
      // nhits = 1 <-> 5
      // nhits = 2 <-> 6
      //
      // Then
      //
      // eff_val_std[][][][][0] : all the particles in the acceptance with at least 1 stub
      // eff_val_std[][][][][1] : all the particles in [0] with at least nhits stubs
      // eff_val_std[][][][][2] : all the particles in [1] with stubs in at least nhits layer/disks
      // eff_val_std[][][][][3] : all the particles in [2] with stubs passing FE in at least nhits layer/disks

    
    
    // -----------------------------------------------------------------------------------------------------------
    // Loop over Events and Particles
    int bin_std[5];
    
    float npart_std_0=0;
    float npart_std_1=0;
    float npart_std_2=0;
    float npart_std_3=0;
    float npart_std_4=0;
    
    // Event Loop
    for (int j=0;j<std::min(n_entries,nevt);++j) {
      if (j%1000==0) cout << j << endl;
      
      int n_comb_tot_std = 0;
        
      FullI->GetEntry(j); // Get the tree info
        
      // First, count the number of patterns in the sector
      int npatt_std=0;
      int npatt_full_std=0;
      int ptbin_std;
      int inacc_std;
      int sectype_std;
      int etabin_std;

      int pdid_std;
      int type,nstubs;
      int ns,nh,nhf;

      // loop over all the particles
      for (int k=0;k<n_part;++k) {
        if (part_rho->at(k)>1) continue;
        if (part_nstubs->at(k)<1) continue;
        if (part_pt->at(k)<min_pt_std) continue;
        if (part_pt->at(k)>max_pt_std) continue;
        if (part_z0->at(k)>max_z0_std) continue;
        if (part_z0->at(k)<min_z0_std) continue;
        if (std::abs(part_eta->at(k))>max_eta_std) continue;

        ptbin_std  = std::min(int((part_pt->at(k)-min_pt_std)/bin_pt_std),nbins_std-1);
        etabin_std = std::min(int((std::abs(part_eta->at(k)))/bin_eta_std),nbins_std-1);
        pdid_std   = std::abs(part_pdg->at(k));
        ns     = std::min(part_nstubs->at(k),6);
        nh     = std::min(part_nhits->at(k),6);
        nhf    = std::min(part_nhits_fe->at(k),6);
        
        type=0;
        if (pdid_std==13)  type=1;
        if (pdid_std==11)  type=2;
        if (pdid_std==211) type=3;
      
        for (int l=0;l<3;++l) {
          ++eff_val_std[ptbin_std][etabin_std][0][l][0];
          if (type!=0) ++eff_val_std[ptbin_std][etabin_std][type][l][0];

          if (part_pt->at(k)<3) continue;
          ++eff_vals_std[0][l][0];
          if (type!=0) ++eff_vals_std[type][l][0];
        }
      
        if (ns<4) continue;
                
        for (int l=0;l<ns-3;++l) {
          ++eff_val_std[ptbin_std][etabin_std][0][l][1];
          if (type!=0) ++eff_val_std[ptbin_std][etabin_std][type][l][1];

          if (part_pt->at(k)<3) continue;
          ++eff_vals_std[0][l][1];
          if (type!=0) ++eff_vals_std[type][l][1];
        }
        
        if (nh<4) continue;            

        for (int l=0;l<nh-3;++l) {
          ++eff_val_std[ptbin_std][etabin_std][0][l][2];
          if (type!=0) ++eff_val_std[ptbin_std][etabin_std][type][l][2];

          if (part_pt->at(k)<3) continue;
          ++eff_vals_std[0][l][2];
          if (type!=0) ++eff_vals_std[type][l][2];
        }
      
        if (nhf<4) continue;            
        
        for (int l=0;l<nhf-3;++l) {
          ++eff_val_std[ptbin_std][etabin_std][0][l][3];
          if (type!=0) ++eff_val_std[ptbin_std][etabin_std][type][l][3];

          if (part_pt->at(k)<3) continue;
          ++eff_vals_std[0][l][3];
          if (type!=0) ++eff_vals_std[type][l][3];
        }
      } // End of loop over particles
    } // End of loop on events


    // -----------------------------------------------------------------------------------------------------------
    // Plots
    std::vector<TH2F*> h_eff_all; 
    std::vector<TH2F*> h_eff_mu; 
    std::vector<TH2F*> h_eff_ele; 
    std::vector<TH2F*> h_eff_pi; 

    for (int i=0;i<3;++i) {
      TH2F *h_eff_all_0  = new TH2F("eff_a0","eff_a0",nbins_std,0,max_eta_std,nbins_std,min_pt_std,max_pt_std);
      TH2F *h_eff_all_1  = new TH2F("eff_a1","eff_a1",nbins_std,0,max_eta_std,nbins_std,min_pt_std,max_pt_std);
      TH2F *h_eff_all_2  = new TH2F("eff_a2","eff_a2",nbins_std,0,max_eta_std,nbins_std,min_pt_std,max_pt_std);
      TH2F *h_eff_all_3  = new TH2F("eff_a3","eff_a3",nbins_std,0,max_eta_std,nbins_std,min_pt_std,max_pt_std);

      h_eff_all.push_back(h_eff_all_0);
      h_eff_mu.push_back(h_eff_all_1);
      h_eff_ele.push_back(h_eff_all_2);
      h_eff_pi.push_back(h_eff_all_3);
    }
      

    for (int j=0;j<nbins_std;++j) {
      for (int k=0;k<nbins_std;++k) {
        for (int kk=0;kk<4;++kk) {
          for (int l=0;l<3;++l) {
            
            if (eff_val_std[j][k][kk][l][0]==0) continue;

            for (int m=1;m<4;++m) eff_val_std[j][k][kk][l][m] /= eff_val_std[j][k][kk][l][0];
      
            if (kk==0) h_eff_all[l]->Fill((k+0.5)*bin_eta_std,min_pt_std+(j+0.5)*bin_pt_std,eff_val_std[j][k][kk][l][2]);
            if (kk==1) h_eff_mu[l]->Fill((k+0.5)*bin_eta_std,min_pt_std+(j+0.5)*bin_pt_std,eff_val_std[j][k][kk][l][2]);
            if (kk==2) h_eff_ele[l]->Fill((k+0.5)*bin_eta_std,min_pt_std+(j+0.5)*bin_pt_std,eff_val_std[j][k][kk][l][2]);
            if (kk==3) h_eff_pi[l]->Fill((k+0.5)*bin_eta_std,min_pt_std+(j+0.5)*bin_pt_std,eff_val_std[j][k][kk][l][2]);
          } // close l loop
        } // close kk loop
      } // close k loop
    } // close j loop

    /*
    cout << eff_vals_std[0][0] << " particles pass the kinematic cuts and have at least one stub" << endl;
    cout << eff_vals_std[0][1] << " of them have at least " << nhits << " stubs" << endl;
    cout << eff_vals_std[0][2] << " of them have stubs in at least " << nhits << " layer disks" << endl;
    cout << eff_vals_std[0][3] << " of them have stubs passing FE cuts in at least " << nhits << " layer disks" << endl;
    */


    // -----------------------------------------------------------------------------------------------------------
    // Canvases
    char buffer_std[80];
    float low_std,up_std;
    
    // -----------------------------------------------------------------------------------------------------------
    // Canvas 1: c_trkeff
    TCanvas *c_trkeff = new TCanvas("c_trkeff","All tracks efficiencies",201,77,870,808);
      c_trkeff->Divide(1,3);
        
      for (int i=0;i<3;++i) {
        c_trkeff->cd(i+1);
        c_trkeff->cd(i+1)->SetFillColor(0);
        c_trkeff->cd(i+1)->SetGridx();
        c_trkeff->cd(i+1)->SetGridy();
        h_eff_all[i]->GetXaxis()->SetLabelSize(0.05);
        h_eff_all[i]->GetXaxis()->SetTitleSize(0.05);
        h_eff_all[i]->GetYaxis()->SetLabelSize(0.05);
        h_eff_all[i]->GetYaxis()->SetTitleSize(0.05);
        h_eff_all[i]->GetYaxis()->SetTitleOffset(0.5);
        h_eff_all[i]->SetAxisRange(0., 1.0,"Z");
        h_eff_all[i]->GetYaxis()->SetTitle("Particle p_{T} (GeV)");
        h_eff_all[i]->GetXaxis()->SetTitle("Particle #eta");
        h_eff_all[i]->Draw("cont4z");
          
        low_std = 100*(myEff->ClopperPearson(eff_vals_std[0][i][0],eff_vals_std[0][i][2],0.68,0)-eff_vals_std[0][i][2]/eff_vals_std[0][i][0]);
        up_std  = 100*(myEff->ClopperPearson(eff_vals_std[0][i][0],eff_vals_std[0][i][2],0.68,1)-eff_vals_std[0][i][2]/eff_vals_std[0][i][0]);
        
        TPaveText *txt_trkeff = new TPaveText(0.25,0.25,0.5,0.5,"br");
        sprintf(buffer_std,"#epsilon_{p_{T}>3}^{N_{h}>%d}=%.2f_{%.2f}^{+%.2f}%%",i+4,100*eff_vals_std[0][i][2]/eff_vals_std[0][i][0],low_std,up_std);
        txt_trkeff->AddText(buffer_std);
        txt_trkeff->SetTextFont(102);
        txt_trkeff->SetTextSize(0.055);
        txt_trkeff->SetFillColor(0);
        txt_trkeff->Draw(); 
      }
      c_trkeff->Update();
      c_trkeff->Write();
      c_trkeff->SaveAs(plotDIR+plotPRE+outputname+"_eff_track.png");


    // -----------------------------------------------------------------------------------------------------------
    // Canvas 2: c_mueff
    TCanvas *c_mueff = new TCanvas("c_mueff","Muons efficiencies",201,77,870,808);
      c_mueff->Divide(1,3);
        
      for (int i=0;i<3;++i) {
        c_mueff->cd(i+1);
        c_mueff->cd(i+1)->SetFillColor(0);
        c_mueff->cd(i+1)->SetGridx();
        c_mueff->cd(i+1)->SetGridy();
        h_eff_mu[i]->GetXaxis()->SetLabelSize(0.05);
        h_eff_mu[i]->GetXaxis()->SetTitleSize(0.05);
        h_eff_mu[i]->GetYaxis()->SetLabelSize(0.05);
        h_eff_mu[i]->GetYaxis()->SetTitleSize(0.05);
        h_eff_mu[i]->GetYaxis()->SetTitleOffset(0.5);
        h_eff_mu[i]->SetAxisRange(0., 1.0,"Z");
        h_eff_mu[i]->GetYaxis()->SetTitle("Particle p_{T} (GeV)");
        h_eff_mu[i]->GetXaxis()->SetTitle("Particle #eta");
        h_eff_mu[i]->Draw("cont4z");
        
        low_std = 100*(myEff->ClopperPearson(eff_vals_std[1][i][0],eff_vals_std[1][i][2],0.68,0)-eff_vals_std[1][i][2]/eff_vals_std[1][i][0]);
        up_std  = 100*(myEff->ClopperPearson(eff_vals_std[1][i][0],eff_vals_std[1][i][2],0.68,1)-eff_vals_std[1][i][2]/eff_vals_std[1][i][0]);
        
        TPaveText *txt_mueff = new TPaveText(0.25,0.25,0.5,0.5,"br");
        sprintf(buffer_std,"#epsilon_{p_{T}>3}^{N_{h}>%d}=%.2f_{%.2f}^{+%.2f}%%",i+4,100*eff_vals_std[1][i][2]/eff_vals_std[1][i][0],low_std,up_std);
        txt_mueff->AddText(buffer_std);
        txt_mueff->SetTextFont(102);
        txt_mueff->SetTextSize(0.055);
        txt_mueff->SetFillColor(0);
        txt_mueff->Draw(); 
      }
      c_mueff->Update();
      c_mueff->Write();
      c_mueff->SaveAs(plotDIR+plotPRE+outputname+"_eff_muon.png");


    // -----------------------------------------------------------------------------------------------------------
    // Canvas 3: c_eleff
    TCanvas *c_eleff = new TCanvas("c_eleff","Electrons efficiencies",201,77,870,808);
      c_eleff->Divide(1,3);
        
      for (int i=0;i<3;++i) {
        c_eleff->cd(i+1);
        c_eleff->cd(i+1)->SetFillColor(0);
        c_eleff->cd(i+1)->SetGridx();
        c_eleff->cd(i+1)->SetGridy();
        h_eff_ele[i]->GetXaxis()->SetLabelSize(0.05);
        h_eff_ele[i]->GetXaxis()->SetTitleSize(0.05);
        h_eff_ele[i]->GetYaxis()->SetLabelSize(0.05);
        h_eff_ele[i]->GetYaxis()->SetTitleSize(0.05);
        h_eff_ele[i]->GetYaxis()->SetTitleOffset(0.5);
        h_eff_ele[i]->SetAxisRange(0., 1.0,"Z");
        h_eff_ele[i]->GetYaxis()->SetTitle("Particle p_{T} (GeV)");
        h_eff_ele[i]->GetXaxis()->SetTitle("Particle #eta");
        h_eff_ele[i]->Draw("cont4z");
        
        low_std = 100*(myEff->ClopperPearson(eff_vals_std[2][i][0],eff_vals_std[2][i][2],0.68,0)-eff_vals_std[2][i][2]/eff_vals_std[2][i][0]);
        up_std  = 100*(myEff->ClopperPearson(eff_vals_std[2][i][0],eff_vals_std[2][i][2],0.68,1)-eff_vals_std[2][i][2]/eff_vals_std[2][i][0]);
        
        TPaveText *txt_eleff = new TPaveText(0.25,0.25,0.5,0.5,"br");
        sprintf(buffer_std,"#epsilon_{p_{T}>3}^{N_{h}>%d}=%.2f_{%.2f}^{+%.2f}%%",i+4,100*eff_vals_std[2][i][2]/eff_vals_std[2][i][0],low_std,up_std);
        txt_eleff->AddText(buffer_std);
        txt_eleff->SetTextFont(102);
        txt_eleff->SetTextSize(0.055);
        txt_eleff->SetFillColor(0);
        txt_eleff->Draw(); 
      }
      c_eleff->Update();
      c_eleff->Write();
      c_eleff->SaveAs(plotDIR+plotPRE+outputname+"_eff_electron.png");
        

    // -----------------------------------------------------------------------------------------------------------
    // Canvas 4: c_pieff
    TCanvas *c_pieff = new TCanvas("c_pieff","Pions efficiencies",201,77,870,808);
      c_pieff->Divide(1,3);  
      
      for (int i=0;i<3;++i) {
        c_pieff->cd(i+1);
        c_pieff->cd(i+1)->SetFillColor(0);
        c_pieff->cd(i+1)->SetGridx();
        c_pieff->cd(i+1)->SetGridy();
        h_eff_pi[i]->GetXaxis()->SetLabelSize(0.05);
        h_eff_pi[i]->GetXaxis()->SetTitleSize(0.05);
        h_eff_pi[i]->GetYaxis()->SetLabelSize(0.05);
        h_eff_pi[i]->GetYaxis()->SetTitleSize(0.05);
        h_eff_pi[i]->GetYaxis()->SetTitleOffset(0.5);
        h_eff_pi[i]->SetAxisRange(0., 1.0,"Z");
        h_eff_pi[i]->GetYaxis()->SetTitle("Particle p_{T} (GeV)");
        h_eff_pi[i]->GetXaxis()->SetTitle("Particle #eta");
        h_eff_pi[i]->Draw("cont4z");
        
        low_std = 100*(myEff->ClopperPearson(eff_vals_std[3][i][0],eff_vals_std[3][i][2],0.68,0)-eff_vals_std[3][i][2]/eff_vals_std[3][i][0]);
        up_std  = 100*(myEff->ClopperPearson(eff_vals_std[3][i][0],eff_vals_std[3][i][2],0.68,1)-eff_vals_std[3][i][2]/eff_vals_std[3][i][0]);
        
        TPaveText *txt_pieff = new TPaveText(0.25,0.25,0.5,0.5,"br");
        sprintf(buffer_std,"#epsilon_{p_{T}>3}^{N_{h}>%d}=%.2f_{%.2f}^{+%.2f}%%",i+4,100*eff_vals_std[3][i][2]/eff_vals_std[3][i][0],low_std,up_std);
        txt_pieff->AddText(buffer_std);
        txt_pieff->SetTextFont(102);
        txt_pieff->SetTextSize(0.055);
        txt_pieff->SetFillColor(0);
        txt_pieff->Draw(); 
      }
      c_pieff->Update();
      c_pieff->Write();
      c_pieff->SaveAs(plotDIR+plotPRE+outputname+"_eff_pion.png");



} // end void DetectorAcceptance()


// ---------------------------------------------------------------------------------------------------------
// Plot Style
void SetPlotStyle() 
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  gStyle->SetNumberContours(90);
  gStyle->SetPalette(kBlueYellow);

  // gStyle->SetCanvasColor(0);
  // gStyle->SetCanvasBorderMode(0);
  // gStyle->SetCanvasBorderSize(2);

  // gStyle->SetFrameBorderMode(0);

  // gStyle->SetLabelFont(42,"xyz");
  // gStyle->SetTitleFont(42,"xyz");

  // gStyle->SetLabelSize(0.035,"xyz");
  // gStyle->SetTitleSize(0.035,"xyz");

  // gStyle->SetLabelOffset(0.004,"y");

  // gStyle->SetPadLeftMargin(0.07692308);
  // gStyle->SetPadTopMargin(0.07124352);
}
