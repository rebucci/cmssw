/*
  ROOT macro for stub efficiencies visualization

Use:
    root[0]-> .L PlotStubEffs.C
    root[1]->PlotStubEffs("sourcefile","outputname",pu)
where 
    - sourcefile is the name of the root output file produced by AM_ana. It does not include the ".root" suffix, but it does have a preceding RootFiles/ path. This can be changed in the Preliminaries section of this macro. sourcefile must be enclosed in quotes.
    - outputname is the name (no type suffix like ".root") of the root file this macro produces. it is also the prefix of all the plot names.
    - pu is the pileup value of the source file. It is used in the title only.


Options
    - ptmax: maximum pT for x axis on plot. Default is 20 GeV/c.
    - doLayerEff: efficiencies calculated for stubs in layer/disk (true) OR in modules (false)
    - PrintEachLayer: save .png files for every layer.
    - SavePlotSource: save the .C macro files for the pT barrel and endcap summaries

Module vs Layer/Disk Efficiency
    - Module efficiency: probability that the particle induces a stub in the module. Much more sensitive to acceptance effects than Layer/Disk efficiency
    - Layer/Disk efficiency: probability that the particle induces a stub in the layer/disk to which belongs the module  

For info about the efficiency definition, have a look at the following presentation:
  https://indico.cern.ch/getFile.py/access?contribId=0&resId=0&materialId=slides&confId=263068
*/


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
#include <iomanip>

void SetPlotStyle();

void PlotStubEffs(TString sourcefile, TString outputname, int pu,  float ptmax=20)
{

  // -----------------------------------------------------------------------------------------------------------
  // Preliminaries
  // -----------------------------------------------------------------------------------------------------------

  gROOT->SetBatch();
  gErrorIgnoreLevel = kError;
  SetPlotStyle();
 
  TString srcDIR  = "RootFiles/"; 
  TString subDIR  = "RootFiles/PlotFiles/";
  TString plotDIR = "Plots/"; 
  TString plotPRE = "stubEff_";

  // Options
  bool doLayerEff     = true;   // module efficiency: 0, layer efficiency: 1 in seb's original
  bool PrintEachLayer = true;  // saves each layer plot
  bool doEtaPlots     = true;  // eta efficiencies
  bool SavePlotSource = false;  // save the macro.C files for the pT barrel and endcap summaries


  // Read in File
  TChain *newtree = new TChain("Efficiencies");
  newtree->Add(srcDIR+sourcefile+".root"); 

  // Output file
  TFile* fout = new TFile(subDIR+"plots_"+plotPRE+outputname+".root","recreate");


  // -----------------------------------------------------------------------------------------------------------
  // Definitions
  // -----------------------------------------------------------------------------------------------------------

  // -----------------------------------------------------------------------------------------------------------
  // Variables

  float pt_val[100];
  float digi_pt[30][100];
  float clus_pt[30][100];
  float stub_pt[30][100];

  float eta_val[50];
  float digi_eta[30][50];    
  float clus_eta[30][50];
  float stub_eta[30][50];


  // -----------------------------------------------------------------------------------------------------------
  // Branch Addresses 

  newtree->SetBranchAddress("pt_val",&pt_val);
  newtree->SetBranchAddress("digi_pt",&digi_pt);
  newtree->SetBranchAddress("clus_pt",&clus_pt);
  if (doLayerEff) newtree->SetBranchAddress("stub_pt_lay",&stub_pt);
  else            newtree->SetBranchAddress("stub_pt",    &stub_pt);


  newtree->SetBranchAddress("eta_val",&eta_val);
  newtree->SetBranchAddress("digi_eta",&digi_eta);
  newtree->SetBranchAddress("clus_eta",&clus_eta);
  if (doLayerEff) newtree->SetBranchAddress("stub_eta_lay",&stub_eta);
  else            newtree->SetBranchAddress("stub_eta",    &stub_eta); 

  // initialize
  newtree->GetEntry(0);


  // -----------------------------------------------------------------------------------------------------------
  // Text in Plots

  char txt_pu[80];
  sprintf (txt_pu, "#sqrt{s}=14TeV, PU %d", pu);

  TLatex TL_PU;
  TL_PU.SetTextSize(0.03);
  TL_PU.SetTextFont(52);

  TLatex TLtitle03;
  TLtitle03.SetTextSize(0.03);

  TLatex TLtitle04;
  TLtitle04.SetTextSize(0.04);


  // ---------------------------------------------------------------------------------------------------------
  //        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
  // ---------------------------------------------------------------------------------------------------------


  // -----------------------------------------------------------------------------------------------------------
  // pT Plots
  // ---------------------------------------------------------------------------------------------------------

  // Loop over Layers
  for (int layer=5;layer<25;++layer)
  {
    int disk =0;
    if (layer>=11  && layer<18) disk= int(layer-11)%8;
    if (layer>=18  && layer<25) disk=-int(layer-18)%8;

    // Histograms
    TH2F *h_eff_pt      = new TH2F("eff_pt",";Particle p_{T} (GeV); Efficiency",300,0.,ptmax,200,0.,1.02);
    TH2F *h_eff_pt_digi = new TH2F("eff_pt_digi",";Particle p_{T} (GeV); Efficiency",300,0.,20.,200,0.,1.02);
    TH2F *h_eff_pt_coff = new TH2F("eff_pt_coff",";Particle p_{T} (GeV); Efficiency",300,0.,20.,200,0.,1.02);
    TH2F *h_eff_pt_soff = new TH2F("eff_pt_soff",";Particle p_{T} (GeV); Efficiency",300,0.,20.,200,0.,1.02);

    // Fill Histograms
    for (int i=0;i<100;++i)
    { 
      h_eff_pt_digi->Fill(pt_val[i],digi_pt[layer][i]);
      h_eff_pt_coff->Fill(pt_val[i],clus_pt[layer][i]);
      h_eff_pt_soff->Fill(pt_val[i],stub_pt[layer][i]);
    }

    TCanvas *c_eff_pt = new TCanvas("c_eff_pt","pT Efficiency",0,0,640,640); //5,75,670,660
      c_eff_pt->Range(-1.2,-0.1,10.6,1.24); //-1.2,-0.1,10.5,1.1
      c_eff_pt->SetGridx();
      c_eff_pt->SetGridy();
      c_eff_pt->SetLeftMargin(0.11); //0.04
      c_eff_pt->SetRightMargin(0.065); //0.05

      h_eff_pt->GetXaxis()->SetLabelSize(0.03);
      h_eff_pt->GetYaxis()->SetLabelSize(0.03);
      h_eff_pt->GetYaxis()->SetTitleOffset(1.25);
      h_eff_pt->GetXaxis()->SetNdivisions(518);
      h_eff_pt->GetYaxis()->SetNdivisions(511);
      h_eff_pt->Draw();
      
      h_eff_pt_digi->SetMarkerStyle(3);
      h_eff_pt_coff->SetMarkerStyle(4);
      h_eff_pt_soff->SetMarkerStyle(20);
      h_eff_pt_digi->Draw("same");
      h_eff_pt_coff->Draw("same");
      h_eff_pt_soff->Draw("same");

      char txt_eff_pt[80];
      if (disk==0) {
        sprintf(txt_eff_pt, "Layer %d",layer);
      }
      else {
        (disk>0)
          ? sprintf(txt_eff_pt, "Disk %d",disk)
          : sprintf(txt_eff_pt, "Disk %d",-disk);
      }

      TLegend *leg_eff_pt = new TLegend(0.7,0.28,0.9,0.45);
        leg_eff_pt->SetTextSize(0.03);
        leg_eff_pt->SetHeader(txt_eff_pt); 
        leg_eff_pt->AddEntry(h_eff_pt_digi,"Digis","p");
        leg_eff_pt->AddEntry(h_eff_pt_coff,"Clusters","p");
        leg_eff_pt->AddEntry(h_eff_pt_soff,"Stubs","p");
        leg_eff_pt->Draw();  

      TLtitle04.DrawLatex(0., 1.03, "CMS Phase-2 Simulation");
      TL_PU.DrawLatex(0.7*ptmax, 1.03, txt_pu);

      c_eff_pt->Modified();
      c_eff_pt->Update();
      c_eff_pt->Write();

      if (PrintEachLayer) {
        char name_eff_pt[100];
        if (doLayerEff) {
          if (disk==0) sprintf (name_eff_pt, "_eff_pt_layer_barrel_%d.png", layer+5);
          if (disk>0)  sprintf (name_eff_pt, "_eff_pt_layer_endcap_p%d.png", disk); 
          if (disk<0)  sprintf (name_eff_pt, "_eff_pt_layer_endcap_m%d.png", abs(disk)); 
        }
        else {
          if (disk==0) sprintf (name_eff_pt, "_eff_pt_module_barrel_%d.png", layer+5);
          if (disk>0)  sprintf (name_eff_pt, "_eff_pt_module_endcap_p%d.png", disk); 
          if (disk<0)  sprintf (name_eff_pt, "_eff_pt_module_endcap_m%d.png", abs(disk)); 
        }
        c_eff_pt->SaveAs(plotDIR+plotPRE+outputname+name_eff_pt);
      }

  }




  // -----------------------------------------------------------------------------------------------------------
  // Eta Plots
  // ---------------------------------------------------------------------------------------------------------
  if (doEtaPlots) {

    // Loop over Layers
    for (int layer=5;layer<18;++layer) // try 25 also?
    {
      int disk =0;
      if (layer>=11  && layer<18) disk= int(layer-11)%8;
      if (layer>=18  && layer<25) disk=-int(layer-18)%8;

      for (int i=0;i<30;++i)
      { 
        for (int j=0;j<50;++j)
        { 
          digi_eta[i][j]  = 0;
          clus_eta[i][j] = 0;
          stub_eta[i][j] = 0;
        }
      }

      TH2F *h_eff_eta      = new TH2F("eff_eta",";Particle #eta; Efficiency",300,-2.5,2.5,200,0.,1.02);
      TH2F *h_eff_eta_digi = new TH2F("eff_eta_digi",";Particle #eta; Efficiency",300,-2.5,2.5,200,0.,1.02);
      TH2F *h_eff_eta_ceff = new TH2F("eff_eta_ceff",";Particle #eta; Efficiency",300,-2.5,2.5,200,0.,1.02);
      TH2F *h_eff_eta_seff = new TH2F("eff_eta_seff",";Particle #eta; Efficiency",300,-2.5,2.5,200,0.,1.02);
        
      for (int i=0;i<50;++i)
      { 
        h_eff_eta_digi->Fill(eta_val[i],digi_eta[layer-5][i]);
        h_eff_eta_ceff->Fill(eta_val[i],clus_eta[layer-5][i]);
        h_eff_eta_seff->Fill(eta_val[i],stub_eta[layer-5][i]);
      }

      TCanvas *c_eff_eta = new TCanvas("c_eff_eta","eta Efficiency",0,0,640,640); //5,75,670,660
        c_eff_eta->Range(-1.2,-0.1,10.6,1.24);
        c_eff_eta->SetGridx();
        c_eff_eta->SetGridy();
        c_eff_eta->SetLeftMargin(0.11); //0.08
        c_eff_eta->SetRightMargin(0.065); //0.05

        h_eff_eta->GetXaxis()->SetLabelSize(0.03);
        h_eff_eta->GetYaxis()->SetLabelSize(0.03);
        h_eff_eta->GetYaxis()->SetTitleOffset(0.8);
        h_eff_eta->GetXaxis()->SetNdivisions(522);
        h_eff_eta->GetYaxis()->SetNdivisions(511);
        h_eff_eta->Draw();
          
        h_eff_eta_digi->SetMarkerStyle(3);
        h_eff_eta_ceff->SetMarkerStyle(4);
        h_eff_eta_seff->SetMarkerStyle(20);
        h_eff_eta_digi->Draw("same");
        h_eff_eta_ceff->Draw("same");
        h_eff_eta_seff->Draw("same");        

        char txt_eff_eta[80];

        if (disk==0)
        {
          sprintf(txt_eff_eta, "Layer %d",layer);
        }
        else
        {
          (disk>0)
            ? sprintf(txt_eff_eta, "Disk %d",disk)
            : sprintf(txt_eff_eta, "Disk %d",-disk);
        }

        TLegend *leg_eff_eta = new TLegend(0.7,0.28,0.9,0.45);
          leg_eff_eta->SetTextSize(0.03);
          leg_eff_eta->SetHeader(txt_eff_eta); 
          leg_eff_eta->AddEntry(h_eff_eta_digi,"Digis","p");
          leg_eff_eta->AddEntry(h_eff_eta_ceff,"Clusters","p");
          leg_eff_eta->AddEntry(h_eff_eta_seff,"Stubs","p");
          leg_eff_eta->Draw();

        TLtitle04.DrawLatex(-2.5, 1.03, "CMS Phase-2 Simulation");
        TL_PU.DrawLatex(0., 1.03, txt_pu);

        c_eff_eta->Modified();
        c_eff_eta->Update();
        c_eff_eta->Write();

        if (PrintEachLayer) {
          char name_eff_eta[100];
          if (doLayerEff) {
            if (disk==0) sprintf (name_eff_eta, "_eff_eta_layer_barrel_%d.png", layer);
            if (disk>0)  sprintf (name_eff_eta, "_eff_eta_layer_endcap_p%d.png", disk); 
            if (disk<0)  sprintf (name_eff_eta, "_eff_eta_layer_endcap_m%d.png", abs(disk)); 
          }
          else {
            if (disk==0) sprintf (name_eff_eta, "_eff_eta_module_barrel_%d.png", layer+5);
            if (disk>0)  sprintf (name_eff_eta, "_eff_eta_module_endcap_p%d.png", disk); 
            if (disk<0)  sprintf (name_eff_eta, "_eff_eta_module_endcap_m%d.png", abs(disk)); 
          }
          c_eff_eta->SaveAs(plotDIR+plotPRE+outputname+name_eff_eta); 
        }

    }

  }


  // -----------------------------------------------------------------------------------------------------------
  // pT Barrel Summary
  // ---------------------------------------------------------------------------------------------------------

  TH2F *h_eff_pt_barrel    = new TH2F("eff_pt_barrel",";Particle p_{T} (GeV); Stub efficiency",300,0.,10.,100,0.,1.15);
  TH2F *h_L1_off   = new TH2F("L1_off","L1_off",300,0.,20.,100,0.,1.02);
  TH2F *h_L2_off   = new TH2F("L2_off","L2_off",300,0.,20.,100,0.,1.02);
  TH2F *h_L3_off   = new TH2F("L3_off","L3_off",300,0.,20.,100,0.,1.02);
  TH2F *h_L4_off   = new TH2F("L4_off","L4_off",300,0.,20.,100,0.,1.02);
  TH2F *h_L5_off   = new TH2F("L5_off","L5_off",300,0.,20.,100,0.,1.02);
  TH2F *h_L6_off   = new TH2F("L6_off","L6_off",300,0.,20.,100,0.,1.02);

  for (int i=0;i<100;++i)
  { 
    if (stub_pt[0][i]!=0) h_L1_off->Fill(pt_val[i],stub_pt[0][i]);
    if (stub_pt[1][i]!=0) h_L2_off->Fill(pt_val[i],stub_pt[1][i]);
    if (stub_pt[2][i]!=0) h_L3_off->Fill(pt_val[i],stub_pt[2][i]);
    if (stub_pt[3][i]!=0) h_L4_off->Fill(pt_val[i],stub_pt[3][i]);
    if (stub_pt[4][i]!=0) h_L5_off->Fill(pt_val[i],stub_pt[4][i]);
    if (stub_pt[5][i]!=0) h_L6_off->Fill(pt_val[i],stub_pt[5][i]);
  }

  TCanvas *c_eff_pt_barrel = new TCanvas("c_eff_pt_barrel","pT Efficiency Barrel",0,0,640,640);
    c_eff_pt_barrel->Range(-1.2,-0.1,10.6,1.24);
    c_eff_pt_barrel->SetGridx();
    c_eff_pt_barrel->SetGridy();
    c_eff_pt_barrel->SetRightMargin(0.07);
    c_eff_pt_barrel->SetTopMargin(0.065);
      
    h_eff_pt_barrel->GetXaxis()->SetTitleSize(0.045);
    h_eff_pt_barrel->GetYaxis()->SetTitleSize(0.035);
    h_eff_pt_barrel->GetYaxis()->SetTitleOffset(0.88);
    h_eff_pt_barrel->GetXaxis()->SetTitleOffset(0.9);
    h_eff_pt_barrel->GetXaxis()->SetNdivisions(518);
    h_eff_pt_barrel->GetYaxis()->SetNdivisions(511);
    h_eff_pt_barrel->Draw();
    
    h_L1_off->SetMarkerStyle(24);
    h_L2_off->SetMarkerStyle(20);
    h_L3_off->SetMarkerStyle(27);
    h_L4_off->SetMarkerStyle(34);
    h_L5_off->SetMarkerStyle(28);
    h_L6_off->SetMarkerStyle(25);
      
    h_L1_off->SetMarkerSize(1.4);
    h_L2_off->SetMarkerSize(1.4);
    h_L3_off->SetMarkerSize(1.4);
    h_L4_off->SetMarkerSize(1.4);
    h_L5_off->SetMarkerSize(1.4);
    h_L6_off->SetMarkerSize(1.4);

    h_L1_off->Draw("same");
    h_L2_off->Draw("same");
    h_L3_off->Draw("same");
    h_L4_off->Draw("same");
    h_L5_off->Draw("same");
    h_L6_off->Draw("same");

    TLegend *leg_eff_pt_bar = new TLegend(0.6,0.2,0.86,0.45);
      leg_eff_pt_bar->SetTextSize(0.035);
      leg_eff_pt_bar->SetFillColor(0);
      leg_eff_pt_bar->AddEntry(h_L1_off,"TBPS layer 1","p");
      leg_eff_pt_bar->AddEntry(h_L2_off,"TBPS layer 2","p");
      leg_eff_pt_bar->AddEntry(h_L3_off,"TBPS layer 3","p");
      leg_eff_pt_bar->AddEntry(h_L4_off,"TB2S layer 1","p");
      leg_eff_pt_bar->AddEntry(h_L5_off,"TB2S layer 2","p");
      leg_eff_pt_bar->AddEntry(h_L6_off,"TB2S layer 3","p");
      leg_eff_pt_bar->Draw();
       
    TLtitle03.DrawLatex(-0.02, 1.18, "CMS Phase-2 simulation");
    TL_PU.DrawLatex(7., 1.18, txt_pu);
    
    c_eff_pt_barrel->Modified();
    c_eff_pt_barrel->Update();
    c_eff_pt_barrel->Write();

    char name_eff_pt_bar[100];
    if (doLayerEff) {
      sprintf (name_eff_pt_bar, "_eff_pt_barrel_layer.png");
      c_eff_pt_barrel->SaveAs(plotDIR+plotPRE+outputname+name_eff_pt_bar);
    }
    else {
      sprintf (name_eff_pt_bar, "_eff_pt_barrel_module.png");
      c_eff_pt_barrel->SaveAs(plotDIR+plotPRE+outputname+name_eff_pt_bar);
    }

    if (SavePlotSource) {
      char name2_eff_pt_bar[100];
      if (doLayerEff) {
        sprintf (name2_eff_pt_bar, "_eff_pt_barrel_layer.C");
        c_eff_pt_barrel->SaveSource(plotDIR+plotPRE+outputname+name2_eff_pt_bar);
      }
      else {
        sprintf (name2_eff_pt_bar, "_eff_pt_barrel_module.C");
        c_eff_pt_barrel->SaveSource(plotDIR+plotPRE+outputname+name2_eff_pt_bar);
      }
    }


  // -----------------------------------------------------------------------------------------------------------
  // pT Endcap Summary
  // ---------------------------------------------------------------------------------------------------------

  TH2F *h_eff_pt_endcap    = new TH2F("eff_pt_endcap",";Particle p_{T} (GeV); Stub efficiency",300,0.,10.,100,0.,1.15);
  TH2F *h_D1_off   = new TH2F("D1_off","D1_off",300,0.,20.,100,0.,1.05);
  TH2F *h_D2_off   = new TH2F("D2_off","D2_off",300,0.,20.,100,0.,1.05);
  TH2F *h_D3_off   = new TH2F("D3_off","D3_off",300,0.,20.,100,0.,1.05);
  TH2F *h_D4_off   = new TH2F("D4_off","D4_off",300,0.,20.,100,0.,1.05);
  TH2F *h_D5_off   = new TH2F("D5_off","D5_off",300,0.,20.,100,0.,1.05);

  for (int i=0;i<100;++i)
  { 
    if (stub_pt[6][i]!=0)  h_D1_off->Fill(pt_val[i],stub_pt[6][i]);
    if (stub_pt[7][i]!=0)  h_D2_off->Fill(pt_val[i],stub_pt[7][i]);
    if (stub_pt[8][i]!=0)  h_D3_off->Fill(pt_val[i],stub_pt[8][i]);
    if (stub_pt[9][i]!=0)  h_D4_off->Fill(pt_val[i],stub_pt[9][i]);
    if (stub_pt[10][i]!=0) h_D5_off->Fill(pt_val[i],stub_pt[10][i]);
  }

  TCanvas *c_eff_pt_endcap = new TCanvas("c_eff_pt_endcap","pT Efficiency Endcap",0,0,640,640);
    c_eff_pt_endcap->Range(-1.2,-0.1,10.6,1.24);
    c_eff_pt_endcap->SetGridx();
    c_eff_pt_endcap->SetGridy();
    c_eff_pt_endcap->SetRightMargin(0.07);
    c_eff_pt_endcap->SetTopMargin(0.065);

    h_eff_pt_endcap->GetXaxis()->SetTitleSize(0.045);
    h_eff_pt_endcap->GetYaxis()->SetTitleSize(0.035);
    h_eff_pt_endcap->GetYaxis()->SetTitleOffset(0.88);
    h_eff_pt_endcap->GetXaxis()->SetTitleOffset(0.9);
    h_eff_pt_endcap->GetXaxis()->SetNdivisions(518);
    h_eff_pt_endcap->GetYaxis()->SetNdivisions(511);
    h_eff_pt_endcap->Draw();

    h_D1_off->SetMarkerStyle(24);
    h_D2_off->SetMarkerStyle(20);
    h_D3_off->SetMarkerStyle(27);
    h_D4_off->SetMarkerStyle(34);
    h_D5_off->SetMarkerStyle(28);
      
    h_D1_off->SetMarkerSize(1.4);
    h_D2_off->SetMarkerSize(1.4);
    h_D3_off->SetMarkerSize(1.4);
    h_D4_off->SetMarkerSize(1.4);
    h_D5_off->SetMarkerSize(1.4);

    h_D1_off->Draw("same");
    h_D2_off->Draw("same");
    h_D3_off->Draw("same");
    h_D4_off->Draw("same");
    h_D5_off->Draw("same");

    TLegend *leg_eff_pt_cap = new TLegend(0.55,0.2,0.9,0.45);
      leg_eff_pt_cap->SetTextSize(0.03);
      leg_eff_pt_cap->SetFillColor(0);
      leg_eff_pt_cap->AddEntry(h_D1_off,"TEDD double-disc 1","p");
      leg_eff_pt_cap->AddEntry(h_D2_off,"TEDD double-disc 2","p");
      leg_eff_pt_cap->AddEntry(h_D3_off,"TEDD double-disc 3","p");
      leg_eff_pt_cap->AddEntry(h_D4_off,"TEDD double-disc 4","p");
      leg_eff_pt_cap->AddEntry(h_D5_off,"TEDD double-disc 5","p");
      leg_eff_pt_cap->Draw();
         
    TLtitle03.DrawLatex(-0.02, 1.18, "CMS Phase-2 simulation");
    TL_PU.DrawLatex(5., 1.18, txt_pu);

    c_eff_pt_endcap->Modified();
    c_eff_pt_endcap->Update();
    c_eff_pt_endcap->Write();
      
    char name_eff_pt_cap[100];
    if (doLayerEff) {
      sprintf (name_eff_pt_cap, "_eff_pt_endcap_layer.png");
      c_eff_pt_endcap->SaveAs(plotDIR+plotPRE+outputname+name_eff_pt_cap);
    }
    else {
      sprintf (name_eff_pt_cap, "_eff_pt_endcap_module.png");
      c_eff_pt_endcap->SaveAs(plotDIR+plotPRE+outputname+name_eff_pt_cap);
    }

    if (SavePlotSource) {
      char name2_eff_pt_cap[100];
      if (doLayerEff) {
        sprintf (name2_eff_pt_cap, "_eff_pt_endcap_layer.C");
        c_eff_pt_endcap->SaveSource(plotDIR+plotPRE+outputname+name2_eff_pt_cap);
      }
      else {
        sprintf (name2_eff_pt_cap, "_eff_pt_endcap_module.C");
        c_eff_pt_endcap->SaveSource(plotDIR+plotPRE+outputname+name2_eff_pt_cap);
      }
    }



} // close PlotStubEffs




// ---------------------------------------------------------------------------------------------------------
// Other Functions
// ---------------------------------------------------------------------------------------------------------


  // ---------------------------------------------------------------------------------------------------------
  // Plot Style
  void SetPlotStyle() 
  {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasBorderSize(2);

    gStyle->SetFrameBorderMode(0);

    gStyle->SetLabelFont(42,"xyz");
    gStyle->SetTitleFont(42,"xyz");

    gStyle->SetLabelSize(0.035,"xyz");
    gStyle->SetTitleSize(0.035,"xyz");
  }
