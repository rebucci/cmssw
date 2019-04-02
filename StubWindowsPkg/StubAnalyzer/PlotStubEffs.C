/*
    ROOT macro for stub efficiencies visualization

  Use:
      root[0]-> .L PlotStubEffs.C
      root[1]->PlotStubEffs("sourcefile","outputname",pu)
  where 
      - sourcefile is the name of the root output file produced by AM_ana. It does not include the ".root" suffix.
      - outputname is the name (no type suffix like ".root") of the root file this macro produces. it is also the prefix of all the plot names.
      - both input and output files have customizable preceding paths and prefixes. These options can be changed in the Preliminaries section of this macro.
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

void PlotStubEffs(TString sourcefile, TString outputname, int pu=0,  float ptmax=20.)
{

  // -----------------------------------------------------------------------------------------------------------
  // Preliminaries
  // -----------------------------------------------------------------------------------------------------------

  gROOT->SetBatch();
  gErrorIgnoreLevel = kError;
  SetPlotStyle();
  
  TString srcDIR  = "";
  TString subDIR  = "plots/efficiencies/";
  TString plotDIR = "plots/efficiencies/"; 
  TString plotPRE = "stubEff_";

  // Options
  bool doLayerEff     = true;   // module efficiency: 0, layer efficiency: 1 in seb's original
  bool PrintEachLayer = true;   // saves each layer plot
  bool doEtaPlots     = true;   // eta efficiencies
  bool SavePlotSource = false;  // save the macro.C files for the pT barrel and endcap summaries


  // Read in File
  TChain *newtree = new TChain("Efficiencies");
  newtree->Add(srcDIR+sourcefile+".root"); 

  // Output file
  TFile* fout = new TFile(subDIR+"plots_"+plotPRE+outputname+".root","recreate");


  // -----------------------------------------------------------------------------------------------------------
  // Definitions
  // -----------------------------------------------------------------------------------------------------------

	int binPtX  = 300; 	// 300
	int binPtY  = 100;  // 100
	int binEtaX = 75;  	// 75
	int binEtaY = 100;  // 100
	
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

  newtree->SetBranchAddress("pt_val", &pt_val);
  newtree->SetBranchAddress("digi_pt",&digi_pt);
  //newtree->SetBranchAddress("clus_pt",&clus_pt); 				// SV
  if (doLayerEff) {
    newtree->SetBranchAddress("clus_pt_lay",&clus_pt); 		// RB
    newtree->SetBranchAddress("stub_pt_lay",&stub_pt);
  }
  else {
    newtree->SetBranchAddress("clus_pt",&clus_pt); 				// RB
    newtree->SetBranchAddress("stub_pt",&stub_pt);
  }

  newtree->SetBranchAddress("eta_val", &eta_val);
  newtree->SetBranchAddress("digi_eta",&digi_eta);
  //newtree->SetBranchAddress("clus_eta",&clus_eta); 			// SV
  if (doLayerEff) {
    newtree->SetBranchAddress("clus_eta_lay",&clus_eta); 	// RB
    newtree->SetBranchAddress("stub_eta_lay",&stub_eta);
  }
  else {
    newtree->SetBranchAddress("clus_eta",&clus_eta); 			// RB
    newtree->SetBranchAddress("stub_eta",&stub_eta); 
  }

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
  // -----------------------------------------------------------------------------------------------------------

    // -----------------------------------------------------------------------------------------------------------
    // pT - Each Layer

      // Loop over Layers
      for (int layer=5;layer<25;++layer) {
        // Name Disk
        int disk =0;
        if (layer>=11  && layer<18) disk= int(layer-11)%8;
        if (layer>=18  && layer<25) disk=-int(layer-18)%8;

        // Histograms
        TH2F *h_eff_pt      = new TH2F("eff_pt",     ";Particle p_{T} (GeV); Efficiency",binPtX,0.,ptmax,binPtY,0.,1.15);
        TH2F *h_eff_pt_digi = new TH2F("eff_pt_digi",";Particle p_{T} (GeV); Efficiency",binPtX,0.,ptmax,binPtY,0.,1.15);
        TH2F *h_eff_pt_clus = new TH2F("eff_pt_clus",";Particle p_{T} (GeV); Efficiency",binPtX,0.,ptmax,binPtY,0.,1.15);
        TH2F *h_eff_pt_stub = new TH2F("eff_pt_stub",";Particle p_{T} (GeV); Efficiency",binPtX,0.,ptmax,binPtY,0.,1.15);

        // Fill Histograms
        for (int i=0;i<100;++i) { 
          h_eff_pt_digi->Fill(pt_val[i],digi_pt[layer][i]);
          h_eff_pt_clus->Fill(pt_val[i],clus_pt[layer][i]);
          h_eff_pt_stub->Fill(pt_val[i],stub_pt[layer][i]);
        }

        TCanvas *c_eff_pt = new TCanvas("c_eff_pt","pT Efficiency",0,0,640,640);
          c_eff_pt->Range(-1.2,-0.1,10.6,1.24);
          c_eff_pt->SetGridx();
          c_eff_pt->SetGridy();

          h_eff_pt->Draw();
          h_eff_pt->Write();
          
          h_eff_pt_digi->SetMarkerStyle(3);
          h_eff_pt_clus->SetMarkerStyle(4);
          h_eff_pt_stub->SetMarkerStyle(20);
          h_eff_pt_digi->Draw("same");
          h_eff_pt_digi->Write();
          h_eff_pt_clus->Draw("same");
          h_eff_pt_clus->Write();
          h_eff_pt_stub->Draw("same");
          h_eff_pt_stub->Write();

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
            leg_eff_pt->SetHeader(txt_eff_pt); 
            leg_eff_pt->AddEntry(h_eff_pt_digi,"Digis","p");
            leg_eff_pt->AddEntry(h_eff_pt_clus,"Clusters","p");
            leg_eff_pt->AddEntry(h_eff_pt_stub,"Stubs","p");
            leg_eff_pt->Draw();  

          TLtitle03.DrawLatex(0., 1.17, "CMS Phase-2 Simulation");
          TL_PU.DrawLatex(0.75*ptmax, 1.17, txt_pu);

          c_eff_pt->Modified();
          c_eff_pt->Update();
          c_eff_pt->Write();

          if (PrintEachLayer) {
            char name_eff_pt[100];
            if (doLayerEff) {
              if (disk==0) sprintf (name_eff_pt, "_eff_pt_layer_barrel_%d.png", layer);
              if (disk>0)  sprintf (name_eff_pt, "_eff_pt_layer_endcap_p%d.png", disk); 
              if (disk<0)  sprintf (name_eff_pt, "_eff_pt_layer_endcap_m%d.png", abs(disk)); 
            }
            else {
              if (disk==0) sprintf (name_eff_pt, "_eff_pt_module_barrel_%d.png", layer);
              if (disk>0)  sprintf (name_eff_pt, "_eff_pt_module_endcap_p%d.png", disk); 
              if (disk<0)  sprintf (name_eff_pt, "_eff_pt_module_endcap_m%d.png", abs(disk)); 
            }
            c_eff_pt->SaveAs(plotDIR+plotPRE+outputname+name_eff_pt);
          }

      }

    // -----------------------------------------------------------------------------------------------------------
    // pT Barrel Summary

      TH2F *h_eff_pt_barrel = new TH2F("eff_pt_barrel",";Particle p_{T} (GeV); Stub efficiency",binPtX,0.,ptmax,binPtY,0.,1.15);
      TH2F *h_L1_pt   = new TH2F("L1_pt","L1_pt",binPtX,0.,ptmax,binPtY,0.,1.15);
      TH2F *h_L2_pt   = new TH2F("L2_pt","L2_pt",binPtX,0.,ptmax,binPtY,0.,1.15);
      TH2F *h_L3_pt   = new TH2F("L3_pt","L3_pt",binPtX,0.,ptmax,binPtY,0.,1.15);
      TH2F *h_L4_pt   = new TH2F("L4_pt","L4_pt",binPtX,0.,ptmax,binPtY,0.,1.15);
      TH2F *h_L5_pt   = new TH2F("L5_pt","L5_pt",binPtX,0.,ptmax,binPtY,0.,1.15);
      TH2F *h_L6_pt   = new TH2F("L6_pt","L6_pt",binPtX,0.,ptmax,binPtY,0.,1.15);

      for (int i=0;i<100;++i)
      { 
        if (stub_pt[0][i]!=0) h_L1_pt->Fill(pt_val[i],stub_pt[0][i]);
        if (stub_pt[1][i]!=0) h_L2_pt->Fill(pt_val[i],stub_pt[1][i]);
        if (stub_pt[2][i]!=0) h_L3_pt->Fill(pt_val[i],stub_pt[2][i]);
        if (stub_pt[3][i]!=0) h_L4_pt->Fill(pt_val[i],stub_pt[3][i]);
        if (stub_pt[4][i]!=0) h_L5_pt->Fill(pt_val[i],stub_pt[4][i]);
        if (stub_pt[5][i]!=0) h_L6_pt->Fill(pt_val[i],stub_pt[5][i]);
      }

      TCanvas *c_eff_pt_barrel = new TCanvas("c_eff_pt_barrel","pT Efficiency Barrel",0,0,640,640);
        c_eff_pt_barrel->Range(-1.2,-0.1,10.6,1.24);
        c_eff_pt_barrel->SetGridx();
        c_eff_pt_barrel->SetGridy();

        h_eff_pt_barrel->Draw();
        
        h_L1_pt->SetMarkerStyle(24);
        h_L2_pt->SetMarkerStyle(20);
        h_L3_pt->SetMarkerStyle(27);
        h_L4_pt->SetMarkerStyle(34);
        h_L5_pt->SetMarkerStyle(28);
        h_L6_pt->SetMarkerStyle(25);
          
        h_L1_pt->SetMarkerSize(1.4);
        h_L2_pt->SetMarkerSize(1.4);
        h_L3_pt->SetMarkerSize(1.4);
        h_L4_pt->SetMarkerSize(1.4);
        h_L5_pt->SetMarkerSize(1.4);
        h_L6_pt->SetMarkerSize(1.4);

        h_L1_pt->Draw("same");
        h_L2_pt->Draw("same");
        h_L3_pt->Draw("same");
        h_L4_pt->Draw("same");
        h_L5_pt->Draw("same");
        h_L6_pt->Draw("same");

        TLegend *leg_eff_pt_bar = new TLegend(0.6,0.2,0.86,0.45);
          leg_eff_pt_bar->AddEntry(h_L1_pt,"TBPS layer 1","p");
          leg_eff_pt_bar->AddEntry(h_L2_pt,"TBPS layer 2","p");
          leg_eff_pt_bar->AddEntry(h_L3_pt,"TBPS layer 3","p");
          leg_eff_pt_bar->AddEntry(h_L4_pt,"TB2S layer 1","p");
          leg_eff_pt_bar->AddEntry(h_L5_pt,"TB2S layer 2","p");
          leg_eff_pt_bar->AddEntry(h_L6_pt,"TB2S layer 3","p");
          leg_eff_pt_bar->Draw();
           
        TLtitle03.DrawLatex(0., 1.17, "CMS Phase-2 simulation");
        TL_PU.DrawLatex(0.75*ptmax, 1.17, txt_pu);
        
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

      TH2F *h_eff_pt_endcap = new TH2F("eff_pt_endcap",";Particle p_{T} (GeV); Stub efficiency",binPtX,0.,ptmax,binPtY,0.,1.15);
      TH2F *h_D1_pt   = new TH2F("D1_pt","D1_pt",binPtX,0.,ptmax,binPtY,0.,1.15);
      TH2F *h_D2_pt   = new TH2F("D2_pt","D2_pt",binPtX,0.,ptmax,binPtY,0.,1.15);
      TH2F *h_D3_pt   = new TH2F("D3_pt","D3_pt",binPtX,0.,ptmax,binPtY,0.,1.15);
      TH2F *h_D4_pt   = new TH2F("D4_pt","D4_pt",binPtX,0.,ptmax,binPtY,0.,1.15);
      TH2F *h_D5_pt   = new TH2F("D5_pt","D5_pt",binPtX,0.,ptmax,binPtY,0.,1.15);

      for (int i=0;i<100;++i) { 
        if (stub_pt[6][i]!=0)  h_D1_pt->Fill(pt_val[i],stub_pt[6][i]);
        if (stub_pt[7][i]!=0)  h_D2_pt->Fill(pt_val[i],stub_pt[7][i]);
        if (stub_pt[8][i]!=0)  h_D3_pt->Fill(pt_val[i],stub_pt[8][i]);
        if (stub_pt[9][i]!=0)  h_D4_pt->Fill(pt_val[i],stub_pt[9][i]);
        if (stub_pt[10][i]!=0) h_D5_pt->Fill(pt_val[i],stub_pt[10][i]);
      }

      TCanvas *c_eff_pt_endcap = new TCanvas("c_eff_pt_endcap","pT Efficiency Endcap",0,0,640,640);
        c_eff_pt_endcap->Range(-1.2,-0.1,10.6,1.24);
        c_eff_pt_endcap->SetGridx();
        c_eff_pt_endcap->SetGridy();

        h_eff_pt_endcap->Draw();

        h_D1_pt->SetMarkerStyle(24);
        h_D2_pt->SetMarkerStyle(20);
        h_D3_pt->SetMarkerStyle(27);
        h_D4_pt->SetMarkerStyle(34);
        h_D5_pt->SetMarkerStyle(28);
          
        h_D1_pt->SetMarkerSize(1.4);
        h_D2_pt->SetMarkerSize(1.4);
        h_D3_pt->SetMarkerSize(1.4);
        h_D4_pt->SetMarkerSize(1.4);
        h_D5_pt->SetMarkerSize(1.4);

        h_D1_pt->Draw("same");
        h_D2_pt->Draw("same");
        h_D3_pt->Draw("same");
        h_D4_pt->Draw("same");
        h_D5_pt->Draw("same");

        TLegend *leg_eff_pt_cap = new TLegend(0.55,0.2,0.9,0.45);
          leg_eff_pt_cap->AddEntry(h_D1_pt,"TEDD double-disc 1","p");
          leg_eff_pt_cap->AddEntry(h_D2_pt,"TEDD double-disc 2","p");
          leg_eff_pt_cap->AddEntry(h_D3_pt,"TEDD double-disc 3","p");
          leg_eff_pt_cap->AddEntry(h_D4_pt,"TEDD double-disc 4","p");
          leg_eff_pt_cap->AddEntry(h_D5_pt,"TEDD double-disc 5","p");
          leg_eff_pt_cap->Draw();
             
        TLtitle03.DrawLatex(0., 1.17, "CMS Phase-2 simulation");
        TL_PU.DrawLatex(0.75*ptmax, 1.17, txt_pu);

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


  // -----------------------------------------------------------------------------------------------------------
  // Eta Plots
  // -----------------------------------------------------------------------------------------------------------
    if (doEtaPlots) {

    // -----------------------------------------------------------------------------------------------------------
    // Eta - Each Layer

      // Loop over Layers
      for (int layer=5;layer<25;++layer) {
        // Name Disk
        int disk =0;
        if (layer>=11  && layer<18) disk= int(layer-11)%8;
        if (layer>=18  && layer<25) disk=-int(layer-18)%8;

        // Histograms
        TH2F *h_eff_eta      = new TH2F("eff_eta",     ";Particle #eta; Efficiency",binEtaX,-2.5,2.5,binEtaY,0.,1.15);
        TH2F *h_eff_eta_digi = new TH2F("eff_eta_digi",";Particle #eta; Efficiency",binEtaX,-2.5,2.5,binEtaY,0.,1.15);
        TH2F *h_eff_eta_clus = new TH2F("eff_eta_clus",";Particle #eta; Efficiency",binEtaX,-2.5,2.5,binEtaY,0.,1.15);
        TH2F *h_eff_eta_stub = new TH2F("eff_eta_stub",";Particle #eta; Efficiency",binEtaX,-2.5,2.5,binEtaY,0.,1.15);
        
        // Fill Histograms
        for (int i=0;i<50;++i) { 
          h_eff_eta_digi->Fill(eta_val[i],digi_eta[layer][i]);
          h_eff_eta_clus->Fill(eta_val[i],clus_eta[layer][i]);
          h_eff_eta_stub->Fill(eta_val[i],stub_eta[layer][i]);
        }

        TCanvas *c_eff_eta = new TCanvas("c_eff_eta","eta Efficiency",0,0,640,640);
          c_eff_eta->Range(-1.2,-0.1,10.6,1.24);
          c_eff_eta->SetGridx();
          c_eff_eta->SetGridy();

          h_eff_eta->Draw();
          h_eff_eta->Write();
            
          h_eff_eta_digi->SetMarkerStyle(3);
          h_eff_eta_clus->SetMarkerStyle(4);
          h_eff_eta_stub->SetMarkerStyle(20);
          h_eff_eta_digi->Draw("same");
          h_eff_eta_digi->Write();
          h_eff_eta_clus->Draw("same");
          h_eff_eta_clus->Write();
          h_eff_eta_stub->Draw("same");        
          h_eff_eta_stub->Write();

          char txt_eff_eta[80];

          if (disk==0) {
            sprintf(txt_eff_eta, "Layer %d",layer);
          }
          else {
            (disk>0)
              ? sprintf(txt_eff_eta, "Disk %d",disk)
              : sprintf(txt_eff_eta, "Disk %d",-disk);
          }

          TLegend *leg_eff_eta = new TLegend(0.7,0.28,0.9,0.45);
            leg_eff_eta->SetHeader(txt_eff_eta); 
            leg_eff_eta->AddEntry(h_eff_eta_digi,"Digis","p");
            leg_eff_eta->AddEntry(h_eff_eta_clus,"Clusters","p");
            leg_eff_eta->AddEntry(h_eff_eta_stub,"Stubs","p");
            leg_eff_eta->Draw();

          TLtitle03.DrawLatex(-2.5, 1.17, "CMS Phase-2 Simulation");
          TL_PU.DrawLatex(1.25, 1.17, txt_pu);

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
              if (disk==0) sprintf (name_eff_eta, "_eff_eta_module_barrel_%d.png", layer);
              if (disk>0)  sprintf (name_eff_eta, "_eff_eta_module_endcap_p%d.png", disk); 
              if (disk<0)  sprintf (name_eff_eta, "_eff_eta_module_endcap_m%d.png", abs(disk)); 
            }
            c_eff_eta->SaveAs(plotDIR+plotPRE+outputname+name_eff_eta); 
          }

      }

    // -----------------------------------------------------------------------------------------------------------
    // Eta Barrel Summary

      TH2F *h_eff_eta_barrel = new TH2F("eff_eta_barrel",";Particle #eta; Stub efficiency",binEtaX,-2.5,2.5,binEtaY,0.,1.15);
      TH2F *h_L1_eta    = new TH2F("L1_eta","L1_eta",binEtaX,-2.5,2.5,binEtaY,0.,1.15);
      TH2F *h_L2_eta    = new TH2F("L2_eta","L2_eta",binEtaX,-2.5,2.5,binEtaY,0.,1.15);
      TH2F *h_L3_eta    = new TH2F("L3_eta","L3_eta",binEtaX,-2.5,2.5,binEtaY,0.,1.15);
      TH2F *h_L4_eta    = new TH2F("L4_eta","L4_eta",binEtaX,-2.5,2.5,binEtaY,0.,1.15);
      TH2F *h_L5_eta    = new TH2F("L5_eta","L5_eta",binEtaX,-2.5,2.5,binEtaY,0.,1.15);
      TH2F *h_L6_eta    = new TH2F("L6_eta","L6_eta",binEtaX,-2.5,2.5,binEtaY,0.,1.15);

      for (int i=0;i<50;++i) { 
        if (stub_eta[0][i]!=0) h_L1_eta ->Fill(eta_val[i],stub_eta[0][i]);
        if (stub_eta[1][i]!=0) h_L2_eta ->Fill(eta_val[i],stub_eta[1][i]);
        if (stub_eta[2][i]!=0) h_L3_eta ->Fill(eta_val[i],stub_eta[2][i]);
        if (stub_eta[3][i]!=0) h_L4_eta ->Fill(eta_val[i],stub_eta[3][i]);
        if (stub_eta[4][i]!=0) h_L5_eta ->Fill(eta_val[i],stub_eta[4][i]);
        if (stub_eta[5][i]!=0) h_L6_eta ->Fill(eta_val[i],stub_eta[5][i]);
      }

      TCanvas *c_eff_eta_barrel = new TCanvas("c_eff_eta_barrel","#eta Efficiency Barrel",0,0,640,640);
        c_eff_eta_barrel->Range(-1.2,-0.1,10.6,1.24);
        c_eff_eta_barrel->SetGridx();
        c_eff_eta_barrel->SetGridy();

        h_eff_eta_barrel->Draw();
        
        h_L1_eta->SetMarkerStyle(24);
        h_L2_eta->SetMarkerStyle(20);
        h_L3_eta->SetMarkerStyle(27);
        h_L4_eta->SetMarkerStyle(34);
        h_L5_eta->SetMarkerStyle(28);
        h_L6_eta->SetMarkerStyle(25);
          
        h_L1_eta->SetMarkerSize(1.4);
        h_L2_eta->SetMarkerSize(1.4);
        h_L3_eta->SetMarkerSize(1.4);
        h_L4_eta->SetMarkerSize(1.4);
        h_L5_eta->SetMarkerSize(1.4);
        h_L6_eta->SetMarkerSize(1.4);

        h_L1_eta->Draw("same");
        h_L2_eta->Draw("same");
        h_L3_eta->Draw("same");
        h_L4_eta->Draw("same");
        h_L5_eta->Draw("same");
        h_L6_eta->Draw("same");

        TLegend *leg_eff_eta_bar = new TLegend(0.6,0.2,0.86,0.45);
          leg_eff_eta_bar->AddEntry(h_L1_eta,"TBPS layer 1","p");
          leg_eff_eta_bar->AddEntry(h_L2_eta,"TBPS layer 2","p");
          leg_eff_eta_bar->AddEntry(h_L3_eta,"TBPS layer 3","p");
          leg_eff_eta_bar->AddEntry(h_L4_eta,"TB2S layer 1","p");
          leg_eff_eta_bar->AddEntry(h_L5_eta,"TB2S layer 2","p");
          leg_eff_eta_bar->AddEntry(h_L6_eta,"TB2S layer 3","p");
          leg_eff_eta_bar->Draw();
           
        TLtitle03.DrawLatex(-2.5, 1.17, "CMS Phase-2 simulation");
        TL_PU.DrawLatex(1.25, 1.17, txt_pu);
        
        c_eff_eta_barrel->Modified();
        c_eff_eta_barrel->Update();
        c_eff_eta_barrel->Write();

        char name_eff_eta_bar[100];
        if (doLayerEff) {
          sprintf (name_eff_eta_bar, "_eff_eta_barrel_layer.png");
          c_eff_eta_barrel->SaveAs(plotDIR+plotPRE+outputname+name_eff_eta_bar);
        }
        else {
          sprintf (name_eff_eta_bar, "_eff_eta_barrel_module.png");
          c_eff_eta_barrel->SaveAs(plotDIR+plotPRE+outputname+name_eff_eta_bar);
        }

        if (SavePlotSource) {
          char name2_eff_eta_bar[100];
          if (doLayerEff) {
            sprintf (name2_eff_eta_bar, "_eff_eta_barrel_layer.C");
            c_eff_eta_barrel->SaveSource(plotDIR+plotPRE+outputname+name2_eff_eta_bar);
          }
          else {
            sprintf (name2_eff_eta_bar, "_eff_eta_barrel_module.C");
            c_eff_eta_barrel->SaveSource(plotDIR+plotPRE+outputname+name2_eff_eta_bar);
          }
        }

    // -----------------------------------------------------------------------------------------------------------
    // Eta Endcap Summary

      TH2F *h_eff_eta_endcap    = new TH2F("eff_eta_endcap",";Particle #eta; Stub efficiency",binEtaX,-2.5,2.5,binEtaY,0.,1.15);
      TH2F *h_D1_eta   = new TH2F("D1_eta","D1_eta",binEtaX,-2.5,2.5,binEtaY,0.,1.15);
      TH2F *h_D2_eta   = new TH2F("D2_eta","D2_eta",binEtaX,-2.5,2.5,binEtaY,0.,1.15);
      TH2F *h_D3_eta   = new TH2F("D3_eta","D3_eta",binEtaX,-2.5,2.5,binEtaY,0.,1.15);
      TH2F *h_D4_eta   = new TH2F("D4_eta","D4_eta",binEtaX,-2.5,2.5,binEtaY,0.,1.15);
      TH2F *h_D5_eta   = new TH2F("D5_eta","D5_eta",binEtaX,-2.5,2.5,binEtaY,0.,1.15);

      for (int i=0;i<50;++i) { 
        if (stub_eta[6][i]!=0)  h_D1_eta->Fill(eta_val[i],stub_eta[6][i]);
        if (stub_eta[7][i]!=0)  h_D2_eta->Fill(eta_val[i],stub_eta[7][i]);
        if (stub_eta[8][i]!=0)  h_D3_eta->Fill(eta_val[i],stub_eta[8][i]);
        if (stub_eta[9][i]!=0)  h_D4_eta->Fill(eta_val[i],stub_eta[9][i]);
        if (stub_eta[10][i]!=0) h_D5_eta->Fill(eta_val[i],stub_eta[10][i]);
      }

      TCanvas *c_eff_eta_endcap = new TCanvas("c_eff_eta_endcap","#eta Efficiency Endcap",0,0,640,640);
        c_eff_eta_endcap->Range(-1.2,-0.1,10.6,1.24);
        c_eff_eta_endcap->SetGridx();
        c_eff_eta_endcap->SetGridy();

        h_eff_eta_endcap->Draw();

        h_D1_eta->SetMarkerStyle(24);
        h_D2_eta->SetMarkerStyle(20);
        h_D3_eta->SetMarkerStyle(27);
        h_D4_eta->SetMarkerStyle(34);
        h_D5_eta->SetMarkerStyle(28);
          
        h_D1_eta->SetMarkerSize(1.4);
        h_D2_eta->SetMarkerSize(1.4);
        h_D3_eta->SetMarkerSize(1.4);
        h_D4_eta->SetMarkerSize(1.4);
        h_D5_eta->SetMarkerSize(1.4);

        h_D1_eta->Draw("same");
        h_D2_eta->Draw("same");
        h_D3_eta->Draw("same");
        h_D4_eta->Draw("same");
        h_D5_eta->Draw("same");

        TLegend *leg_eff_eta_cap = new TLegend(0.55,0.2,0.9,0.45);
          leg_eff_eta_cap->AddEntry(h_D1_eta,"TEDD double-disc 1","p");
          leg_eff_eta_cap->AddEntry(h_D2_eta,"TEDD double-disc 2","p");
          leg_eff_eta_cap->AddEntry(h_D3_eta,"TEDD double-disc 3","p");
          leg_eff_eta_cap->AddEntry(h_D4_eta,"TEDD double-disc 4","p");
          leg_eff_eta_cap->AddEntry(h_D5_eta,"TEDD double-disc 5","p");
          leg_eff_eta_cap->Draw();
             
        TLtitle03.DrawLatex(-2.5, 1.17, "CMS Phase-2 simulation");
        TL_PU.DrawLatex(1.25, 1.17, txt_pu);

        c_eff_eta_endcap->Modified();
        c_eff_eta_endcap->Update();
        c_eff_eta_endcap->Write();
          
        char name_eff_eta_cap[100];
        if (doLayerEff) {
          sprintf (name_eff_eta_cap, "_eff_eta_endcap_layer.png");
          c_eff_eta_endcap->SaveAs(plotDIR+plotPRE+outputname+name_eff_eta_cap);
        }
        else {
          sprintf (name_eff_eta_cap, "_eff_eta_endcap_module.png");
          c_eff_eta_endcap->SaveAs(plotDIR+plotPRE+outputname+name_eff_eta_cap);
        }

        if (SavePlotSource) {
          char name2_eff_eta_cap[100];
          if (doLayerEff) {
            sprintf (name2_eff_eta_cap, "_eff_eta_endcap_layer.C");
            c_eff_eta_endcap->SaveSource(plotDIR+plotPRE+outputname+name2_eff_eta_cap);
          }
          else {
            sprintf (name2_eff_eta_cap, "_eff_eta_endcap_module.C");
            c_eff_eta_endcap->SaveSource(plotDIR+plotPRE+outputname+name2_eff_eta_cap);
          }
        }


    } // close doEtaPlots



} // close PlotStubEffs




// ---------------------------------------------------------------------------------------------------------
// Other Functions
// ---------------------------------------------------------------------------------------------------------


  // ---------------------------------------------------------------------------------------------------------
  // Plot Style
  void SetPlotStyle() 
  {
    // Canvas and Pad Settings
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasBorderSize(2);
    gStyle->SetFrameBorderMode(0);
    
    gStyle->SetPadTopMargin(0.065);
    gStyle->SetPadRightMargin(0.07);

    // Axis Label Settings
    gStyle->SetLabelFont(42,"xyz");
    gStyle->SetLabelSize(0.03,"xyz");
    gStyle->SetNdivisions(518,"x");
    gStyle->SetNdivisions(511,"y");

    // Axis Title Settings
    gStyle->SetTitleFont(42,"xyz");
    gStyle->SetTitleSize(0.035,"xyz");
    gStyle->SetTitleXOffset(0.9);
    gStyle->SetTitleYOffset(0.8);

    // Legend Options
    gStyle->SetLegendTextSize(0.03);
    gStyle->SetFillColor(0);
  }

