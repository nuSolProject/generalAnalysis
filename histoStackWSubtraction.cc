// compare.C
// Calibrated comparison of two runs with background subtraction.
// Run with:  root -l -q compare.C
//
// Adjust run times & energy range below.

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TError.h"
#include "TLine.h"
#include <iostream>
#include <algorithm>
#include <cmath>

gROOT->SetBatch(kTRUE);

void histoStackWSubtraction() {

  // ---------------- USER SETTINGS ----------------
  // Live times (seconds)
  double runTime1 = 359296.0; // background (AmBeNFI2)
  double runTime2 = 580960.0; // signal     (AmBeNFI3)

  // Energy histogram parameters (keV)
  int    nBins    = 400;
  double eMin_keV = 600.0;
  double eMax_keV = 1000.0;

  // Output scaling mode
  bool countsPerSecondPerKeV = false;
  // If true: scales additionally by (1 / binWidth) to give counts/s/keV.
  // If false: counts/s per bin (bin width shown implicitly).

  // Calibration coefficients (integral = a E^2 + b E + c0)
  const double a  = -7.049e-5;
  const double b  =  0.2975;
  const double c0 = -5.923;

  // File / tree / branch names
  const char* fileName1 = "AmBeNFI4.root"; //signal
  const char* fileName2 = "AmBeNFI2.root"; //bg
  const char* treeName  = "ntuple";
  const char* leafName  = "fullInt";

  // Subtraction order (signal - background):
  bool subtract_run1_from_run2 = true;
  // -----------------------------------------------

  if(runTime1 <= 0 || runTime2 <= 0) {
    std::cerr << "ERROR: Run times must be > 0.\n";
    return;
  }

  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kWarning;

  // Inverse calibration:
  // integral = aE^2 + bE + c0
  // Solve aE^2 + bE + (c0 - I) = 0 → E = (-b + sqrt(b^2 - 4 a (c0 - I))) / (2 a)
  // (We assume the physically meaningful root is with +sqrt.)
  TString energyExpr = Form("(-(%g) + sqrt((%g)*(%g) - 4*(%g)*((%g)-%s)))/(2*(%g))",
                            b, b, b, a, c0, leafName, a);

  // Open files
  TFile *f1 = TFile::Open(fileName1,"READ");
  TFile *f2 = TFile::Open(fileName2,"READ");
  if(!f1 || f1->IsZombie()){ std::cerr<<"Failed to open "<<fileName1<<"\n"; return; }
  if(!f2 || f2->IsZombie()){ std::cerr<<"Failed to open "<<fileName2<<"\n"; return; }

  TTree *t1 = dynamic_cast<TTree*>(f1->Get(treeName));
  TTree *t2 = dynamic_cast<TTree*>(f2->Get(treeName));
  if(!t1){ std::cerr<<"Missing tree '"<<treeName<<"' in "<<fileName1<<"\n"; return; }
  if(!t2){ std::cerr<<"Missing tree '"<<treeName<<"' in "<<fileName2<<"\n"; return; }

  // Clean any existing leftover histograms with these names
  if(gDirectory->FindObject("hE1")) gDirectory->Remove(gDirectory->FindObject("hE1"));
  if(gDirectory->FindObject("hE2")) gDirectory->Remove(gDirectory->FindObject("hE2"));

  // Create calibrated energy histograms
  TH1D *hE1 = new TH1D("hE1", ";Energy (keV);Counts / s / bin", nBins, eMin_keV, eMax_keV);
  TH1D *hE2 = new TH1D("hE2", ";Energy (keV);Counts / s / bin", nBins, eMin_keV, eMax_keV);

  // Draw command with calibration
  TString drawCmd1 = energyExpr + ">>hE1";
  TString drawCmd2 = energyExpr + ">>hE2";

  Long64_t nDraw1 = t1->Draw(drawCmd1, "", "goff");
  Long64_t nDraw2 = t2->Draw(drawCmd2, "", "goff");

  if(nDraw1 <= 0) std::cerr << "Warning: no entries filled for run 1.\n";
  if(nDraw2 <= 0) std::cerr << "Warning: no entries filled for run 2.\n";

  // Ensure sumw2
  if(!hE1->GetSumw2N()) hE1->Sumw2();
  if(!hE2->GetSumw2N()) hE2->Sumw2();

  // Normalize to counts per second per bin
  hE1->Scale(1.0 / runTime1);
  hE2->Scale(1.0 / runTime2);

  // Optional: convert to counts/s/keV
  if(countsPerSecondPerKeV) {
    double bw = hE1->GetXaxis()->GetBinWidth(1);
    hE1->Scale(1.0 / bw);
    hE2->Scale(1.0 / bw);
    hE1->GetYaxis()->SetTitle("Counts / s / keV");
    hE2->GetYaxis()->SetTitle("Counts / s / keV");
  }

  // Style
  hE1->SetLineColor(kBlue+1);
  hE1->SetLineWidth(2);

  hE2->SetLineColor(kRed+1);
  hE2->SetLineWidth(2);
  hE2->SetLineStyle(2);

  // Choose max
  double maxY = std::max(hE1->GetMaximum(), hE2->GetMaximum());
  hE1->SetMaximum(maxY * 1.4);
  hE1->SetMinimum(1e-5);

  // Difference histogram
  TH1D *hDiff = nullptr;
  hDiff = (TH1D*)hE1->Clone("hDiff");
  hDiff->SetTitle(";Energy (keV);(Source - Background)");
  hDiff->Add(hE2, -1.0);
    
  hDiff->SetLineColor(kBlack);
  hDiff->SetLineWidth(2);
  //hDiff->SetMarkerStyle(20);
  hDiff->SetMarkerColor(kBlack);

  // Determine symmetric Y range
  double maxAbs = 0.0;
  for(int b=1; b<=hDiff->GetNbinsX(); ++b)
    maxAbs = std::max(maxAbs, std::abs(hDiff->GetBinContent(b)));
  if(maxAbs <= 0) maxAbs = 1.0;
  hDiff->SetMaximum( 1.25 * maxAbs);
  hDiff->SetMinimum(-.25 * maxAbs);

  // Canvas & pads
  TCanvas *c = new TCanvas("c","Calibrated Spectrum Comparison",1000,800);
  c->Clear();
  TPad *pad1 = new TPad("pad1","pad1",0.0,0.28,1.0,1.0);
  TPad *pad2 = new TPad("pad2","pad2",0.0,0.0, 1.0,0.28);
  pad1->SetBottomMargin(0.02);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.30);
  pad1->Draw();
  pad2->Draw();

  // Top pad (log Y)
  pad1->cd();
  pad1->SetLogy();
  hE1->Draw("HIST");
  hE2->Draw("HIST SAME");

  TLegend *leg = new TLegend(0.63,0.68,0.88,0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hE1, Form("Source %s", fileName1), "l");
  leg->AddEntry(hE2, Form("Background %s", fileName2), "l");
  leg->Draw();

  double rawRate1 = (runTime1>0)? (double)nDraw1/runTime1 : 0;
  double rawRate2 = (runTime2>0)? (double)nDraw2/runTime2 : 0;

  TPaveText *pt = new TPaveText(0.15,0.68,0.40,0.88,"NDC");
  pt->SetFillColor(0);
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetTextFont(42);
  pt->SetTextSize(0.032);
  //pt->AddText(Form("Run1 time: %.1f s", runTime1));
  //pt->AddText(Form("Run2 time: %.1f s", runTime2));
  //pt->AddText(Form("Run1 entries: %lld (%.3g /s)", nDraw1, rawRate1));
  //pt->AddText(Form("Run2 entries: %lld (%.3g /s)", nDraw2, rawRate2));
  pt->Draw();

  // Bottom pad (difference)
  pad2->cd();
  hDiff->GetYaxis()->SetTitleOffset(0.55);
  hDiff->GetYaxis()->SetTitleSize(0.085);
  hDiff->GetYaxis()->SetLabelSize(0.075);
  hDiff->GetXaxis()->SetTitleSize(0.10);
  hDiff->GetXaxis()->SetLabelSize(0.085);
  hDiff->Draw("E1");

  TLine *zero = new TLine(eMin_keV,0,eMax_keV,0);
  zero->SetLineColor(kGray+2);
  zero->SetLineStyle(2);
  zero->Draw("SAME");

  c->SaveAs("comparison.png");

  // Save histograms
  TFile fout("comparison_subtracted.root","RECREATE");
  hE1->Write("hRun1_energy_norm");
  hE2->Write("hRun2_energy_norm");
  hDiff->Write(subtract_run1_from_run2 ? "hDiff_Run2_minus_Run1" : "hDiff_Run1_minus_Run2");
  fout.Close();

  f1->Close(); f2->Close();
  delete f1; delete f2;
}
