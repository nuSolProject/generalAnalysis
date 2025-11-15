#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include <TStyle.h>

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>  // For exit()
#include <fstream>
#include <sstream>
#include <string>

struct histNFitData{

  //Histogram Options
  int nBins = 1000;
  double histMin = 0;
  double histMax = 350e3;
  std::string histName = "oldRP - 5.56+/-0.06% Resolution";
  bool setLogy = 1;

  //Fit Options
  bool attemptFit = 1;
  int nGaus = 1;
  double fitMin = 168e3;
  double fitMax = 191e3;

  //Exponential
  double expTerm1 = 0.;//5.
  double expTerm2 = -0.;//-10.

  //Gaussian 1
  double gaus1Amp = 1e3;
  double gaus1Mean = 165e3;
  double gaus1SD = 165e2;

  //Gaussian 2
  double gaus2Amp = 0;
  double gaus2Mean = 0;
  double gaus2SD = 0;

  //Gaussian 3
  double gaus3Amp = 0;
  double gaus3Mean = 0;
  double gaus3SD = 0;

  //Gaussian 4~
  double gaus4Amp = 0;
  double gaus4Mean = 0;
  double gaus4SD = 0;
 
};

histNFitData globalData;


void spectrumFitting(){


  std::string rootFilename = "data/oldRP/CH1137Cs.root";
  std::string outFilename = "137CsOld.png";
  
  TFile *file = TFile::Open(rootFilename.c_str());
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Could not open the ROOT file." << std::endl;
    return;
  }

  // Access the ntuple directly
  TTree* tree = (TTree*)file->Get("ntuple");
  if (!tree) {
    std::cerr << "Error: 'ntuple' not found in the ROOT file." << std::endl;
    file->Close();
    return;
  }

  // Branch variables
  Double_t peakVal;

  // Set the branch address
  tree->SetBranchAddress("fullInt1", &peakVal); //fullInt

  // Read entries from the tree and process them as needed
  Long64_t nEntries = tree->GetEntries();
  std::vector<double> peakValues;
  peakValues.reserve(nEntries);

  for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
    tree->GetEntry(iEntry);
    peakValues.push_back(peakVal);
  }

  std::string histogramTitle = globalData.histName
    + ";Waveform Integral (V*s);Entries";

  // Process peak values, create histogram, fit, etc.
  // This is an example of creating a simple histogram with these values
  TH1D *histogram =
    new TH1D("histogram",
	     histogramTitle.c_str(),
	     globalData.nBins, globalData.histMin, globalData.histMax);
    
  for (double value : peakValues) {
    histogram->Fill(value);
  }

  TCanvas *canvas = new TCanvas("canvas", "Histogram Canvas", 800, 600);
  histogram->Draw();

  //-----------------------------------------
  if(globalData.attemptFit){

    if(globalData.nGaus == 1){
      TF1 *fit = new TF1("fit", "gaus(0) + expo(3) + [5]",
			 globalData.fitMin, globalData.fitMax);

      fit-> SetParName(0, "gaus1Amp");
      fit-> SetParName(1, "gaus1Mean");
      fit-> SetParName(2, "gaus1SD");
      fit-> SetParName(3, "expLnScale");
      fit-> SetParName(4, "expConstant");
      fit-> SetParName(5, "constant");

      fit->SetParameters(globalData.gaus1Amp, globalData.gaus1Mean, globalData.gaus1SD,
			 globalData.expTerm1, globalData.expTerm2);
      histogram->Fit("fit", "", "", globalData.fitMin, globalData.fitMax);
    }
    else if(globalData.nGaus == 2){
      TF1 *fit = new TF1("fit", "gaus(0) + gaus(3) + expo(6) + [8]",
			 globalData.fitMin, globalData.fitMax);

      fit-> SetParName(0, "gaus1Amp");
      fit-> SetParName(1, "gaus1Mean");
      fit-> SetParName(2, "gaus1SD");
      fit-> SetParName(3, "gaus2Amp");
      fit-> SetParName(4, "gaus2Mean");
      fit-> SetParName(5, "gaus2SD");
      fit-> SetParName(6, "expLnScale");
      fit-> SetParName(7, "expConstant");
      fit-> SetParName(8, "constant");

      //fit-> SetParLimits(1, 12e-9, 12.5e-9);
      //fit-> SetParLimits(4, .09, .094);

      //fit-> SetParLimits(0, 2000, 2800);
      //fit-> SetParLimits(3, 50, 200);

      fit->SetParameters(globalData.gaus1Amp, globalData.gaus1Mean, globalData.gaus1SD,
			 globalData.gaus2Amp, globalData.gaus2Mean, globalData.gaus2SD,
			 globalData.expTerm1, globalData.expTerm2);
      histogram->Fit("fit", "", "", globalData.fitMin, globalData.fitMax);
    }
    else if(globalData.nGaus == 3){
      TF1 *fit = new TF1("fit", "gaus(0) + gaus(3) + gaus(6) + expo(9) + [11]",
			 globalData.fitMin, globalData.fitMax);

      fit-> SetParName(0, "gaus1Amp");
      fit-> SetParName(1, "gaus1Mean");
      fit-> SetParName(2, "gaus1SD");
      fit-> SetParName(3, "gaus2Amp");
      fit-> SetParName(4, "gaus2Mean");
      fit-> SetParName(5, "gaus2SD");
      fit-> SetParName(6, "gaus3Amp");
      fit-> SetParName(7, "gaus3Mean");
      fit-> SetParName(8, "gaus3SD");
      fit-> SetParName(9, "expLnScale");
      fit-> SetParName(10, "expConstant");
      fit-> SetParName(11, "constant");

      //amps
      //fit->SetParLimits(0, 1600, 1800);
      //fit->SetParLimits(3, 400, 800);
      //fit->SetParLimits(6, 200, 500);
      

      //fit->SetParLimits(4, .1025e-6, .105e-6);

      fit->SetParameters(globalData.gaus1Amp, globalData.gaus1Mean, globalData.gaus1SD,
			 globalData.gaus2Amp, globalData.gaus2Mean, globalData.gaus2SD,
			 globalData.gaus3Amp, globalData.gaus3Mean, globalData.gaus3SD,
			 globalData.expTerm1, globalData.expTerm2);
      histogram->Fit("fit", "", "", globalData.fitMin, globalData.fitMax);
    }
    else if(globalData.nGaus == 4){
      TF1 *fit =
	new TF1("fit", "gaus(0) + gaus(3) + gaus(6) + gaus(9) + expo(12) + [14]",
		globalData.fitMin, globalData.fitMax);

      fit-> SetParName(0, "gaus1Amp");
      fit-> SetParName(1, "gaus1Mean");
      fit-> SetParName(2, "gaus1SD");
      fit-> SetParName(3, "gaus2Amp");
      fit-> SetParName(4, "gaus2Mean");
      fit-> SetParName(5, "gaus2SD");
      fit-> SetParName(6, "gaus3Amp");
      fit-> SetParName(7, "gaus3Mean");
      fit-> SetParName(8, "gaus3SD");
      fit-> SetParName(9, "gaus4Amp");
      fit-> SetParName(10, "gaus4Mean");
      fit-> SetParName(11, "gaus4SD");
      fit-> SetParName(12, "expLnScale");
      fit-> SetParName(13, "expConstant");
      fit-> SetParName(14, "constant");

      //calculate abs from rel amps/means/sds
      globalData.gaus1Amp *= globalData.gaus3Amp;
      globalData.gaus1Mean *= globalData.gaus3Mean;
      globalData.gaus1SD *= globalData.gaus3SD;

      globalData.gaus2Amp *= globalData.gaus3Amp;
      globalData.gaus2Mean *= globalData.gaus3Mean;
      globalData.gaus2SD *= globalData.gaus3SD;

      globalData.gaus4Amp *= globalData.gaus3Amp;
      globalData.gaus4Mean *= globalData.gaus3Mean;
      globalData.gaus4SD *= globalData.gaus3SD;

      fit-> SetParLimits(0, globalData.gaus1Amp*.90, globalData.gaus1Amp*1.1);
      fit-> SetParLimits(1, globalData.gaus1Mean*.95, globalData.gaus1Mean*1.05);
      fit-> SetParLimits(2, globalData.gaus1SD*.95, globalData.gaus1SD*1.05);

      fit-> SetParLimits(3, globalData.gaus2Amp*.90, globalData.gaus2Amp*1.5);
      fit-> SetParLimits(4, globalData.gaus2Mean*.95, globalData.gaus2Mean*1.05);
      //fit-> SetParLimits(5, globalData.gaus2SD*.95, globalData.gaus2SD*1.05);

      fit-> SetParLimits(9, globalData.gaus4Amp*.90, globalData.gaus4Amp*1.1);
      fit-> SetParLimits(10, globalData.gaus4Mean*.95, globalData.gaus4Mean*1.05);
      //fit-> SetParLimits(11, globalData.gaus4SD*.95, globalData.gaus4SD*1.05);

      fit-> SetParLimits(6, 2000, 2500);

      //ROOT shits it's pants if you have 12 or more arguments for parameters. Why?
      //because fuck you, that's always why with ROOT
      double stupidFuckingFuckYouROOTArray[14] = {globalData.gaus1Amp, globalData.gaus1Mean,
	globalData.gaus1SD,globalData.gaus2Amp, globalData.gaus2Mean, globalData.gaus2SD,
	globalData.gaus3Amp, globalData.gaus3Mean, globalData.gaus3SD, globalData.gaus4Amp,
	globalData.gaus4Mean, globalData.gaus4SD, globalData.expTerm1, globalData.expTerm2};

      fit->SetParameters(stupidFuckingFuckYouROOTArray);
      histogram->Fit("fit", "", "", globalData.fitMin, globalData.fitMax);
    }
    else{
      std::cout << "1 <= nGaus <= 4. Wanna do more? Code it yourself bozo\n";
    }

    gStyle -> SetOptFit(1);
  }

  if(globalData.setLogy){
    canvas->SetLogy();
  }
  
  canvas->SaveAs(outFilename.c_str());

  // Cleanup
  delete canvas;
  delete histogram;
  file->Close();
}
