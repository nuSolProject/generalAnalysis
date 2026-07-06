#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <ctime>

//root includes
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TAxis.h"
#include "Math/MinimizerOptions.h"
#include "TNtupleD.h"
#include "TFile.h"
#include "TError.h"

const int eventsPerFile = 1e5;//1e5
const int samplesPerEvent = 1024;
double RPCHData[eventsPerFile][2][samplesPerEvent];
bool savePngs = 0;

struct outData{
  std::string timeStamp;
  double fullInt1;
  double trigTime1;
  double smartInt1;
  
  double fullInt2;
  double trigTime2;
  double smartInt2;
  
  outData(std::string tS, double fI1, double tT1, double sI1,
	  double fI2, double tT2, double sI2):
    timeStamp(tS), fullInt1(fI1), trigTime1(tT1), smartInt1(sI1),
    fullInt2(fI2), trigTime2(tT2), smartInt2(sI2) {}
};

/*
  f1.SetParameter(0, p0Guess1);
  f1.SetParameter(1, tauRiseGuess);
  f1.SetParameter(2, tauDecay1Guess);
  f1.SetParameter(3, p0Guess1/6.);
  f1.SetParameter(4, tauDecay2Guess);
  f1.SetParameter(5, triggerCH1Guess);
  f1.SetParameter(6, 0.0);*/
struct fitData{
  double ampFastCH1;
  double tauRiseCH1;
  double tauDecay1CH1;
  double ampSlowCH1;
  double tauDecay2CH1;
  double trigCH1;
  double bslnCH1;

  double ampFastCH2;
  double tauRiseCH2;
  double tauDecay1CH2;
  double ampSlowCH2;
  double tauDecay2CH2;
  double trigCH2;
  double bslnCH2;

  fitData(double aF1, double tR1, double tD11, double aS1, double tD21, double t1, double b1,
	  double aF2, double tR2, double tD12, double aS2, double tD22, double t2, double b2):
    ampFastCH1(aF1), tauRiseCH1(tR1), tauDecay1CH1(tD11), ampSlowCH1(aS1), tauDecay2CH1(tD21), trigCH1(t1), bslnCH1(b1),
    ampFastCH2(aF2), tauRiseCH2(tR2), tauDecay1CH2(tD12), ampSlowCH2(aS2), tauDecay2CH2(tD22), trigCH2(t2), bslnCH2(b2)
  {}
};



std::vector<outData> outDataList;
std::vector<fitData> fitDataList;

static std::string nsToIso8601(uint64_t ns_since_epoch) {
  time_t sec = static_cast<time_t>(ns_since_epoch / 1000000000ull);
  uint32_t nsec = static_cast<uint32_t>(ns_since_epoch % 1000000000ull);
  std::tm tm_utc{};
#if defined(_WIN32)
  gmtime_s(&tm_utc, &sec);
#else
  gmtime_r(&sec, &tm_utc);
#endif
  std::ostringstream oss;
  oss << std::put_time(&tm_utc, "%Y-%m-%dT%H:%M:%S") << '.'
      << std::setw(9) << std::setfill('0') << nsec << "Z";
  return oss.str();
}

bool readStoreData(std::string curFilename){
#pragma pack(push,1)
  struct FileHeader {
    char     magic[8];          // "RP2CHV10"
    uint16_t version;
    uint16_t header_size;       // 256
    uint8_t  endianness;        // 1 = little
    uint8_t  channels;          // 2
    uint8_t  sample_format;     // 2 = float32 interleaved
    uint8_t  bits_per_sample;   // 32
    uint32_t nsamples;          // 1024
    uint32_t presamples;        // 64
    uint32_t decimation;        // 1
    uint32_t sample_rate_hz;
    uint32_t sample_period_ps;
    uint32_t trigger_src;       // RP_TRIG_SRC_CHA/CHB_PE
    float    trigger_level_v;
    uint64_t file_start_time_ns;
    char     run_note[128];
    uint8_t  reserved[76];
  };
  struct EventHeader {
    uint64_t timestamp_ns;
    uint64_t seq_no;
    uint32_t tpos;
    uint16_t flags;
    uint32_t payload_bytes;     // 8192
    uint8_t  reserved[38];
  };
#pragma pack(pop)
  static_assert(sizeof(FileHeader)==256, "fh size");
  static_assert(sizeof(EventHeader)==64, "eh size");

  std::ifstream in(curFilename, std::ios::binary);
  if (!in) return false;

  FileHeader fh{};
  if (!in.read(reinterpret_cast<char*>(&fh), sizeof(fh))) return false;
  if (std::memcmp(fh.magic, "RP2CHV10", 8) != 0) return false;
  if (fh.channels != 2 || fh.sample_format != 2 || fh.bits_per_sample != 32) return false;
  if (fh.nsamples != samplesPerEvent) return false;

  const size_t floatsPerEvent = static_cast<size_t>(fh.nsamples) * 2u;
  std::vector<float> payload(floatsPerEvent);

  int e = 0;
  while (e < eventsPerFile) {
    if (in.peek() == std::char_traits<char>::eof()) break;

    EventHeader eh{};
    if (!in.read(reinterpret_cast<char*>(&eh), sizeof(eh))) break;
    if (eh.payload_bytes != fh.nsamples * 2u * sizeof(float)) break;
    if (!in.read(reinterpret_cast<char*>(payload.data()), eh.payload_bytes)) break;

    // Demux: CH1 -> RPCHData[e][0][i], CH2 -> RPCHData[e][1][i]
    for (int i = 0; i < samplesPerEvent; ++i) {
      RPCHData[e][0][i] = static_cast<double>(payload[2*i + 0]); // CH1
      RPCHData[e][1][i] = static_cast<double>(payload[2*i + 1]); // CH2
    }

    outDataList.emplace_back(nsToIso8601(eh.timestamp_ns),
			     /*fullInt=*/0.0, /*trigTime=*/0.0, /*smartInt=*/0.0, 0.0, 0.0, 0.0);
    ++e;
  }
  return (e > 0);
}


void coutRPCHData(){

  std::cout << outDataList[0].timeStamp << ":\n";
  for(int i = 0; i < eventsPerFile; ++i){
    std::cout << RPCHData[0][0][i] << "\n";
  }

  
}


void computeFullInts(){


  for(int eventNo = 0; eventNo < eventsPerFile; ++eventNo){

    //CH1
    double fullInt = 0.0;
    for(int sampleNo = 0; sampleNo < samplesPerEvent-1; ++sampleNo){
      fullInt += RPCHData[eventNo][0][sampleNo];
      fullInt += RPCHData[eventNo][0][sampleNo+1];
    }
    fullInt *= 0.5 * 8.0; //1/2 * 8ns * total

    outDataList[eventNo].fullInt1 = fullInt;

    //CH2
    fullInt = 0.0;
    for(int sampleNo = 0; sampleNo < samplesPerEvent-1; ++sampleNo){
      fullInt += RPCHData[eventNo][1][sampleNo];
      fullInt += RPCHData[eventNo][1][sampleNo+1];
    }
    fullInt *= 0.5 * 8.0; //1/2 * 8ns * total

    outDataList[eventNo].fullInt2 = fullInt;
    
  }

}


void coutFullInts(){


  for(int i = 0; i < eventsPerFile; ++i){
    std::cout << outDataList[i].fullInt1 << ","
	      << outDataList[i].fullInt2 << "\n";
  }

  
}


int findPeakIndex(int curEvent, int channel){


  int peakIndex = 0;
  
  double maxSample = RPCHData[curEvent][channel][0];
  for(int sampleNo = 1; sampleNo < samplesPerEvent; ++sampleNo){
    if(RPCHData[curEvent][channel][sampleNo] > maxSample){
      maxSample = RPCHData[curEvent][channel][sampleNo];
      peakIndex = sampleNo;
    }
  }

  //std::cout << peakIndex << ",";
  

  return peakIndex;
}


int findTrigIndex(int curEvent, int channel){


  int trigIndex = findPeakIndex(curEvent, channel);
  double trigger = 0.005; //5 mv
  
  while(RPCHData[curEvent][channel][trigIndex] > trigger){
    --trigIndex;
  }

  //std::cout << trigIndex << "\n";
  
  if(channel == 0){
    outDataList[curEvent].trigTime1 = trigIndex;
  }
  else{
    outDataList[curEvent].trigTime2 = trigIndex;
  }
  

  return trigIndex;
}


void computeSmartInts(){


  int preSamples = 8;
  int totalSamples = 512;
  
  for(int eventNo = 0; eventNo < eventsPerFile; ++eventNo){

    //CH1
    double fullInt = 0.0;
    int trigIndex = findTrigIndex(eventNo, 0);

    if(trigIndex >= preSamples && (trigIndex + totalSamples) < samplesPerEvent){
    
      for(int sampleNo = trigIndex-preSamples; sampleNo < trigIndex + totalSamples - 1; ++sampleNo){
	fullInt += RPCHData[eventNo][0][sampleNo];
	fullInt += RPCHData[eventNo][0][sampleNo+1];
      }
      fullInt *= 0.5 * 8.0; //1/2 * 8ns * total

    }

    outDataList[eventNo].smartInt1 = fullInt;

    //CH2
    fullInt = 0.0;
    trigIndex = findTrigIndex(eventNo, 1);

    if(trigIndex >= preSamples && (trigIndex + totalSamples) < samplesPerEvent){
    
      for(int sampleNo = trigIndex-preSamples; sampleNo < trigIndex + totalSamples - 1; ++sampleNo){
	fullInt += RPCHData[eventNo][1][sampleNo];
	fullInt += RPCHData[eventNo][1][sampleNo+1];
      }
      fullInt *= 0.5 * 8.0; //1/2 * 8ns * total

    }

    outDataList[eventNo].smartInt2 = fullInt;
    
  }
  

}


void coutSmartInts(){


  for(int i = 0; i < eventsPerFile; ++i){
    std::cout << outDataList[i].smartInt1 << ","
	      << outDataList[i].smartInt2 << "\n";
  }

  
}





void storeDataInNtuple(std::string outFilename, const std::vector<outData>& dataList, const std::vector<fitData>& fitList){

  
  // Create a TNtupleD
  TNtupleD* nt1 = new TNtupleD("nt1", "Data Ntuple", "timeStamp:fullInt1:trigTime1:smartInt1:fullInt2:trigTime2:smartInt2");
  // Fill the ntuple
  for(const auto& data : dataList){
    Double_t values[] = {std::stod(data.timeStamp), data.fullInt1, data.trigTime1,
      data.smartInt1, data.fullInt2, data.trigTime2, data.smartInt2};
    nt1->Fill(values);
  }

  // Create a TNtupleD
  TNtupleD* nt2 = new TNtupleD("nt2", "Data Ntuple", "ampFastCH1:tauRiseCH1:tauDecay1CH1:ampSlowCH1:tauDecay2CH1:trigCH1:bslnCH1:ampFastCH2:tauRiseCH2:tauDecay1CH2:ampSlowCH2:tauDecay2CH2:trigCH2:bslnCH2");
  // Fill the ntuple
  for(const auto& data : fitList){
    Double_t values[] = {data.ampFastCH1, data.tauRiseCH1, data.tauDecay1CH1, data.ampSlowCH1, data.tauDecay2CH1, data.trigCH1, data.bslnCH1, data.ampFastCH2, data.tauRiseCH2, data.tauDecay1CH2, data.ampSlowCH2, data.tauDecay2CH2, data.trigCH2, data.bslnCH2};
    nt2->Fill(values);
  }
  
  // Save the ntuple to a ROOT file using the global outputFilename
  TFile* outFile = new TFile(outFilename.c_str(), "recreate"); //recreate
  nt1->Write();
  nt2->Write();
  outFile->Close(); // Ensure we properly close and release the file
  // Clean up
  delete outFile; // Make sure to delete the TFile object to free memory
  delete nt1; // Also delete the TNtupleD to prevent memory leaks
  delete nt2;


}


void computeFitEvent(int eventNo){

  std::vector<double> RPTimeVec(samplesPerEvent);
  std::vector<double> RPCH1DataVec(samplesPerEvent);
  std::vector<double> RPCH2DataVec(samplesPerEvent);

  for(int i = 0; i < samplesPerEvent; ++i){
    RPTimeVec[i] = i*8e-9;
    RPCH1DataVec[i] = RPCHData[eventNo][0][i];
    RPCH2DataVec[i] = RPCHData[eventNo][1][i];
  }

  // parameter guesses
  int CH1PeakIndex = findPeakIndex(eventNo, 0);
  double p0Guess1 = std::abs(RPCHData[eventNo][0][CH1PeakIndex])*2.0;

  int CH2PeakIndex = findPeakIndex(eventNo, 1);
  double p0Guess2 = std::abs(RPCHData[eventNo][1][CH2PeakIndex])*2.0;

  double p0Min = 0.0;
  double p0Max = 20.0;

  double tauRiseGuess = 45e-9;
  double tauRiseMin   = 10e-9;
  double tauRiseMax   = 80e-9;

  double tauDecay1Guess = 92e-9;
  double tauDecay1Min   = 50e-9;
  double tauDecay1Max   = 200e-9;

  double tauDecay2Guess = 750e-9;
  double tauDecay2Min   = 500e-9;
  double tauDecay2Max   = 1000e-9;

  double triggerCH1Guess = (outDataList[eventNo].trigTime1 * 8e-9) - (8e-9);
  double triggerCH2Guess = (outDataList[eventNo].trigTime2 * 8e-9) - (8e-9);

  double triggerMin = -8e-9;
  double triggerMax = 8e-9;

  const double xMin = RPTimeVec.front();
  const double xMax = RPTimeVec.back();

  const char *fitExpr =
    "(x>[5]) ? [6] + [0]*(exp(-(x-[5])/[2]) - exp(-(x-[5])/[1]))"
    " + [3]*(exp(-(x-[5])/[4]) - exp(-(x-[5])/[1])) : [6]";

  // -----------------------
  // CH1
  // -----------------------
  TGraph g1(samplesPerEvent, RPTimeVec.data(), RPCH1DataVec.data());

  TF1 f1(Form("f1_event_%d", eventNo), fitExpr, xMin, xMax);

  f1.SetParameter(0, p0Guess1);
  f1.SetParameter(1, tauRiseGuess);
  f1.SetParameter(2, tauDecay1Guess);
  f1.SetParameter(3, p0Guess1/6.);
  f1.SetParameter(4, tauDecay2Guess);
  f1.SetParameter(5, triggerCH1Guess);
  f1.SetParameter(6, 0.0);

  f1.SetParLimits(0, p0Min, p0Max);
  f1.SetParLimits(1, tauRiseMin, tauRiseMax);
  f1.SetParLimits(2, tauDecay1Min, tauDecay1Max);
  f1.SetParLimits(3, p0Min, p0Max/6.);
  f1.SetParLimits(4, tauDecay2Min, tauDecay2Max);
  f1.SetParLimits(5, triggerCH1Guess - triggerMin, triggerCH1Guess + triggerMax);

  g1.Fit(&f1, "QN0", "", xMin, xMax);

  // -----------------------
  // CH2
  // -----------------------
  TGraph g2(samplesPerEvent, RPTimeVec.data(), RPCH2DataVec.data());

  TF1 f2(Form("f2_event_%d", eventNo), fitExpr, xMin, xMax);

  f2.SetParameter(0, p0Guess2);
  f2.SetParameter(1, tauRiseGuess);
  f2.SetParameter(2, tauDecay1Guess);
  f2.SetParameter(3, p0Guess2/6.);
  f2.SetParameter(4, tauDecay2Guess);
  f2.SetParameter(5, triggerCH2Guess);
  f2.SetParameter(6, 0.0);

  f2.SetParLimits(0, p0Min, p0Max);
  f2.SetParLimits(1, tauRiseMin, tauRiseMax);
  f2.SetParLimits(2, tauDecay1Min, tauDecay1Max);
  f2.SetParLimits(3, p0Min, p0Max/6.);
  f2.SetParLimits(4, tauDecay2Min, tauDecay2Max);
  f2.SetParLimits(5, triggerCH2Guess - triggerMin, triggerCH2Guess + triggerMax);

  g2.Fit(&f2, "QN0", "", xMin, xMax);

  //std::cout << "Storing parameters into fitDataList\n";

  //output parameters to fitDataList
  fitDataList.emplace_back(f1.GetParameter(0), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  fitDataList[eventNo].tauRiseCH1 = f1.GetParameter(1);
  fitDataList[eventNo].tauDecay1CH1 = f1.GetParameter(2);
  fitDataList[eventNo].ampSlowCH1 = f1.GetParameter(3);
  fitDataList[eventNo].tauDecay2CH1 = f1.GetParameter(4);
  fitDataList[eventNo].trigCH1 = f1.GetParameter(5);
  fitDataList[eventNo].bslnCH1 = f1.GetParameter(6);

  fitDataList[eventNo].ampFastCH2 = f2.GetParameter(0);
  fitDataList[eventNo].tauRiseCH2 = f2.GetParameter(1);
  fitDataList[eventNo].tauDecay1CH2 = f2.GetParameter(2);
  fitDataList[eventNo].ampSlowCH2 = f2.GetParameter(3);
  fitDataList[eventNo].tauDecay2CH2 = f2.GetParameter(4);
  fitDataList[eventNo].trigCH2 = f2.GetParameter(5);
  fitDataList[eventNo].bslnCH2 = f2.GetParameter(6);

  //std::cout << "Finished storing parameters into fitDataList\n";

  if(savePngs){
    
    TCanvas *c = new TCanvas(Form("c_event_%d", eventNo), "", 1000, 700);

    g1.SetTitle(Form("Event %d;Time (s);Voltage (V)", eventNo));

    g1.SetMarkerStyle(20);
    g1.SetMarkerSize(0.5);
    g1.SetLineColor(kBlue);
    g1.SetMarkerColor(kBlue);

    g2.SetMarkerStyle(20);
    g2.SetMarkerSize(0.5);
    g2.SetLineColor(kRed);
    g2.SetMarkerColor(kRed);

    f1.SetLineColor(kBlue + 2);
    f1.SetLineWidth(2);

    f2.SetLineColor(kRed + 2);
    f2.SetLineWidth(2);

    g1.Draw("AP");
    g2.Draw("P SAME");

    f1.Draw("SAME");
    f2.Draw("SAME");

    TLegend *leg = new TLegend(0.65, 0.72, 0.88, 0.88);
    leg->AddEntry(&g1, "CH1 data", "p");
    leg->AddEntry(&f1, "CH1 fit", "l");
    leg->AddEntry(&g2, "CH2 data", "p");
    leg->AddEntry(&f2, "CH2 fit", "l");
    leg->Draw();

    c->SaveAs(Form("%d.png", eventNo));

    delete leg;
    delete c;
  }
  
}


int RPBinDecoder(){


  if((eventsPerFile > 1e3) && savePngs){
    std::cout << "Why are you saving " << eventsPerFile << " pngs?\n";
    return 1;
  }

  for(int runNo = 1; runNo < 2; ++runNo){// runNo > 851

    outDataList.clear();
    fitDataList.clear();

    std::string curFilename = "/home/bhartsock/RP/RPBinDecoder/data/22Na/";
    curFilename += "run";
    curFilename += std::to_string(runNo);
    curFilename += ".bin";
  
    if(!readStoreData(curFilename)){
      std::cout << "Can't open " << curFilename << "\n";
      return 1;
    }
    //coutRPCHData();

    computeFullInts();
    //coutFullInts();

    /*
      for(int i = 0; i < eventsPerFile; ++i){
      findTrigIndex(i, 0);
      }*/
    computeSmartInts();
    //coutSmartInts();

    for(int i = 0; i < eventsPerFile; ++i){
      computeFitEvent(i);
    }

    //std::cout << "Finished fits, outputting data\n";

    std::string outFilename = "22Na";
    outFilename += std::to_string(runNo) + ".root";

    storeDataInNtuple(outFilename, outDataList, fitDataList);

  }
  std::cout << "Hey good job, you requested enough memory :)\n";
 
  return 0;
}
