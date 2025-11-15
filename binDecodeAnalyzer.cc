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

std::vector<outData> outDataList;

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
    fullInt *= 0.5 * 8.0 * (samplesPerEvent-1); //1/2 * 8ns * total

    outDataList[eventNo].fullInt1 = fullInt;

    //CH2
    fullInt = 0.0;
    for(int sampleNo = 0; sampleNo < samplesPerEvent-1; ++sampleNo){
      fullInt += RPCHData[eventNo][1][sampleNo];
      fullInt += RPCHData[eventNo][1][sampleNo+1];
    }
    fullInt *= 0.5 * 8.0 * (samplesPerEvent-1); //1/2 * 8ns * total

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
      fullInt *= 0.5 * 8.0 * (totalSamples-1); //1/2 * 8ns * total

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
      fullInt *= 0.5 * 8.0 * (totalSamples-1); //1/2 * 8ns * total

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


void storeDataInNtuple(const std::vector<outData>& dataList){

  
  // Create a TNtupleD
  TNtupleD* ntuple = new TNtupleD("ntuple", "Data Ntuple", "timeStamp:fullInt1:trigTime1:smartInt1:fullInt2:trigTime2:smartInt2");
  // Fill the ntuple
  for(const auto& data : dataList){
    Double_t values[] = {std::stod(data.timeStamp), data.fullInt1, data.trigTime1,
      data.smartInt1, data.fullInt2, data.trigTime2, data.smartInt2};
    ntuple->Fill(values);
  }
  // Save the ntuple to a ROOT file using the global outputFilename
  TFile* outFile = new TFile("137Cs.root", "recreate"); //recreate
  ntuple->Write();
  outFile->Close(); // Ensure we properly close and release the file
  // Clean up
  delete outFile; // Make sure to delete the TFile object to free memory
  delete ntuple; // Also delete the TNtupleD to prevent memory leaks

  
}


int binDecodeAnalyzer(){


  std::string curFilename = "/homes/h293d863/Documents/lab/newAnalysisizer/data/newRP/";
  curFilename += "137Cs.bin";
  
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

  storeDataInNtuple(outDataList);
  std::cout << "Hey good job, you requested enough memory :)\n";
 
  return 0;
}
