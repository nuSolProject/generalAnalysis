#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <vector>
#include <filesystem>

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

//globals
//input csv data
const int dataSize = 2000;
double* data = nullptr;
double timeData[dataSize];
double dataCh1[dataSize];
double dataCh2[dataSize];
double dataChSum[dataSize];
bool tekScope = 1; //tek scope by default

//outData
struct outData{
  double peak;
  double fullInt;
  double RAPeak; //rolling average peak
  double landauPeak;
  double intPeak;
  double peakTime;
  double decayTime;
  outData(double P, double FI, double RAP, double lP, double iP, double pT,
	  double dT) : peak(P), fullInt(FI), RAPeak(RAP), landauPeak(lP),
		       intPeak(iP), peakTime(pT), decayTime(dT) {}
};

std::vector<outData> outDataList;

//options
std::string inputDirectory;
std::string outputFilename;
bool savePicture; //1 for saving pngs, 0 for no pngs saved
bool pTube; // 0 for SiPM (positive landau), 1 for Ptube (negative)


void updateProgressBar(int currentValue, int maxValue){


  //you can change these
  int maxLength = 30; //# of subdivisions of bar
  std::string programName = "dataAnalysisizer.cc"; //optional

  //you can't (shouldn't change these)
  ++currentValue; //0 based index or < in for loops used to run this
  double currentLength;
  if(currentValue > 0){ //div by 0 check
    currentLength = ((double)currentValue / (double)maxValue) * (double)maxLength;
  }
  else{
    currentLength = 0;
  }

  //if you need to update progress bar
  if(currentLength != (((double)currentValue-1)/(double)maxValue)*(double)maxLength){
  
    int percentage = (currentLength * (100 / (double)maxLength));
 
    std::cout << "\r" <<  "Running " << programName << " |";
    for(int i = 0; i < (int)currentLength; ++i){
      std::cout << '#';
    }
    for(int i = (maxLength - (int)currentLength); i > 0; --i){
      std::cout << ' ';
    }
    if(percentage != 69){
      std::cout << "| @ " << percentage << "%  " << std::flush;
    }
    else{
      std::cout << "| @ " << "nice" <<  "% " << std::flush;  
    }

    if(currentValue == maxValue){ //new line for end of program
      std::cout << std::endl;
    }

  }

 
}


bool readOptions(){

  
  //read in peakIn.txt and get input direc and output name
  std::ifstream inputFile("options.txt");
  
  if(!inputFile){
    std::cout << "Error: Unable to open options.txt file.\n";
    return 0;
  }
  
  std::string line;
  while(std::getline(inputFile, line)){
    
    std::istringstream iss(line);
    std::string parameter;
    std::string value;
    
    if(std::getline(iss, parameter, ':') && std::getline(iss, value)){
      // Trim leading and trailing whitespace from the value
      value.erase(0, value.find_first_not_of(" \t"));
      value.erase(value.find_last_not_of(" \t") + 1);
      if(parameter == "Input (directory)"){
	inputDirectory = value;
      }
      else if(parameter == "Output (file)"){
	outputFilename = value;
      }
      else if(parameter == "Save pictures (bool)"){
	savePicture = (value == "1" || value == "true"); // Set true if '1' or 'true'
      }
      else if(parameter == "Invert signal (bool)"){
	pTube = (value == "1" || value == "true"); // Set true if '1' or 'true'
      }
    }
  }
  
  inputFile.close();
  if(inputDirectory.empty() || outputFilename.empty()){
    std::cout << "Error: Invalid input in options.txt file.\n";

    
    return 0;
  }

  
  return 1;
}


bool getCSVData(std::string filename){

  
  std::ifstream file(filename.c_str());
  std::string line;

  if(!file.is_open()){
    std::cout << "Can't open file: " << filename << "\n";
    return false;
  }

  // Skip first 21 lines (Tek)
  for(int i = 0; i < 21; ++i){
    std::getline(file, line);
  }
  
  // Read and store time and data values
  for(int i = 0; i < dataSize && std::getline(file, line); ++i){
    
    std::stringstream ss(line);
    std::string col1, col2, col3;

    if(std::getline(ss, col1, ',') && std::getline(ss, col2, ',')){
      timeData[i] = std::atof(col1.c_str());
      dataCh1[i] = std::atof(col2.c_str());

      // If there's a third column, read it, else set it to 0.
      if(std::getline(ss, col3, ',')){
        dataCh2[i] = std::atof(col3.c_str());
      }
      else{
        dataCh2[i] = 0.0; // No third column
      }

      // Invert if pTube
      if(pTube){
        dataCh1[i] *= -1;
        dataCh2[i] *= -1;
      }
    }
    else{
      // If line doesn't have at least two columns, break or handle error
      timeData[i] = 0.0;
      dataCh1[i] = 0.0;
      dataCh2[i] = 0.0;
    }
  }

  file.close();

  
  return true;
}


void sumChannels(){


  for(int i = 0; i < dataSize; ++i){
    dataChSum[i] = dataCh1[i] + dataCh2[i];
  }

  
}


void coutCSVData(){

  
  for(int i = 0; i < dataSize; ++i){
    std::cout << data[i] << "," << timeData[i] << "\n";
  }
  
  
}


void zeroData(){

  
  double zeroVal = 0.;
  int nMax = 25;
  for(int n = 0; n < nMax; ++n){ //find offset
    zeroVal += data[n];
  }
  zeroVal /= (double) nMax;
  for(int n = 0; n < dataSize; ++n){ //account for offset
    data[n] -= zeroVal;
  }

  
}


void getPeak(){

  
  double peak = data[0];
  for(int csvIndex = 0; csvIndex < dataSize; ++ csvIndex){
    if(data[csvIndex] > peak){
      peak = data[csvIndex];
    }
  }
  
  outDataList.emplace_back(peak, 0, 0, 0, 0, 0, 0);
  
  
}


void getRAPeak(int fileIndex){

  
  double peakTimeEst;
  
  int offset = 5 * dataSize/1000; //set 0 for 1, 1 for 3, 2 for 5 points...
  int nPoints = (2*offset) + 1;
  int peakIndex, decayTimeIndex;
  double RAPeak = data[0];
  
  for(int csvIndex = offset; csvIndex < (dataSize-offset); ++csvIndex){
    
    double rollingAv = 0;
    for(int j = (csvIndex-offset); j <= (csvIndex+offset); ++j){
      rollingAv += data[j];
    }
    
    rollingAv /= nPoints;
    if(rollingAv > RAPeak){
      RAPeak = rollingAv;
      peakIndex = csvIndex;
    }
	
  }
  
  for(int csvIndex = (peakIndex+1); csvIndex < (dataSize-offset); ++csvIndex){
    
    double rollingAv = 0;
    for(int j = (csvIndex-offset); j <= (csvIndex+offset); ++j){
      rollingAv += data[j];
    }
    rollingAv /= nPoints;
    
    if(rollingAv < (RAPeak/2.718281) || (csvIndex == (dataSize-offset-1))){
      decayTimeIndex = csvIndex;
      break;
    }
    
  }
  
  outDataList[fileIndex].RAPeak = RAPeak;
  outDataList[fileIndex].peakTime = timeData[peakIndex];
  outDataList[fileIndex].decayTime = timeData[decayTimeIndex] - timeData[peakIndex];

  
}


void coutOutListData(){

  
  for(int i = 0; i < outDataList.size(); ++i){
    std::cout << i << "," << outDataList[i].peak << "," << outDataList[i].fullInt
	      << "," << outDataList[i].RAPeak << "," << outDataList[i].landauPeak
	      << "," << outDataList[i].intPeak << ","
	      << outDataList[i].peakTime << "," << outDataList[i].decayTime
	      << "\n";
  }
  
  
}


void getLandauPeak(int fileIndex){
  

  gErrorIgnoreLevel = kFatal;
  
  std::string pngFilename = std::to_string(fileIndex) + ".png";
  /*
    std::vector<double> timeVec(timeData, timeData + sizeof(timeData)
    / sizeof(timeData[0]));
    std::vector<double> dataVec(data, data + sizeof(data) / sizeof(data[0]));*/
  std::vector<double> timeVec;
  std::vector<double> dataVec;
  for(int i = 0; i < dataSize; ++i){
    timeVec.push_back(timeData[i]);
    dataVec.push_back(data[i]);
  }
  
  //std::cout << timeVec.size() << "," << dataVec.size() << "\n";
  int j = 0; //find j for initial time of peak
  while((timeVec[j] < 0) && j < dataSize){
    j += 1;
  }
  int k = 0; //find k corresponding to time esitmate, "won't" break 
  while((timeVec[k] < outDataList[fileIndex].decayTime) && (k < dataSize)){
    k +=1;
  }
  
  // Create a ROOT TCanvas to display the fit result
  TCanvas *c1 = new TCanvas("c1", "Landau Fit", 800, 600);
  // Create a TGraph to hold the data points
  TGraph *graph = new TGraph(timeVec.size(), &timeVec[0], &dataVec[0]);
 
  // Create a Landau distribution function with initial parameter values
  TF1 *landauFunc = new TF1("landauFunc","landau",
			    timeVec[j - (10 * dataSize/1000)],
			    timeVec[k + (160 * dataSize/1000)]);
  // Set the minimizer to Minuit (the original one)
  //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit");
  
  // Fit the Landau function to the data
  graph->Fit("landauFunc", "Q", "", timeVec[j - (10 * dataSize/1000)],
	     timeVec[k + (160 * dataSize/1000)]); 
  //change to "R" to see things 0_0
  // Get the fit results
  double amp = landauFunc->GetParameter(0);
  double mpv = landauFunc->GetParameter(1);
  double width = landauFunc->GetParameter(2);
  // Calculate the maximum value of the fit
  double maxFitValue = landauFunc->GetMaximum();
  double timeOfMaxFitValue = landauFunc->GetX(maxFitValue*.95,timeVec.front(),
					      timeVec.back());
  double peakOverE = maxFitValue / TMath::E();
  double timeOfPeakOverE = landauFunc->GetX(peakOverE, timeOfMaxFitValue,
					    timeVec.back());
  double timeDiff = timeOfPeakOverE - timeOfMaxFitValue;
  
  if(savePicture){
    // Draw the data points and the fit function
    graph->SetTitle("Landau Fit");
    graph->GetXaxis()->SetTitle("Time");
    graph->GetYaxis()->SetTitle("Amplitude");
    graph->Draw("AP");
    landauFunc->Draw("SAME");
    // Save the canvas as a PDF
    c1->SaveAs(pngFilename.c_str());
  }
  
  // Clean up
  delete c1;
  delete graph;
  delete landauFunc;
  
  outDataList[fileIndex].landauPeak = maxFitValue;
  outDataList[fileIndex].peakTime = timeOfMaxFitValue;
  outDataList[fileIndex].decayTime = timeDiff;

  
}


void getIntPeak(int fileIndex){

  
  int intEndIndex;
  double intTime = outDataList[fileIndex].peakTime
    + (3 * outDataList[fileIndex].decayTime);
       
  for(int csvIndex = 0; csvIndex < dataSize; ++csvIndex){
    
    if(timeData[csvIndex] > intTime || csvIndex == (dataSize -1)){
      intEndIndex = csvIndex;
      break;
    }
    
  }
  
  double intPeak = 0;
  for(int csvIndex = 0; csvIndex < intEndIndex; ++csvIndex){
    intPeak += .5 * (data[csvIndex] + data[csvIndex+1])
      * (timeData[csvIndex+1] - timeData[csvIndex]);
  }
  
  outDataList[fileIndex].intPeak = intPeak;
  
  
}


void getFullInt(int fileIndex){

  
  double fullInt = 0;
  for(int csvIndex = 0; csvIndex < (dataSize-1); ++csvIndex){
    fullInt += .5 * (data[csvIndex] + data[csvIndex+1])
      * (timeData[csvIndex+1] - timeData[csvIndex]);
  }
  
  outDataList[fileIndex].fullInt = fullInt;

  
}

void storeDataInNtuple(const std::vector<outData>& dataList){

  
  // Create a TNtupleD
  TNtupleD* ntuple = new TNtupleD("ntuple", "Data Ntuple", "peak:fullInt:RAPeak:landauPeak:intPeak:peakTime:decayTime");
  // Fill the ntuple
  for(const auto& data : dataList){
    Double_t values[] = {data.peak, data.fullInt, data.RAPeak, data.landauPeak,
      data.intPeak, data.peakTime, data.decayTime};
    ntuple->Fill(values);
  }
  // Save the ntuple to a ROOT file using the global outputFilename
  TFile* outFile = new TFile(outputFilename.c_str(), "recreate"); //recreate
  ntuple->Write();
  outFile->Close(); // Ensure we properly close and release the file
  // Clean up
  delete outFile; // Make sure to delete the TFile object to free memory
  delete ntuple; // Also delete the TNtupleD to prevent memory leaks

  
}


int main(){

  
  if(!readOptions()){ //if invalid options.txt
    return 0;
  }

  int file_count = 0;
  for(const auto& entry : std::filesystem::directory_iterator(inputDirectory)){
    if(std::filesystem::is_regular_file(entry.status())){
      ++file_count;
    }
  }
  //std::cout << file_count << "\n";

  data = dataCh2;
  
  std::string currentCSVFileName;
  int maxCSVFiles = 1e6;//1e6

  if(maxCSVFiles < file_count){
    file_count = maxCSVFiles;
  }
  
  for(int fileIndex = 0; fileIndex < maxCSVFiles; ++fileIndex){
    
    std::ostringstream oss;
    if(tekScope){
      oss << inputDirectory << "/tek" << fileIndex << "ALL.csv";
    }
    else{
      oss << inputDirectory << "/SDS" << fileIndex << "ALL.csv";
    }
    currentCSVFileName = oss.str();
    if(!getCSVData(currentCSVFileName)){ //can't read/find filename
      if(!tekScope){
	break;
      }
      else{
	--fileIndex;
	tekScope = 0;
	continue;
      }
    }

    sumChannels();
    
    //coutCSVData();
    //zeroData();
    getPeak();
    //getFullInt(fileIndex);
    getRAPeak(fileIndex);
    getLandauPeak(fileIndex);
    getIntPeak(fileIndex);
    getFullInt(fileIndex);

    updateProgressBar(fileIndex, file_count);

    //std::cout << fileIndex << "\n";
    
  }
  //coutOutListData();
  storeDataInNtuple(outDataList);
  
  
  return 0;
}
