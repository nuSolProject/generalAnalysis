void exp3PulseFit(){

  //globals
  const int nSamples = 1000;
  double timeData[nSamples];
  double dataCh1[nSamples];
  double dataCh2[nSamples];

  //start store csv into RAM
  std::string filename = "/home/bhartsock/Desktop/root/exp3PulseFit/copy/tek0ALL.csv";

  std::ifstream file(filename.c_str());
  std::string line;

  if(!file.is_open()){
    std::cout << "Can't open file: " << filename << "\n";
    return;
  }

  // Skip first 21 lines (Tek)
  for(int i = 0; i < 21; ++i){
    std::getline(file, line);
  }
  
  // Read and store time and data values
  for(int i = 0; i < nSamples && std::getline(file, line); ++i){
    
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

    }
    else{
      // If line doesn't have at least two columns, break or handle error
      timeData[i] = 0.0;
      dataCh1[i] = 0.0;
      dataCh2[i] = 0.0;
    }
  }

  file.close(); //end storing csv into RAM

  //we have data :)
  for(int sample = 0; sample < nSamples; ++sample){
    //std::cout << timeData[sample] << "," << dataCh1[sample] << "," << dataCh2[sample] << "\n";
  }


  //lets start rootifying it :...(
  std::vector<double> timeVec;
  std::vector<double> dataVec;
  for(int i = 0; i < nSamples; ++i){
    timeVec.push_back(timeData[i]);
    dataVec.push_back(dataCh1[i]);
  }

  // Create a ROOT TCanvas to display the fit result
  TCanvas *c1 = new TCanvas("c1", "exp3 Fit", 1800, 1200);
  // Create a TGraph to hold the data points
  TGraph *graph = new TGraph(timeVec.size(), &timeVec[0], &dataVec[0]);
  //create fit function
  TF1 *fit = new TF1("fit","(x>[5]) ? [6] + [0]*(exp(-(x-[5])/[2]) - exp(-(x-[5])/[1])) + [3]*(exp(-(x-[5])/[4]) - exp(-(x-[5])/[1])) : [6]",-2.8e-8, 8.34e-7);

  fit->SetParameter(0, 0.25);      // main component scale
  fit->SetParameter(1, 500e-12);     // rise tau
  fit->SetParameter(2, 250e-9);     // main decay tau

  fit->SetParameter(3, 0.05);      // slow component scale
  fit->SetParameter(4, 500e-9);    // slow decay tau

  fit->SetParameter(5, 0.0);       // t0
  fit->SetParameter(6, 0.0);       // baseline

  fit->SetParNames("A_{f}", "#tau_{r}", "#tau_{d1}", "A_{s}", "#tau_{d2}", "t_{0}", "V_{0}");

  graph->SetTitle("");
  graph->GetXaxis()->SetTitle("Time (s)");
  graph->GetYaxis()->SetTitle("Voltage (V)");

  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(0.5);

  graph->Fit("fit", "Q", "", -3e-8, 8.34e-7);
  graph->Draw("AP");
  fit->Draw("SAME");

  gStyle -> SetOptFit(1);

  c1->SaveAs("exp3FitTek0.png");




  std::cout << "Good job Bozo\n";
}
