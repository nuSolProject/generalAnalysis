//options then run with 'root -l -b -q SDSParserV1.cc'
struct Options {
  std::string inputDirectory = "/homes/h293d863/Public/labData/LiSTest/run1/Autosave"; //where the files are (use 'pwd' on data directory and copy here)
  std::string outputFilename = "outTest.root"; //name of output file
  int channelNo = 3;                   //channel number, check filename C1auto... C2auto...
  bool invertSignal = true;           //set true for PMT, false for SiPM
  bool savePictures = false;           //false by default... USE WITH SMALL maxFiles!!!
  std::string pictureDirectory = "png";
  int headerLines = 21;                //SDS csv format
  double triggerThreshold = 0.001;     //volts after inversion
  int rollingHalfWindow = 10;          //rolling average width in samples
  int baselineSamples = 0;             //
  int maxFiles = 0;                    //0 = process every matching CSV
};

struct Waveform {
  std::vector<double> time;
  std::vector<double> voltage;
};

struct Measurements {
  double peak = std::numeric_limits<double>::quiet_NaN();
  double trigTime = std::numeric_limits<double>::quiet_NaN();
  double fullInt = std::numeric_limits<double>::quiet_NaN();
  double RAPeak = std::numeric_limits<double>::quiet_NaN();
  double intPeak = std::numeric_limits<double>::quiet_NaN();
  double peakTime = std::numeric_limits<double>::quiet_NaN();
  double decayTime = std::numeric_limits<double>::quiet_NaN();
};

//only C<channel>autosave<five-or-more-digits>.csv are included
std::vector<std::pair<long long, std::filesystem::path>> findCSVFiles(const Options& opt) {
  namespace fs = std::filesystem;
  std::vector<std::pair<long long, fs::path>> files;
  const std::string prefix = "C" + std::to_string(opt.channelNo) + "autosave";

  for(const fs::directory_entry& entry : fs::directory_iterator(opt.inputDirectory)) {
    if(!entry.is_regular_file()) continue;
    const std::string name = entry.path().filename().string();
    if(name.size() < prefix.size() + 5 + 4 ||
       name.compare(0, prefix.size(), prefix) != 0 ||
       name.compare(name.size() - 4, 4, ".csv") != 0) continue;

    const std::string digits = name.substr(prefix.size(), name.size() - prefix.size() - 4);
    if(!std::all_of(digits.begin(), digits.end(), [](unsigned char c) {
      return std::isdigit(c) != 0;
    })) continue;

    try {
      files.emplace_back(std::stoll(digits), entry.path());
    }
    catch(const std::exception&) {
      std::cerr << "Skipping filename with invalid sequence number: " << name << '\n';
    }
  }

  std::sort(files.begin(), files.end(),
    [](const std::pair<long long, fs::path>& a,
       const std::pair<long long, fs::path>& b) {
      if(a.first != b.first) return a.first < b.first;
      return a.second.string() < b.second.string();
    });

  if(opt.maxFiles > 0 && files.size() > static_cast<std::size_t>(opt.maxFiles)) {
    files.resize(static_cast<std::size_t>(opt.maxFiles));
  }
  return files;
}

bool parseDouble(const std::string& cell, double& value) {
  try {
    std::size_t consumed = 0;
    value = std::stod(cell, &consumed);
    while(consumed < cell.size() &&
          std::isspace(static_cast<unsigned char>(cell[consumed]))) ++consumed;
    return consumed == cell.size() && std::isfinite(value);
  }
  catch(const std::exception&) {
    return false;
  }
}

//readCSV files... ty ChatGPT for this YEARS ago
bool readCSV(const std::filesystem::path& path, const Options& opt, Waveform& wave) {
  std::ifstream file(path);
  if(!file) return false;

  std::string line;
  for(int i = 0; i < opt.headerLines; ++i) {
    if(!std::getline(file, line)) return false;
  }

  // SDS CSV: time in column 1, voltage for the selected channel in column 2.
  // The channel is selected by the *filename*, not by a third CSV column.
  while(std::getline(file, line)) {
    std::istringstream row(line);
    std::string timeCell, voltageCell;
    if(!std::getline(row, timeCell, ',') || !std::getline(row, voltageCell, ',')) continue;
    double t = 0.0, v = 0.0;
    if(!parseDouble(timeCell, t) || !parseDouble(voltageCell, v)) continue;
    wave.time.push_back(t);
    wave.voltage.push_back(opt.invertSignal ? -v : v);
  }

  if(wave.time.size() < 2) return false;

  if(opt.baselineSamples > 0) {
    const int count = std::min(opt.baselineSamples, static_cast<int>(wave.voltage.size()));
    double baseline = 0.0;
    for(int i = 0; i < count; ++i) baseline += wave.voltage[i];
    baseline /= count;
    for(double& voltage : wave.voltage) voltage -= baseline;
  }
  return true;
}

double integrate(const Waveform& wave, int begin, int end) {
  double integral = 0.0;
  for(int i = begin; i < end; ++i) {
    integral += 0.5 * (wave.voltage[i] + wave.voltage[i + 1]) *
                (wave.time[i + 1] - wave.time[i]);
  }
  return integral;
}

Measurements analyze(const Waveform& wave, const Options& opt) {
  Measurements result;
  const int n = static_cast<int>(wave.voltage.size());

  const int rawPeakIndex = static_cast<int>(
    std::max_element(wave.voltage.begin(), wave.voltage.end()) - wave.voltage.begin());
  result.peak = wave.voltage[rawPeakIndex];
  result.fullInt = integrate(wave, 0, n - 1);

  //last below-threshold sample before the raw peak
  int triggerIndex = -1;
  for(int i = rawPeakIndex; i >= 0; --i) {
    if(wave.voltage[i] < opt.triggerThreshold) {
      triggerIndex = i;
      result.trigTime = wave.time[i];
      break;
    }
  }

  //rolling average peak, followed by the first 1/e crossing
  const int halfWindow = std::min(opt.rollingHalfWindow, (n - 1) / 2);
  const int windowSize = 2 * halfWindow + 1;
  double sum = 0.0;
  for(int i = 0; i < windowSize; ++i) sum += wave.voltage[i];

  int rollingPeakIndex = halfWindow;
  double rollingPeak = sum / windowSize;
  for(int center = halfWindow + 1; center < n - halfWindow; ++center) {
    sum += wave.voltage[center + halfWindow] - wave.voltage[center - halfWindow - 1];
    const double mean = sum / windowSize;
    if(mean > rollingPeak) {
      rollingPeak = mean;
      rollingPeakIndex = center;
    }
  }
  result.RAPeak = rollingPeak;
  result.peakTime = wave.time[rollingPeakIndex];

  //recompute a rolling mean
  //if pulse never decays to 1/e within the waveform, leave decayTime as NaN.
  const double cutoff = rollingPeak / std::exp(1.0);
  int decayIndex = -1;
  for(int center = rollingPeakIndex + 1; center < n - halfWindow; ++center) {
    double mean = 0.0;
    for(int j = center - halfWindow; j <= center + halfWindow; ++j) {
      mean += wave.voltage[j];
    }
    mean /= windowSize;
    if(mean < cutoff) {
      decayIndex = center;
      result.decayTime = wave.time[center] - wave.time[rollingPeakIndex];
      break;
    }
  }

  if(triggerIndex >= 0 && decayIndex >= triggerIndex) {
    result.intPeak = integrate(wave, triggerIndex, decayIndex);
  }
  return result;
}

void saveWaveformPNG(const std::filesystem::path& csvPath,
                     const std::filesystem::path& outputDir,
                     const Waveform& wave) {
  TCanvas canvas("waveformCanvas", "Waveform", 1000, 650);
  TGraph graph(static_cast<int>(wave.time.size()), wave.time.data(), wave.voltage.data());
  const std::string title = csvPath.filename().string() + ";Time (s);Voltage (V)";
  graph.SetTitle(title.c_str());
  graph.SetLineWidth(2);
  graph.Draw("AL");
  canvas.SetGrid();
  canvas.SaveAs((outputDir / (csvPath.stem().string() + ".png")).string().c_str());
}

//for root -l -b -q SDSParserV1.cc
void SDSParserV1() {
  const Options opt;  // Edit the defaults in Options above.

  if(opt.channelNo <= 0 || opt.headerLines < 0 || opt.rollingHalfWindow < 0 ||
     opt.baselineSamples < 0 || opt.maxFiles < 0) {
    std::cerr << "Invalid Options values.\n";
    return;
  }

  gROOT->SetBatch(kTRUE);
  try {
    const std::vector<std::pair<long long, std::filesystem::path>> csvFiles = findCSVFiles(opt);
    if(csvFiles.empty()) {
      std::cerr << "No C" << opt.channelNo << "autosave#####.csv files found in "
                << opt.inputDirectory << '\n';
      return;
    }

    if(opt.savePictures) std::filesystem::create_directories(opt.pictureDirectory);

    TFile output(opt.outputFilename.c_str(), "RECREATE");
    if(output.IsZombie()) {
      std::cerr << "Unable to create ROOT file: " << opt.outputFilename << '\n';
      return;
    }

    TNtupleD ntuple("ntuple", "SDS waveform analysis",
      "fileNo:peak:trigTime:fullInt:RAPeak:intPeak:peakTime:decayTime");
    ntuple.SetDirectory(nullptr);

    int processed = 0;
    for(std::size_t i = 0; i < csvFiles.size(); ++i) {
      const std::filesystem::path& path = csvFiles[i].second;
      Waveform wave;
      if(!readCSV(path, opt, wave)) {
        std::cerr << "\nSkipping unreadable/empty CSV: " << path << '\n';
        continue;
      }

      const Measurements m = analyze(wave, opt);
      Double_t values[] = {
        static_cast<Double_t>(csvFiles[i].first), m.peak, m.trigTime,
        m.fullInt, m.RAPeak, m.intPeak, m.peakTime, m.decayTime
      };
      ntuple.Fill(values);
      ++processed;

      if(opt.savePictures) saveWaveformPNG(path, opt.pictureDirectory, wave);
      const int percent = static_cast<int>(100 * (i + 1) / csvFiles.size());
      std::cout << "\rProcessing SDS CSVs: " << (i + 1) << '/' << csvFiles.size()
                << " (" << percent << "%)" << std::flush;
    }
    std::cout << '\n';

    output.cd();
    ntuple.Write();
    output.Close();
    std::cout << "Saved " << processed << " events to " << opt.outputFilename << '\n';
  }
  catch(const std::exception& e) {
    std::cerr << "Error: " << e.what() << '\n';
  }
}
