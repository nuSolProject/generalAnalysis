#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <cstdint>
#include <vector>
#include <iomanip>
#include <cmath>
#include <map>
#include <filesystem>
#include <sstream>


struct DataWithUnit {
  double value;       // 64-bit float
  uint32_t magnitude; // 32-bit int representing the magnitude index
  // Next 0x1c(=28) bytes describe the unit in a complex structure
  // For simplicity, we store them as raw bytes. You would parse them as needed.
  uint8_t unit_data[0x1c];  
};

// A helper template to read a single value from the file at a given offset.
template<typename T>
T readValue(std::ifstream &ifs, std::streamoff offset) {
  T val{};
  ifs.seekg(offset, std::ios::beg);
  if (!ifs.read(reinterpret_cast<char*>(&val), sizeof(T))) {
    throw std::runtime_error("Failed to read data at offset " + std::to_string(offset));
  }
  return val;
}

// A helper to read a DataWithUnit structure.
DataWithUnit readDataWithUnit(std::ifstream &ifs, std::streamoff offset) {
  DataWithUnit dwu;
  ifs.seekg(offset, std::ios::beg);
  if (!ifs.read(reinterpret_cast<char*>(&dwu.value), sizeof(dwu.value))) {
    throw std::runtime_error("Failed to read DataWithUnit.value");
  }
  if (!ifs.read(reinterpret_cast<char*>(&dwu.magnitude), sizeof(dwu.magnitude))) {
    throw std::runtime_error("Failed to read DataWithUnit.magnitude");
  }
  if (!ifs.read(reinterpret_cast<char*>(dwu.unit_data), sizeof(dwu.unit_data))) {
    throw std::runtime_error("Failed to read DataWithUnit.unit_data");
  }
  return dwu;
}

// Helper to read an array of 32-bit integers
std::vector<int32_t> readInt32Array(std::ifstream &ifs, std::streamoff offset, size_t count) {
  std::vector<int32_t> arr(count);
  ifs.seekg(offset, std::ios::beg);
  if (!ifs.read(reinterpret_cast<char*>(arr.data()), count * sizeof(int32_t))) {
    throw std::runtime_error("Failed to read int32 array at offset " + std::to_string(offset));
  }
  return arr;
}

// Helper to read an array of 64-bit floats
std::vector<double> readDoubleArray(std::ifstream &ifs, std::streamoff offset, size_t count) {
  std::vector<double> arr(count);
  ifs.seekg(offset, std::ios::beg);
  if (!ifs.read(reinterpret_cast<char*>(arr.data()), count * sizeof(double))) {
    throw std::runtime_error("Failed to read double array at offset " + std::to_string(offset));
  }
  return arr;
}

struct Header {
  uint32_t version;               // 0x00-0x03
  uint32_t data_offset_byte;      // 0x04-0x07
  int32_t ch_on[4];               // ch1_on, ch2_on, ch3_on, ch4_on at 0x08-0x17
  DataWithUnit ch_volt_div_val[4];  // 0x18-0x3f (CH1), 0x40-0x67 (CH2), 0x68-0x8f (CH3), 0x90-0xb7 (CH4)
  DataWithUnit ch_vert_offset[4];   // similarly at 0xb8-0xdf, 0xe0-0x107, 0x108-0x12f, 0x130-0x157
  int32_t digital_on;             // 0x158-0x15b
  std::vector<int32_t> d0_d15_on; // 0x15c-0x19b (16 ints)
  DataWithUnit time_div;          // 0x19c-0x1c3
  DataWithUnit time_delay;        // 0x1c4-0x1eb
  uint32_t wave_length;           // 0x1ec-0x1ef (analog wave length)
  DataWithUnit sample_rate;       // 0x1f0-0x217
  uint32_t digital_wave_length;   // 0x218-0x21b
  DataWithUnit digital_sample_rate; // 0x21c-0x243
  double ch_probe[4];             // 0x244-0x263 (64-bit floats each)
  uint8_t data_width;             // 0x264
  uint8_t byte_order;             // 0x265
  int32_t hori_div_num;           // 0x26c-0x26f
  int32_t ch_vert_code_per_div[4]; // 0x270-0x27f
  int32_t math_switch[4];         // 0x280-0x28f
  DataWithUnit math_vdiv_val[4];  // math vdiv values 0x290-0x32f
  DataWithUnit math_vpos_val[4];  // math vpos values 0x330-0x3cf
  uint32_t math_store_len[4];     // 0x3d0-0x3df
  double math_f_time[4];          // 0x3e0-0x3ff (64-bit doubles)
  int32_t math_vert_code_per_div; // 0x400-0x403

  // Due to the vast number of parameters, we won't implement all fields here.
  // This should demonstrate the approach. You can follow the same pattern for the rest.
};

// A function to read the header from the file
Header readHeader(std::ifstream &ifs) {
  Header h{};

  h.version = readValue<uint32_t>(ifs, 0x00);
  h.data_offset_byte = readValue<uint32_t>(ifs, 0x04);
  for (int i = 0; i < 4; i++) {
    h.ch_on[i] = readValue<int32_t>(ifs, 0x08 + i*4);
  }
  // CH volt/div values
  h.ch_volt_div_val[0] = readDataWithUnit(ifs, 0x18);
  h.ch_volt_div_val[1] = readDataWithUnit(ifs, 0x40);
  h.ch_volt_div_val[2] = readDataWithUnit(ifs, 0x68);
  h.ch_volt_div_val[3] = readDataWithUnit(ifs, 0x90);

  // CH vert offset
  h.ch_vert_offset[0] = readDataWithUnit(ifs, 0xb8);
  h.ch_vert_offset[1] = readDataWithUnit(ifs, 0xe0);
  h.ch_vert_offset[2] = readDataWithUnit(ifs, 0x108);
  h.ch_vert_offset[3] = readDataWithUnit(ifs, 0x130);

  h.digital_on = readValue<int32_t>(ifs, 0x158);
  h.d0_d15_on = readInt32Array(ifs, 0x15c, 16);

  h.time_div = readDataWithUnit(ifs, 0x19c);
  h.time_delay = readDataWithUnit(ifs, 0x1c4);

  h.wave_length = readValue<uint32_t>(ifs, 0x1ec);
  h.sample_rate = readDataWithUnit(ifs, 0x1f0);

  h.digital_wave_length = readValue<uint32_t>(ifs, 0x218);
  h.digital_sample_rate = readDataWithUnit(ifs, 0x21c);

  // Probes
  for (int i = 0; i < 4; i++) {
    h.ch_probe[i] = readValue<double>(ifs, 0x244 + i*8);
  }

  h.data_width = readValue<uint8_t>(ifs, 0x264);
  h.byte_order = readValue<uint8_t>(ifs, 0x265);

  h.hori_div_num = readValue<int32_t>(ifs, 0x26c);

  // ch_vert_code_per_div
  for (int i = 0; i < 4; i++) {
    h.ch_vert_code_per_div[i] = readValue<int32_t>(ifs, 0x270 + i*4);
  }

  // math_switch
  for (int i = 0; i < 4; i++) {
    h.math_switch[i] = readValue<int32_t>(ifs, 0x280 + i*4);
  }

  // math vdiv val
  h.math_vdiv_val[0] = readDataWithUnit(ifs, 0x290);
  h.math_vdiv_val[1] = readDataWithUnit(ifs, 0x2b8);
  h.math_vdiv_val[2] = readDataWithUnit(ifs, 0x2e0);
  h.math_vdiv_val[3] = readDataWithUnit(ifs, 0x308);

  // math vpos val
  h.math_vpos_val[0] = readDataWithUnit(ifs, 0x330);
  h.math_vpos_val[1] = readDataWithUnit(ifs, 0x358);
  h.math_vpos_val[2] = readDataWithUnit(ifs, 0x380);
  h.math_vpos_val[3] = readDataWithUnit(ifs, 0x3a8);

  for (int i = 0; i < 4; i++) {
    h.math_store_len[i] = readValue<uint32_t>(ifs, 0x3d0 + i*4);
  }

  for (int i = 0; i < 4; i++) {
    h.math_f_time[i] = readValue<double>(ifs, 0x3e0 + i*8);
  }

  h.math_vert_code_per_div = readValue<int32_t>(ifs, 0x400);

  // Additional fields can be read similarly...

  return h;
}

// A function to write the header fields to a CSV file.
void writeHeaderToCSV(const Header &h, const std::string &outputPath) {
  std::ofstream ofs(outputPath);
  if (!ofs.is_open()) {
    throw std::runtime_error("Failed to open output CSV file");
  }

  // Write header line
  ofs << "Parameter,Value\n";

  ofs << "version," << h.version << "\n";
  ofs << "data_offset_byte," << h.data_offset_byte << "\n";
  for (int i = 0; i < 4; i++) {
    ofs << "ch" << (i+1) << "_on," << h.ch_on[i] << "\n";
  }

  // For DataWithUnit fields, just print the raw value and magnitude for now
  for (int i = 0; i < 4; i++) {
    ofs << "ch" << (i+1) << "_volt_div_val.value," << h.ch_volt_div_val[i].value << "\n";
    ofs << "ch" << (i+1) << "_volt_div_val.magnitude," << h.ch_volt_div_val[i].magnitude << "\n";
  }

  for (int i = 0; i < 4; i++) {
    ofs << "ch" << (i+1) << "_vert_offset.value," << h.ch_vert_offset[i].value << "\n";
    ofs << "ch" << (i+1) << "_vert_offset.magnitude," << h.ch_vert_offset[i].magnitude << "\n";
  }

  ofs << "digital_on," << h.digital_on << "\n";
  for (int i = 0; i < (int)h.d0_d15_on.size(); i++) {
    ofs << "d" << i << "_on," << h.d0_d15_on[i] << "\n";
  }

  ofs << "time_div.value," << h.time_div.value << "\n";
  ofs << "time_div.magnitude," << h.time_div.magnitude << "\n";
  ofs << "time_delay.value," << h.time_delay.value << "\n";
  ofs << "time_delay.magnitude," << h.time_delay.magnitude << "\n";

  ofs << "wave_length," << h.wave_length << "\n";
  ofs << "sample_rate.value," << h.sample_rate.value << "\n";
  ofs << "sample_rate.magnitude," << h.sample_rate.magnitude << "\n";
  ofs << "digital_wave_length," << h.digital_wave_length << "\n";
  ofs << "digital_sample_rate.value," << h.digital_sample_rate.value << "\n";
  ofs << "digital_sample_rate.magnitude," << h.digital_sample_rate.magnitude << "\n";

  for (int i = 0; i < 4; i++) {
    ofs << "ch" << (i+1) << "_probe," << h.ch_probe[i] << "\n";
  }

  ofs << "data_width," << (int)h.data_width << "\n";
  ofs << "byte_order," << (int)h.byte_order << "\n";
  ofs << "hori_div_num," << h.hori_div_num << "\n";
  for (int i = 0; i < 4; i++) {
    ofs << "ch" << (i+1) << "_vert_code_per_div," << h.ch_vert_code_per_div[i] << "\n";
  }

  for (int i = 0; i < 4; i++) {
    ofs << "math" << (i+1) << "_switch," << h.math_switch[i] << "\n";
  }

  for (int i = 0; i < 4; i++) {
    ofs << "math" << (i+1) << "_vdiv_val.value," << h.math_vdiv_val[i].value << "\n";
    ofs << "math" << (i+1) << "_vdiv_val.magnitude," << h.math_vdiv_val[i].magnitude << "\n";
  }

  for (int i = 0; i < 4; i++) {
    ofs << "math" << (i+1) << "_vpos_val.value," << h.math_vpos_val[i].value << "\n";
    ofs << "math" << (i+1) << "_vpos_val.magnitude," << h.math_vpos_val[i].magnitude << "\n";
  }

  for (int i = 0; i < 4; i++) {
    ofs << "math" << (i+1) << "_store_len," << h.math_store_len[i] << "\n";
  }

  for (int i = 0; i < 4; i++) {
    ofs << "math" << (i+1) << "_f_time," << h.math_f_time[i] << "\n";
  }

  ofs << "math_vert_code_per_div," << h.math_vert_code_per_div << "\n";

  // You can continue printing more parameters similarly.

  ofs.close();
}

// A helper function to get the scaling factor from the magnitude index.
double magnitudeFactor(uint32_t magnitudeIndex) {
  // Table from the specification:
  // 0 YOCTO (1e-24)
  // 1 ZEPTO (1e-21)
  // 2 ATTO  (1e-18)
  // 3 FEMTO (1e-15)
  // 4 PICO  (1e-12)
  // 5 NANO  (1e-9)
  // 6 MICRO (1e-6)
  // 7 MILLI (1e-3)
  // 8 IU (assume 1 for no scaling)
  // 9 KILO  (1e3)
  // 10 MEGA (1e6)
  // 11 GIGA (1e9)
  // 12 TERA (1e12)
  // 13 PETA (1e15)
  // 14 EXA  (1e18)
  // 15 ZETTA(1e21)
  // 16 YOTTA(1e24)

  static const double factors[] = {
    1e-24, 1e-21, 1e-18, 1e-15, 1e-12, 1e-9, 1e-6, 1e-3, 1.0,
    1e3, 1e6, 1e9, 1e12, 1e15, 1e18, 1e21, 1e24
  };

  if (magnitudeIndex > 16) {
    // Unknown magnitude, assume no scaling
    return 1.0;
  }
  return factors[magnitudeIndex];
}

// Convert a DataWithUnit to a properly scaled SI value.
// In a real scenario, you'd parse unit_data to determine if it's V, s, etc.
// For this example, assume the "value" is already given in base units and just apply magnitude.
// If you know the parameter represents time, you get seconds. For voltage, you get volts.
double convertDataWithUnitToSI(const DataWithUnit &dwu) {
  double factor = magnitudeFactor(dwu.magnitude);
  double scaledValue = dwu.value * factor;
  return scaledValue;
}


// This function reads one channel from a binary file given by `h` (the header) and returns both
// the time and voltage arrays. `channelIndex` selects which channel (0-based: 0=CH1,1=CH2,etc.)
std::pair<std::vector<double>, std::vector<double>> decodeChannelData(
								      std::ifstream &ifs, const Header &h, int channelIndex)
{
  // Move file pointer to data start
  ifs.seekg(h.data_offset_byte, std::ios::beg);

  // Determine which channels are on and what order data is stored.
  // The file format suggests that all enabled channels are stored sequentially:
  // CH1 data block first (if on), then CH2 (if on), etc.
  // Math and digital channels follow after analog channels.
  // You must determine the offset for the selected channel from the order and h.wave_length.

  // For simplicity, assume channels are in order CH1, CH2, CH3, CH4.
  // Skip the data for channels before the selected one.
  int analogChannelsOnCount = 0;
  int channelPositionInFile = 0;
  for (int i = 0; i < 4; i++) {
    if (h.ch_on[i] == 1) {
      if (i == channelIndex) {
	channelPositionInFile = analogChannelsOnCount;
      }
      analogChannelsOnCount++;
    }
  }

  if (h.ch_on[channelIndex] != 1) {
    // The requested channel is not on in this file.
    // Return empty vectors or throw an exception.
    throw std::runtime_error("Requested channel is not enabled in this file.");
  }

  // Each enabled analog channel has `h.wave_length` samples.
  // Each sample is either 1 byte (if data_width=0) or 2 bytes (if data_width=1).
  size_t sampleSize = (h.data_width == 0) ? 1 : 2;
  std::streamoff channelDataOffset = h.data_offset_byte + channelPositionInFile * (h.wave_length * sampleSize);

  ifs.seekg(channelDataOffset, std::ios::beg);

  // Convert units for channel parameters
  double ch_volt_div_val = convertDataWithUnitToSI(h.ch_volt_div_val[channelIndex]);
  double ch_vert_offset_val = convertDataWithUnitToSI(h.ch_vert_offset[channelIndex]);
  int code_per_div = h.ch_vert_code_per_div[channelIndex];

  // Determine center_code based on data width
  int center_code = (h.data_width == 0) ? 128 : 32768;

  // Time parameters
  double time_div = convertDataWithUnitToSI(h.time_div);           // in seconds
  double time_delay = convertDataWithUnitToSI(h.time_delay);       // in seconds
  double sample_rate = convertDataWithUnitToSI(h.sample_rate);     // samples/second
  double sample_interval = 1.0 / sample_rate;

  int grid = h.hori_div_num;
  if (grid <= 0) grid = 10; // default to 10 if invalid

  std::vector<double> times(h.wave_length);
  std::vector<double> voltages(h.wave_length);

  // Read samples
  for (uint32_t i = 0; i < h.wave_length; i++) {
    int raw_code = 0;

    if (h.data_width == 0) { 
      // 8-bit data: unsigned char
      uint8_t sample = 0;
      ifs.read(reinterpret_cast<char*>(&sample), 1);
      raw_code = static_cast<int>(sample);
    } else if (h.data_width == 1) {
      // 16-bit data: use uint16_t
      uint16_t sample = 0;
      ifs.read(reinterpret_cast<char*>(&sample), 2);
      raw_code = static_cast<int>(sample);
    }

    int center_code = (h.data_width == 0) ? 128 : 32768;
    double voltage = (raw_code - center_code) * (ch_volt_div_val / code_per_div) - ch_vert_offset_val;
    double time_val = -(time_div * grid / 2.0) - time_delay + i * sample_interval;

    voltages[i] = voltage;
    times[i] = time_val;
  }

  return {times, voltages};
}


int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <inputDir>\n";
    return 1;
  }

  std::string inputDir = argv[1];
  std::string inputPathCh1;
  std::string inputPathCh2;
  std::string outputPath = "output";

  int maxFiles = 11000; //1e5
  for(int file = 0; file < maxFiles; ++file){

    std::ostringstream file1OSS;
    file1OSS << inputDir << "/C1autosave";
    file1OSS << std::setfill('0') << std::setw(5) << file;
    file1OSS << ".bin";
    inputPathCh1 = file1OSS.str();

    //std::cout << inputPathCh1 << "\n";

    std::ostringstream file2OSS;
    file2OSS << inputDir << "/C2autosave";
    file2OSS << std::setfill('0') << std::setw(5) << file;
    file2OSS << ".bin";
    inputPathCh2 = file2OSS.str();

  try {
    // Read CH1 file
    std::ifstream ifsCh1(inputPathCh1, std::ios::binary);
    if (!ifsCh1.is_open()) {
      throw std::runtime_error("Failed to open CH1 input binary file.");
    }
    Header h1 = readHeader(ifsCh1);

    // Read CH2 file
    std::ifstream ifsCh2(inputPathCh2, std::ios::binary);
    if (!ifsCh2.is_open()) {
      throw std::runtime_error("Failed to open CH2 input binary file.");
    }
    Header h2 = readHeader(ifsCh2);

    // For simplicity, assume both files have the same waveform length and timing
    if (h1.wave_length != h2.wave_length) {
      throw std::runtime_error("CH1 and CH2 files have different wave lengths. Cannot combine easily.");
    }

    // Decode CH1 data (channelIndex=0)
    auto [times, ch1_voltages] = decodeChannelData(ifsCh1, h1, 0);

    // Decode CH2 data (channelIndex=1)
    auto [times_ch2, ch2_voltages] = decodeChannelData(ifsCh2, h2, 1);

    // Since we assumed the files are compatible, times_ch2 should match times.
    // Otherwise, you'd need to handle differences (e.g., interpolate or discard extra samples).

    // Use parameters from CH1's header (h1)
    double horizontal_scale = convertDataWithUnitToSI(h1.time_div);    // seconds/div
    double horizontal_delay = convertDataWithUnitToSI(h1.time_delay);  // seconds
    double sample_rate = convertDataWithUnitToSI(h1.sample_rate);      // samples/second
    double sample_interval = 1.0 / sample_rate;
    uint32_t record_length = h1.wave_length;

    // For vertical units and scale:
    double ch1_volt_div = convertDataWithUnitToSI(h1.ch_volt_div_val[0]);
    double ch1_vert_offset = convertDataWithUnitToSI(h1.ch_vert_offset[0]);
    double ch2_volt_div = convertDataWithUnitToSI(h2.ch_volt_div_val[1]);
    double ch2_vert_offset = convertDataWithUnitToSI(h2.ch_vert_offset[1]);

    // Vertical position (in divisions) = offset / volts_per_div
    double ch1_vert_pos = ch1_volt_div != 0.0 ? (ch1_vert_offset / ch1_volt_div) : 0.0;
    double ch2_vert_pos = ch2_volt_div != 0.0 ? (ch2_vert_offset / ch2_volt_div) : 0.0;

    // Probe Attenuation:
    double ch1_probe = h1.ch_probe[0]; // For CH1
    double ch2_probe = h2.ch_probe[1]; // For CH2

    int file_count = 0;
    for(const auto& entry : std::filesystem::directory_iterator(outputPath)){
      if(std::filesystem::is_regular_file(entry.status())){
	++file_count;
      }
    }

    std::string outFilename = outputPath + "/SDS"
      + std::to_string(file_count) + "ALL.csv";

    // Write combined CSV: Time, CH1, CH2
    std::ofstream ofs(outFilename);
    if (!ofs.is_open()) {
      throw std::runtime_error("Failed to open output CSV file.");
    }

    ofs << "Model,SDS6054A\n";
    ofs << "Firmware Version,IDK\n";
    ofs << "Binary Decoded into CSV\n";  // Blank line

    ofs << "Waveform Type,ANALOG,\n";
    ofs << "Point Format,Y,\n";
    ofs << "Horizontal Units,s,\n";
    ofs << "Horizontal Scale," << std::scientific << horizontal_scale << ",\n";
    ofs << "Horizontal Delay," << horizontal_delay << ",\n";
    ofs << "Sample Interval," << sample_interval << ",\n";
    ofs << "Record Length," << record_length << ",\n";
    ofs << "Gating,0.1% to 100.0%,\n";
    ofs << "Probe Attenuation," << ch1_probe << "," << ch2_probe << "\n";
    ofs << "Vertical Units,V,V\n";
    ofs << "Vertical Offset," << ch1_vert_offset << "," << ch2_vert_offset << "\n";
    ofs << "Vertical Scale," << ch1_volt_div << "," << ch2_volt_div << "\n";
    ofs << "Vertical Position," << ch1_vert_pos << "," << ch2_vert_pos << "\n";
    ofs << ",,\n";
    ofs << "bozo\n";
    ofs << ",,\n";
    ofs << "Label,GAGG,Veto\n";
    
    ofs << "Time(s),CH1(V),CH2(V)\n";
    for (size_t i = 0; i < times.size(); i++) {
      ofs << std::setprecision(4) << times[i] << "," << ch1_voltages[i] << "," << ch2_voltages[i] << "\n";
    }
    ofs << "\n";

    std::cout << "Decoding complete. Results written to " << outFilename << "\n";
  } catch (const std::exception &ex) {
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
  }

  }

  return 0;
}
