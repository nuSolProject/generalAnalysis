cmake_minimum_required(VERSION 3.10)

project(generalAnalysis CXX)

find_package(ROOT)

add_executable(CSVAnalysisBrooks CSVAnalysisBrooks.cc)
target_link_libraries(CSVAnalysisBrooks PUBLIC ${ROOT_LIBRARIES})

add_executable(SDSbin2TekCSV SDSbin2TekCSV.cc)
target_link_libraries(SDSbin2TekCSV PUBLIC ${ROOT_LIBRARIES})
