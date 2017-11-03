//
// Created by brad on 1/4/16.
//

#pragma once

#include <vector>
#include <string>
//#include <log4cxx/logger.h>

class MFCCs {
public:
  MFCCs();
  ~MFCCs();
  void loadConfig(std::string path);
  std::vector<std::vector<double>> getFeatures(std::vector<double> PCMData);

private:
  //log4cxx::LoggerPtr logger;

  //configurable params:
  bool normVolume; //noramlize volume so loudest value is MXPCM
  int numDiffs; //number of differentials to take
  double smoothFactor;

  int minDataSize;

  std::vector<double> normalizeVolume(std::vector<double> raw);

  std::vector<std::vector<double>> smooth(std::vector<std::vector<double>> rough);
  std::vector<std::vector<double>> calculateMFCCs(std::vector<double> signal);
  std::vector<std::vector<double>> getDiffs(std::vector<std::vector<double>> original);

  static inline int NextPowerOf2(int val);
  inline std::vector<double> ceplifter(double n, double l);
  std::vector<double> filter(std::vector<double> b, double a, std::vector<double> x);
  std::vector<double> preEmphasis(double coeff, std::vector<double> x);
  std::vector<double> window(int windowSize);

};
