//
// Created by brad on 3/9/16.
//

#pragma once

#include "MFCCs.h"
#include <utility>
#include "DataTypes.h"

class Preprocessor {
public:
  Preprocessor(int lang);
  ~Preprocessor();

  //TODO: might want to consider a different vector for each speaker?
  std::vector<segment> preProcessDirectory(std::string path);

private:

  int language;

//store to an internal data structure or
  std::vector<std::pair<finfo, std::vector<double>>> loadFiles(std::string path);

  //segment file
  std::vector<segment> segmentFile(finfo name, std::vector<double> data);

  //extract features
  void extractFeatures(segment& data);

  MFCCs raw;
  MFCCs smooth;

  //util
  std::vector<double> quantizePCM(std::vector<char> raw);
};
