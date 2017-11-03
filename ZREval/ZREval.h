//
// Created by brad on 3/9/16.
//

#pragma once

#include "Preprocessing.h"
#include "DataTypes.h"
#include "ComparisonMethod.h"

class ZREval {
public:

  ZREval(int lang);
  ~ZREval();

  void runEvalFromCfg(std::string path);

private:

  int overlapRange;

  //does amplitude envelope segementation
  //does MFCC feature extraction
  Preprocessor preprocessor;

  //subsequence comparison functions
  Comparison DPNgramComparison;

  //initial segmentation
  std::vector<segment> prepocessed;

  //clustered discovered subsequences
  std::vector<Cluster> classes;

  //preprocessed - classes
  std::vector<Class> singletons;

  std::vector<std::vector<std::vector<Match>>> findSubSequences(std::vector<segment> segments);

  //clustering
  void generateClusters(std::vector<std::vector<std::vector<Match>>>& pairs);
  std::vector<Cluster> groupByFinfo1(std::vector<std::vector<Match>>& matches);

  //complete transcription
  void generateTranscription();
  void generateSingletons();

  //mainly used for debugging /testing
  void printSegments(std::vector<segment> segments);
  void printMatches(std::vector<std::vector<std::vector<Match>>> matches);
};
