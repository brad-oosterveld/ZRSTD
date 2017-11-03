//
// Created by brad on 3/11/16.
//

#pragma once

#include <vector>
#include <string>
#include "DataTypes.h"

struct qmCell {
  //quality
  double q;
  //previous i and j
  int pi;
  int pj;

  qmCell() {
    q = 0;
    pi = 0;
    pj = 0;
  }
};

//DPNgrams comparison
class Comparison {
public:
  Comparison();
  Comparison( const double& _qT, const double& _btminF, const int& _maxL,
	      const int& _minl, const double& _bon, const double& _penal, const double& _btm )
    : qThreshold( _qT ),
    btminFactor(_btminF),
    maxLines( _maxL ),
    minlen( _minl ),
    bonus( _bon ),
    penalty( _penal ),
    btmin( _btm )
    {
      //REV: Nothing to do.
    }
  
  Comparison(bool vis);
  ~Comparison();

  void loadConfig(std::string);

  std::vector<Match> compare(segment novel, segment base, bool useSmooth);
  //std::vector<match> compareNGram(std::vector<std::vector<double>> novel,
  //                                std::vector<std::vector<double>> base,
  //                                int novelIndex, int baseIndex);

private:
  bool visualize;

  //calculates distance measure between two points
  double cosineDistance(std::vector<double> p1, std::vector<double> p2);
  //creates quality matrix for two sub-sequences
  std::vector<std::vector<qmCell>> getQuality(std::vector<std::vector<double>> novel,
                                              std::vector<std::vector<double>> base);

  //visualization
  void visualizeLine(std::vector<std::vector<qmCell>> quality,
                     std::vector<std::vector<std::vector<int>>> alignments,
                     std::string fname
  );
  void visualizeD(std::vector<std::vector<double>> distance);


  double qThreshold = 19;
  double btminFactor = 0.05;
  int maxLines = 10;
  int minlen = 3;

  double bonus = 1.0;
  double penalty = -1;
  double btmin = 0.1;
};


