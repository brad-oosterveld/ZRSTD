//
// Created by brad on 3/11/16.
//

#pragma once

#include <vector>
#include <string>
#include <map>

struct finfo {
  std::string name;
  std::string file;
  int index;

  //TODO: assignment operators
  //finfo& operator=(finfo& other){
  //  name=other.name;
  //  file=other.file;
  //  index=other.index;
  //  return *this;
  //}

  finfo() :
    name(""),
    file(""),
    index(-1) {}
};

//we can leave this as a single class now and worry about efficiency later if it becomes a problem
//maybe we can intentionally slice some of it off


struct segment {
  finfo fname;
  int start;
  int end;
  std::vector<double> PCM;
  std::vector <std::vector<double>> smooth;
  //std::vector<std::vector<double>> raw;
};

struct Match {
  double quality;
  segment s1;
  segment s2;
  int classification;

  Match() :
    quality(0),
    //start1(0),
    //end1(0),
    //chunk1(),
    //start2(0),
    //end2(0),
    //chunk2(),
    s1(),
    s2(),
    classification(-1) {}
};

struct Cluster {
  Cluster();
  Cluster(Match &seed, int clusterIndex, int range);
  ~Cluster();

  int id;
  int overlapRange;

  //bool is true if connected componets for s2 have been searcehd;
  std::vector <std::pair<bool, Match &>> matches;
  //vector <segment> refined;
  std::map<int, segment> refined;
  std::map<int, std::vector<segment>> exemplars;

  //re averages exemplars of given index
  void refine(int index);

  int s1InRange(segment s);
  //TODO: do we want this to be a reference?
  void addMatch(Match &m);

  void addConnected(std::vector <std::vector<std::vector < Match>>
  > &pairs);

};

struct Class {
//start and end are in frames
  Class() {

  }

  Class(segment seed, bool type) :
    start(seed.start),
    end(seed.end) {
    fromCluster = type;
    fname = seed.fname;
  }

  Class(int s, int e, int type) :
    start(s),
    end(e),
    fromCluster(type) {

  }

  int start;
  int end;
  bool fromCluster;
  finfo fname;
};

struct TranscriptionUnit {
  TranscriptionUnit() :
    subSequences() {
  }

  int segmentID;
  std::vector <Class> subSequences;
  void insertClass(Class newMember);
};