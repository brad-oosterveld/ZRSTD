//
// Created by brad on 3/18/16.
//

#include "DataTypes.h"
#include <math.h>
#include <iostream>

using namespace std;

Cluster::Cluster() {

}

//assumes cluster is always seeded with s1
Cluster::Cluster(Match &seed, int clusterIndex, int range) {
  id = clusterIndex;
  overlapRange = range;
  seed.classification = id;
  matches.push_back(std::pair<bool, Match &>(false, seed));
  refined.insert(std::pair<int, segment>(seed.s1.fname.index, seed.s1));
  exemplars.insert(
    std::pair<int, std::vector<segment>>(seed.s1.fname.index, std::vector<segment>(1, seed.s1)));

  refined.insert(std::pair<int, segment>(seed.s2.fname.index, seed.s2));
  exemplars.insert(
    std::pair<int, std::vector<segment>>(seed.s2.fname.index, std::vector<segment>(1, seed.s2)));
}

Cluster::~Cluster() {

}

int Cluster::s1InRange(segment s) {
  for (auto r: refined) {
    if (abs(r.second.start - s.start) < overlapRange && abs(r.second.end - s.end) < overlapRange) {
      //std::cout<<abs(r.second.start - s.start)<<" "<<abs(r.second.end - s.end);
      return r.first;
    }
  }
  return -1;
}

//adds new match and updates refined
void Cluster::addMatch(Match &m) {

  //TODO: make sure bad things aren't happening with classifications being over written

  if (m.classification == -1) {
    m.classification = id;
  }
  else {
    std::cout<<"BAD! match getting reclassified!\n";
  }

  //s1 first
  //if there is no entry for index add a new one, if there is add and refine
  if (exemplars.find(m.s1.fname.index) == exemplars.end()) {
    std::vector<segment> newExemplar(1, m.s1);
    exemplars.insert(std::pair<int, std::vector<segment>>(m.s1.fname.index, newExemplar));
    refined.insert(std::pair<int, segment>(m.s1.fname.index, m.s1));
  } else {
    exemplars[m.s1.fname.index].push_back(m.s1);
    refine(m.s1.fname.index);
  }

  //do the same thing for s2
  if (exemplars.find(m.s2.fname.index) == exemplars.end()) {
    std::vector<segment> newExemplar(1, m.s2);
    exemplars.insert(std::pair<int, std::vector<segment>>(m.s2.fname.index, newExemplar));
    refined.insert(std::pair<int, segment>(m.s2.fname.index, m.s2));
  } else {
    exemplars[m.s2.fname.index].push_back(m.s2);
    refine(m.s2.fname.index);
  }
}

void Cluster::refine(int index) {

  int avgStart = 0;
  int avgEnd = 0;

  for (segment s: exemplars[index]) {
    avgStart += s.start;
    avgEnd += s.end;
  }

  refined[index].start = avgStart / (int) exemplars[index].size();
  refined[index].end = avgEnd / (int) exemplars[index].size();

}

//adds connected components
void Cluster::addConnected(vector<vector<vector<Match>>> &pairs) {

  //TODO: make sure this works with adding new entries before the loop ends
  for (auto &m: matches) {
    if (!m.first) {
      //search through pairs[index of s2] for overlapping segments that might match other segments
      //if they overlap add them
      for (vector<Match> &pair: pairs[m.second.s2.fname.index]) {
        for (Match p: pair) {
          if (p.classification == -1 && s1InRange(p.s1) > -1) addMatch(p);
          //else {
          //  cout<<"tried to add an overlaping match , but it was already classified\n";
          //  cout<<"trying to add to class: "<<id<<" already classified as: "<<p.classification<<
          //  " "<<s1InRange(p.s1)<<endl;
          //}
        }
      }
      //search the other part of the triangle
      /*
       * for every index less than m.second.s2.fname.index
       * search through matches[that index][start + m.second.s2.fname.index- that index]
       */
      //for (int i = 0; i < m.second.s2.fname.index; i++) {
      //  for (Match p:pairs[m.second.s2.fname.index][m.second.s2.fname.index - i] ) {
      //    if (p.classification == -1 && s1InRange(p.s1) > -1) addMatch(p);
      //
      //  }
      //}

      m.first = true;
    }
  }

}