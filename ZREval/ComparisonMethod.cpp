//
// Created by brad on 3/11/16.
//

#include "ComparisonMethod.h"
#include <math.h>
#include <fstream>
#include <iostream>

using namespace std;

Comparison::Comparison() {
  visualize = false;
}

Comparison::Comparison(bool vis) {
  visualize = vis;
}

Comparison::~Comparison() {

}

//Comparison ******************************************************************
double Comparison::cosineDistance(vector<double> p1, vector<double> p2) {
  double dot = 0.0;
  double mag1 = 0.0;
  double mag2 = 0.0;
  for (int i = 0; i < p1.size(); i++) {
    dot += p1[i] * p2[i];
    mag1 += p1[i] * p1[i];
    mag2 += p2[i] * p2[i];
  }

  double distance = dot / (sqrt(mag1) * sqrt(mag2));
  return distance;
}

//calculate quality matrix
vector<vector<qmCell> > Comparison::getQuality(vector<vector<double> > novel,
                                               vector<vector<double> > base) {
  //intializing dtw matrix
  int n = novel.size();
  int b = base.size();

  vector<vector<double> > distance(n, vector<double>(b, 0));
  vector<vector<qmCell> > quality(n, vector<qmCell>(b));

  //weights to favor local alignment
  //set by hand for this task
  //double bonus = 1;
  //double penalty = -1;
  //double btmin = .1;
  //calculate dtw distance
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < b; j++) {
      double cost = cosineDistance(novel[i], base[j]);
      distance[i][j] = cost;
    }
  }

#ifdef DEBUG10
  if(distance.size() > 0 && distance[0].size() > 0)
    {
fprintf(stdout, "SEG1: ");
for(size_t x=0; x<7; ++x)
  {
    fprintf(stdout, " [%lf] ", novel[novel.size()-1][x]);
  }
fprintf(stdout, "\nSEG2: ");
for(size_t x=0; x<7; ++x)
  {
    fprintf(stdout, " [%lf] ", base[base.size()-1][x]);
  }
fprintf(stdout, "\n");

fprintf(stdout, "Top-right dist: [%lf] (len1 = [%ld], len2 = [%ld])\n", distance[ distance.size()-1 ][ distance[0].size()-1 ], distance.size(), distance[0].size());
  }
#endif

  for (int i = 0; i < n - 1; i++) {
    for (int j = 0; j < b - 1; j++) {
      if (distance[i + 1][j + 1] > 0) {
        double diag = distance[i][j];
        double horiz = distance[i + 1][j];
        double vert = distance[i][j + 1];

        if (diag >= vert && diag >= horiz) {
          quality[i + 1][j + 1].q = diag * bonus + quality[i][j].q;
          quality[i + 1][j + 1].pi = i;
          quality[i + 1][j + 1].pj = j;
        }
        else if (horiz > vert && horiz > diag) {
          quality[i + 1][j + 1].q =
            (((1 - horiz) * quality[i + 1][j].q) * penalty) + quality[i + 1][j].q;
          quality[i + 1][j + 1].pi = i + 1;
          quality[i + 1][j + 1].pj = j;
        }
        else {
          quality[i + 1][j + 1].q =
            (((1 - vert) * quality[i][j + 1].q) * penalty) + quality[i][j + 1].q;
          quality[i + 1][j + 1].pi = i;
          quality[i + 1][j + 1].pj = j + 1;
        }
        if (quality[i + 1][j + 1].q < btmin) quality[i + 1][j + 1].q = 0;
      }
    }
  }
  visualizeD(distance);
  return quality;
}

//remove neighboring quality scores until they start increasing
void clearAlignment(vector<vector<int>> alignment, vector<vector<qmCell>> &qm) {
  for (vector<int> p: alignment) {
    double origin = qm[p[0]][p[1]].q;
    int after = p[1];
    double next = origin;
    while (after + 1 <= qm[0].size() && next >= qm[p[0]][after + 1].q && next != 0) {
      qm[p[0]][after].q = 0;
      after++;
      next = qm[p[0]][after].q;
    }
    int before = p[1];
    next = origin;
    while (before - 1 > 0 && next >= qm[p[0]][before - 1].q && next != 0) {
      qm[p[0]][before].q = 0;
      before--;
      next = qm[p[0]][before].q;
    }
    qm[p[0]][p[1]].q = -1;
  }
}

//remove neighboring quality scores in a box bounded by the alignment, prevents any overlap
void clearAlignmentBox(vector<vector<int>> alignment, vector<vector<qmCell>> &qm) {
  //assert(alignment.size() >= 1);
  for (int i = alignment[alignment.size() - 1][0]; i <= alignment[0][0]; i++) {
    for (int j = alignment[alignment.size() - 1][1]; j <= alignment[0][1]; j++) {
      qm[i][j].q = -5;
    }
  }
}

//genrateates and traverses quality matrix and returns list of alignment pairs
vector<Match> Comparison::compare(segment novel, segment base, bool useSmooth) {
  //TODO: figure out what to do about raw vs smooth
  //TODO: uughhh a bool flag or something is probably the quickest, but we might want different params

  vector<vector<qmCell>> QualityMatrix;
  if (useSmooth) {
    QualityMatrix = getQuality(novel.smooth, base.smooth);
  }
  else {
  //  QualityMatrix = getQuality(novel.raw, base.raw);
  }

  vector<vector<qmCell>> Original = QualityMatrix;

  int numLines = 0;

  vector<Match> matches;

  //dummy qmCell used to track alignment path
  qmCell best;
  best.q = qThreshold;

  //get starting position of best alignment
  for (int i = 0; i < QualityMatrix.size(); i++) {
    for (int j = 0; j < QualityMatrix[0].size(); j++) {
      if (QualityMatrix[i][j].q > best.q) {
        best.q = QualityMatrix[i][j].q;
        best.pi = i;
        best.pj = j;
      }
    }
  }

  //find alignments above quality threshold
  //faux recursion through qm finding optimal alignments and setting "used" cells to 0

  vector<vector<vector<int>>> alignments;
  int ntries=0;
  while (best.q > qThreshold && numLines < maxLines) {
    ++ntries;
    qmCell nextBest = best;
    vector<vector<int>> alignment;
    while (QualityMatrix[nextBest.pi][nextBest.pj].q > best.q * btminFactor) {
      vector<int> point{nextBest.pi, nextBest.pj};
      alignment.push_back(point);
      //hack: current not previous indicies in nextBest
      nextBest.q = QualityMatrix[nextBest.pi][nextBest.pj].q;
      //go to next cell
      int tempi = QualityMatrix[nextBest.pi][nextBest.pj].pi;
      int tempj = QualityMatrix[nextBest.pi][nextBest.pj].pj;
      nextBest.pi = tempi;
      nextBest.pj = tempj;
    }
    alignments.push_back(alignment);
    clearAlignmentBox(alignment, QualityMatrix);
    //clearAlignment(alignment,QualityMatrix);

    Match temp;
    //TODO: is it worth making this into a ctor
    temp.quality = best.q;
    temp.s1.start = alignment[alignment.size() - 1][0];
    temp.s1.end = alignment[0][0]; //best.pi;
    temp.s1.fname = novel.fname;
    temp.s2.start = alignment[alignment.size() - 1][1];
    temp.s2.end = alignment[0][1];//best.pj;
    temp.s2.fname = base.fname;
    temp.classification = -1;

    if (temp.s1.end - temp.s1.start > minlen && temp.s2.end - temp.s2.start > minlen) {

      //scale to global scope
      temp.s1.start = 160 * temp.s1.start + novel.start;
      temp.s2.start = 160 * temp.s2.start + base.start;
      temp.s1.end = 160 * temp.s1.end + novel.start;
      temp.s2.end = 160 * temp.s2.end + base.start;
      //temp.s1.start = temp.s1.start;
      //temp.s1.end = temp.s1.end;
      //temp.s2.start = temp.s2.start;
      //temp.s2.end = temp.s2.end;
      cout<<temp.s1.fname.index<<" "<<temp.s2.fname.index<<endl;
      cout<<temp.quality<<" "<<temp.s1.start/16000.0<<" "<<temp.s1.end/16000.0<<" "<<temp.s2.start/16000.0<<" "<<temp.s2.end/16000.0<<endl;
      matches.push_back(temp);
    }
    else numLines--;


    //get next best Q/
    best.q = 0;
    //TODO: you can start at 1 here since it wont ever be in the first row or col
    for (int i = 0; i < QualityMatrix.size(); i++) {
      for (int j = 0; j < QualityMatrix[0].size(); j++) {
        if (QualityMatrix[i][j].q > best.q) {
          best.q = QualityMatrix[i][j].q;
          best.pi = i;
          best.pj = j;
        }
      }
    }
    numLines++;
  }

  if (visualize) { //generate matlab plots of subsequence alignemnts
    cout<<"plotting\n";
    string fname = "../plots/plot";
    fname += to_string(novel.fname.index);//novel.fname.name + novel.fname.file;
    fname += "and";
    fname += to_string(base.fname.index);//base.fname.name + base.fname.file;
    fname += ".m";
    visualizeLine(Original, alignments, fname);
  }
  return matches;
}

//Visulization ****************************************************************
//generate matlab plot of quality matrix overlayed with subsequence alignments
void Comparison::visualizeLine(vector<vector<qmCell>> quality,
                               vector<vector<vector<int>>> alignments,
                               string fname) {
  string out;

  string mat = "m = [ ";

  for (int i = 0; i < quality.size(); i++) {
    for (int j = 0; j < quality[i].size(); j++) {
      mat += to_string(quality[i][j].q);
      mat += " ";
    }
    mat += "\n";
  }

  mat += " ];\n";

  mat += "imagesc(m);\n";
  mat += "set(gca,'DataAspectRatio',[1 1 1]);\n";

  out += mat;

  for (int i = 0; i < alignments.size(); i++) {
    string x = "x" + to_string(i) + " = [ ";
    string y = "y" + to_string(i) + " = [ ";
    for (vector<int> xy: alignments[i]) {
      x += to_string(xy[1] + 1);
      x += " ";
      y += to_string(xy[0] + 1);
      y += " ";
    }
    x += "];\n";
    y += "];\n";

    string line = "line(x" + to_string(i) + ",y" + to_string(i) + ",'Color',[1,1,1]);\n";
    out += x;
    out += y;
    out += line;
  }
  out += "title(\'" + fname + "\');\n";
  ofstream fout(fname);
  fout<<out;
  fout.close();
}

//create matlab plot of distance matrix
void Comparison::visualizeD(vector<vector<double> > distance) {
  string out = "d = [ ";

  for (int i = 0; i < distance.size(); i++) {
    for (int j = 0; j < distance[i].size(); j++) {
      out += to_string(distance[i][j]);
      out += " ";
    }
    out += "\n";
  }

  out += " ];\n";

  out += "imagesc(d);\n";
  out += "set(gca,'DataAspectRatio',[1 1 1]);";
  ofstream fout("../plots/dmat.m");
  fout<<out;
  fout.close();
}