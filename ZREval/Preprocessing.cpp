//
// Created by brad on 3/9/16.
//

#include "Preprocessing.h"
#include "dirent.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include "Timer.h"

using namespace std;

/* Overhead *************************************************************/

Preprocessor::Preprocessor(int lang) :
//logger(log4cxx::Logger::getLogger("astd.preprocessing")),
  raw(),
  smooth() {
  language = lang;
  // raw.loadConfig("rawConfigpath");
  // smooth.loadConfig("smoothConfigpath");
}

Preprocessor::~Preprocessor() {

}

/* Main Method ***************************************************************/

std::vector<segment> Preprocessor::preProcessDirectory(std::string path) {

  vector<pair<finfo, vector<double>>> files = loadFiles(path);

  vector<segment> segmented;

  Timer segtimer;
  //TODO: parallelize this
  for (auto f: files) {
    vector<segment> currentFile = segmentFile(f.first, f.second);
    segmented.insert(segmented.end(), currentFile.begin(), currentFile.end());
  }
  cout<<"FINISHED files segmented. Took "<<segtimer.elapsed()<<" seconds"<<std::endl;

  Timer feattimer;

#pragma omp parallel for
  for (size_t s = 0; s < segmented.size(); ++s)
    //  for (segment &s: segmented)
  {
    extractFeatures(segmented[s]);
    segmented[s].PCM.clear();
    segmented[s].PCM.shrink_to_fit();
  }

  cout<<"FINISHED features extracted. Took "<<feattimer.elapsed()<<" seconds"<<std::endl;

  return segmented;
}


/* Loading files *************************************************************/

//loads a directory's worth of .wav files
vector<pair<finfo, vector<double>>> Preprocessor::loadFiles(string path) {


  DIR* dir;
  struct dirent* ent;
  std::vector<std::string> fnames;
  if ((dir = opendir(path.c_str())) != NULL) {
    while ((ent = readdir(dir)) != NULL) {
      //checking file names for the right format
      if (ent->d_name[0] != '.' && ent->d_name[1] != '.') {
        //needs to be the same dir name as path
        fnames.push_back(ent->d_name);
      }
    }
    closedir(dir);
  }
  else {
    // LOG4CXX_ERROR(logger, "could not open directory "<<path);
    cout<<"could not open directory "<<path;
  }

  sort(fnames.begin(), fnames.end());

  vector<pair<finfo, vector<double>>> output(fnames.size());

//  for (std::string fname:fnames) {
#pragma omp parallel for
  for (int i = 0; i < output.size(); i++) {
    std::string full = path + std::string("/") + fnames[i];
    std::ifstream fin;
    fin.open(full);
    std::vector<char> raw;
    char bits = fin.get();

    while (fin.good()) {
      raw.push_back(bits);
      bits = fin.get();
    }
    fin.close();


    finfo fileinfo;
    //diffrent corpora have diffrent file naming conventions
    if (language == 0) {
      fileinfo.name = fnames[i].substr(0, 5);
      fileinfo.file = fnames[i][5];
    }
    else if (language == 1) {
      fileinfo.name = fnames[i].substr(0, 15);
      fileinfo.file = fnames[i].substr(15, 4);
    }
    output[i] = (make_pair(fileinfo, quantizePCM(raw)));
  }

  //LOG4CXX_TRACE(logger, fnames.size()<<" files read");
  cout<<fnames.size()<<" files read \n";
  return output;
}

/* Segmentation *************************************************************/

vector<segment> Preprocessor::segmentFile(finfo name, vector<double> data) {

  vector<segment> foundSegments;
  double runningAvg = 0.0;
  int remainingFrames = 0;

  double runningAvgPct = .0066;
  int windowSize = 960;
  double envelopeThreshold = .01;
  int frontPad = 1200;
  int backPad = 480; //really the opposite of padding
  int decayFrames = 3200;

  segment current;
  current.fname = name;

  for (int i = 0; i < data.size(); i++) {
    //full wave rectify current frame
    double thisFrame = fabs(data[i]);
    //calculate running average
    runningAvg = runningAvg * (1 - runningAvgPct) + thisFrame * runningAvgPct;

    //if above threshold reset number of decay frames
    if (runningAvg > envelopeThreshold) {
      //TODO:make sure this doesn't overlap with other existing segments
      //pad the beginning of the segment with some of the preceeding frames
      if (current.PCM.size() == 0) {
        int padding = min(frontPad, i);
        int startPos = max(0, i - frontPad);
        current.PCM.assign(data.begin() + startPos, data.begin() + startPos + padding);
      }
      remainingFrames = decayFrames;
    }

    //if there are still remaining frames add to current segment
    if (remainingFrames > 0) {
      current.PCM.push_back(data[i]);
    }
    else if (current.PCM.size() > 0) {
      //cout<<"segment found! "<<current.PCM.size()<<" ";
      current.end = i - backPad;
      current.start = max(0, int(i - current.PCM.size()));
      current.fname.index = (int)foundSegments.size();
      //copy current into found segments
      foundSegments.push_back(current);
      //resets current to an empty segment;
      current.PCM.resize(0);
    }
    remainingFrames--;
  }
  if (current.PCM.size() > 0) {
    current.end= (int)data.size()-1;
    current.start=max(0, int(current.end - current.PCM.size()));
    current.fname.index = (int)foundSegments.size();
    foundSegments.push_back(current);
  }
  return foundSegments;
}

/* Feature Extraction ********************************************************/

//TODO: is doing this as a reference worth it?
void Preprocessor::extractFeatures(segment &data) {

  //TODO: check if it would matter that we're doing a bunch of the same computation twice
  //data.raw = raw.getFeatures(data.PCM);
  data.smooth = smooth.getFeatures(data.PCM);
}

/* Util *************************************************************/
std::vector<double> Preprocessor::quantizePCM(std::vector<char> raw) {
  int i;
  int numFrames;
  std::vector<double> quantized;

  i = 44; //remove .wav header
  numFrames =
    ((0xff & raw.data()[43])<<24) + ((0xff & raw.data()[42])<<16) +
    ((0xff & raw.data()[41])<<8) + (0xff & raw.data()[40]);
  //cout<<numFrames<<" "<<raw.size()-44<<endl;

  //cout<<raw.size()<<endl;
  for (i; i < numFrames; i += 2) {
    short firstByte = 0xff & raw[i];
    short secondByte = 0xff & raw[i + 1];
    short out = (secondByte<<8) + firstByte;
    double q = out / 32768.0f;
    //cout<<q<<endl;
    quantized.push_back(q);
  }
  //cout<<numFrames<<" "<<quantized.size()<<endl;
  return quantized;
}
