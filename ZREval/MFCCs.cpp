//
// Created by brad on 1/4/16.
//

#include "MFCCs.h"
#include "armadillo"
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
//#include <log4cxx/logger.h>
#include <iostream>
//#include <boost/foreach.hpp>

using namespace std;

//initialize everything to default params
MFCCs::MFCCs(){// :
//  logger(log4cxx::Logger::getLogger("ade.asr.MFCCs")) {
  normVolume = false;
  numDiffs = 2;
  smoothFactor = .25;
}

MFCCs::~MFCCs() {

}

vector<vector<double>> MFCCs::getFeatures(std::vector<double> PCMData) {
  vector<vector<double>> out;

  //TODO: is it worth getting rid of this redundancy?
  //calculating minimum possible size for PCMData
  int fs = 16000;
  int Tw = 25;
  int Ts = 10;

  int numWindows = (int) round((1E-3) * Tw * fs);  //frame duration in samples
  int numShift = (int) round((1E-3) * Ts * fs);

  int minSegmentSize = -numShift + numWindows;

  if ((((int) PCMData.size() - numWindows) / numShift + 1) > 0) {

    //normalize volume
    if (normVolume) PCMData = normalizeVolume(PCMData);

    //caluclate MFCCS
    //will always smooth, can set smooth PCT to 0 for no smoothing
    vector<vector<double>> mfccs = smooth(calculateMFCCs(PCMData));

    if (numDiffs > 0) {

      vector<vector<double>> firstDiffs;
      vector<vector<double>> secondDiffs;

      if (numDiffs >= 1) {
        firstDiffs = smooth(getDiffs(mfccs));
      }
      if (numDiffs >= 2) {
        secondDiffs = smooth(getDiffs(firstDiffs));
      }
      if (numDiffs >= 3) {
        //LOG4CXX_WARN(logger, "we only support 2 diffs, and that's probably all you want");
        cout<< "we only support 2 diffs, and that's probably all you want";
      }

      //combine into single feature vector
      for (int i = 0; i < mfccs.size(); i++) {
        vector<double> combinedFeature;
        combinedFeature.insert(combinedFeature.end(), mfccs[i].begin(), mfccs[i].end());
        if (numDiffs >= 1) {
          combinedFeature.insert(combinedFeature.end(), firstDiffs[i].begin(),
                                 firstDiffs[i].end());
        }
        if (numDiffs >= 2) {
          combinedFeature.insert(combinedFeature.end(), secondDiffs[i].begin(),
                                 secondDiffs[i].end());
        }
        out.push_back(combinedFeature);
      }
    }
    else {
      out = mfccs;
    }
  }
  else {
    //LOG4CXX_WARN(logger,
    //             "No MFCCs extracted because input segment is too small.\n segment size:"<<
    //             PCMData.size()<<" min segment size: "<<minSegmentSize<<
    //             ". Try changing filter params");
  }
  //LOG4CXX_DEBUG(logger, "MFCCs generated: "<<out.size());
 // cout<<"MFCCs generated: "<<out.size();
  return out;
}

vector<double> MFCCs::normalizeVolume(vector<double> raw) {
  vector<double> norm(raw.size(), 0);
  double maxPCM = 0;
  for (double r: raw) {
    if (r > maxPCM)maxPCM = r;
  }
  for (int i = 0; i < raw.size(); i++) {
    norm[i] = raw[i] / maxPCM;
  }
  return norm;
}

vector<vector<double>> MFCCs::smooth(vector<vector<double>> rough) {
  vector<double> padding(rough[0].size(), 0.0);
  vector<vector<double>> smoothed;
  smoothed.push_back(padding);

  //double smoothFactor = .1;

  for (int i = 0; i < rough.size(); i++) {
    vector<double> smoothedPoint(rough[i].size(), 0);
    for (int k = 0; k < smoothedPoint.size(); k++) {
      smoothedPoint[k] +=
        smoothed[i][k] * (1.0 - smoothFactor) + rough[i][k] * smoothFactor;
    }
    smoothed.push_back(smoothedPoint);
  }
  smoothed.erase(smoothed.begin());
  return smoothed;
}

vector<vector<double>> MFCCs::getDiffs(vector<vector<double>> original) {
  vector<vector<double>> diffs(original.size(),
                               vector<double>(original[0].size(), 0));
  //insert copies of first and last time step
  vector<double> first = original[0];
  original.insert(original.begin(), first);
  vector<double> last = original[original.size() - 1];
  original.insert(original.end(), last);

  for (int i = 1; i < original.size() - 1; i++) {
    for (int j = 0; j < original[i].size(); j++) {
      diffs[i - 1][j] = (original[i + 1][j] - original[i - 1][j]) / 2.0;
    }
  }
  return diffs;
}

//cepstral lifter routine
inline vector<double> MFCCs::ceplifter(double n, double l) {
  vector<double> lift;
  for (int i = 0; i < n; i++) {
    lift.push_back(1 + 0.5 * l * sin(M_PI * i / l));
  }
  return lift;
}

inline int MFCCs::NextPowerOf2(int val) {
  val--;
  val = (val>>1) | val;
  val = (val>>2) | val;
  val = (val>>4) | val;
  val = (val>>8) | val;
  val = (val>>16) | val;
  return ++val;
}

//htk pre-emphasis function
vector<double> MFCCs::preEmphasis(double coeff, vector<double> x) {
  vector<double> y(x.size(), 0);
  y[0] = x[0];
  for (int n = 1; n < x.size(); n++) {
    y[n] = x[n] - coeff * x[n - 1];
  }
  return y;
}

//hamming window a la matlab
vector<double> MFCCs::window(int windowSize) {
  vector<double> windowed(windowSize, 0);
  for (int n = 0; n < windowSize; n++) {
    windowed[n] = 0.54 - .46 * cos((2 * M_PI) * ((double) n / windowSize - 1));
  }
  return windowed;
}

//based on mike's adaptation of  Kamil Wojcicki's "HTK MFCC MATLAB CODE"
//http://www.mathworks.com/matlabcentral/fileexchange/32849-htk-mfcc-matlab/content/mfcc/mfcc.m
vector<vector<double>> MFCCs::calculateMFCCs(vector<double> signal) {
  vector<vector<double>> out;
  //ofstream fout;
  //fout.open("mfccstuff.txt");

  //params hard coded for now
  int fs = 16000;
  int Tw = 25;
  int Ts = 10;
  double alpha = .97;
  int numFilterBankChannels = 26;
  int numCC = 13;//12 + 1;// but don't use the first
  int cepstralSineLifter = 22;

  double fMin = 0; // filter coefficients start at this frequency (Hz)
  double fLow = 20; //minFreq;//       % 20 lower cutoff frequency (Hz) for the filterbank
  double fHigh = 4400; //maxFreq;//      % 4400 upper cutoff frequency (Hz) for the filterbank
  double fMax =
    0.5 * fs; //     % filter coefficients end at this frequency (Hz)

  int numWindows = (int) round(
    (1E-3) * Tw * fs);  // % frame duration (in samples)
  int numShift = (int) round((1E-3) * Ts * fs);   //% frame shift (in samples)

  double power = log2(NextPowerOf2(numWindows));
  int nfft = pow(2, power);     // %length of FFT analysis
  int uniqueFFT = (nfft / 2) + 1;    //% length of the unique part of the FFT

  //preemphasize signal (high pass filter)

  //setting up frame buckets for fft
  //int numFrames = int((filteredSpeech.size() - numWindows) / (numShift) + 1);
  //vector<vector<double>> frames;
  //for (int i = 0; i < numFrames; i++) {
  //  vector<double> nextFrame(filteredSpeech.begin() + (i * numShift),
  //                           filteredSpeech.begin() + (i * numShift) + numWindows);
  //  frames.push_back(nextFrame);
  //}
  int numFrames = ((int) (signal.size()) - numWindows) / (numShift) + 1;
  vector<vector<double>> frames;
  for (int i = 0; i < numFrames; i++) {
    vector<double> nextFrame(signal.begin() + (i * numShift),
                             signal.begin() + (i * numShift) + numWindows);
    //frames.push_back(filter(filterVec,1,nextFrame));
    frames.push_back(preEmphasis(alpha, nextFrame));
  }

  //create matrix for application of Hamming Window
  vector<double> windowed = window(numWindows);
  arma::vec aWindowed(windowed);
  arma::mat windowMat = arma::diagmat(aWindowed);

  //convert to armadillo for fft and matrix multiplication
  arma::mat mFrames(frames[0].size(), frames.size());
  for (int i = 0; i < frames[0].size(); i++) {
    for (int j = 0; j < frames.size(); j++) {
      mFrames(i, j) = frames[j][i];
    }
  }

  //Apply hamming window to frame buckets
  arma::mat wFrames = windowMat * mFrames;

  //get energy
  double silenceFloor = 50; //50db defaults from HTK
  double energyScale = .1; //default vlaue from HTK
  vector<double> energy(wFrames.n_cols, 0);

  for (int i = 0; i < wFrames.n_cols; i++) {
    double e = 0;
    for (int j = 0; j < wFrames.n_rows; j++) {
      e += pow(wFrames(j, i), 2);
    }
    energy[i] = log(e);
  }

  //normalize energy
  double maxEnergy = 0;
  for (double e: energy) {
    if (e > maxEnergy)maxEnergy = e;
  }
  double minEnergy = maxEnergy - (silenceFloor * log(10.0)) / 10.0;
  for (int i = 0; i < energy.size(); i++) {
    if (energy[i] < minEnergy) energy[i] = minEnergy;
    energy[i] = 1.0 - (maxEnergy - energy[i]) * energyScale;
  }

  // Magnitude spectrum computation
  // could implement HTK's c++ fft here to get results closer to theirs, not sure if its worth it
  arma::mat MAG = arma::abs(arma::fft(wFrames, nfft));

  //calculate frequency range in Hz, based on sized of unique part fo fft
  vector<double> fRange(uniqueFFT, 0);
  for (int i = 0; i < uniqueFFT; i++) {
    fRange[i] = i * (fMax / (uniqueFFT - 1));
  }

  double hz2melf_low = 2595 * log10(1 + fLow / 700);
  double hz2melf_high = 2595 * log10(1 + fHigh / 700);

  vector<double> melChannels;
  for (int i = 0; i <= numFilterBankChannels + 1; i++) {
    double temp = hz2melf_low + i * ((hz2melf_high - hz2melf_low) /
                                     (numFilterBankChannels + 1));
    double x = 700 * (pow(10, (temp / 2595)) - 1);
    melChannels.push_back(x);
  }

  // create filter bank matrix to put magnitude spectrum into Mel Frequency Bands
  // super contrived, based on doing it in MATLAB but probably doesn't need to be a matrix;
  arma::mat H = arma::zeros(numFilterBankChannels, uniqueFFT);
  for (int m = 0; m < numFilterBankChannels; m++) {
    for (int k = 0; k < uniqueFFT; k++) {
      if (fRange[k] >= melChannels[m] && fRange[k] < melChannels[m + 1]) {
        H(m, k) =
          (fRange[k] - melChannels[m]) / (melChannels[m + 1] - melChannels[m]);
      }
      if (fRange[k] >= melChannels[m + 1] && fRange[k] <= melChannels[m + 2]) {
        H(m, k) = (melChannels[m + 2] - fRange[k]) /
                  (melChannels[m + 2] - melChannels[m + 1]);
      }
    }
  }

  // GET FREQUENCY BAND ENERGY
  // Filterbank application to unique part of the magnitude spectrum
  arma::mat FBE =
    H * MAG(arma::span(0, uniqueFFT - 1), arma::span(0, frames.size() - 1));

  // DCT Matrix Computation
  arma::mat DCT = arma::zeros(numCC, numFilterBankChannels);
  for (int i = 0; i < numCC; i++) {
    for (int j = 0; j < numFilterBankChannels; j++) {
      DCT(i, j) = sqrt(2.0 / numFilterBankChannels) *
                  cos((i) * (M_PI * (j + 1 - .5) / numFilterBankChannels));
    }
  }

  // Conversion of logFBEs to cepstral coefficients through DCT
  arma::mat CC = DCT * arma::log(FBE);

  // Cepstral lifter computation
  arma::vec lifter(ceplifter(numCC, cepstralSineLifter));
  arma::mat lifterMat = arma::diagmat(lifter);

  // Cepstral liftering gives liftered cepstral coefficients
  arma::mat MFCCs = lifterMat * CC;

  //convert back to std::vector
  typedef vector<double> stdvec;
  for (int i = 0; i < MFCCs.n_cols; i++) {
    //vector<double> features(MFCCs.begin_col(i) + 1, MFCCs.end_col(i));
    //features.push_back(energy[i]);

    //this order for compatability with Mike's code
    vector<double> features;
    features.push_back(energy[i]);
    features.insert(features.end(), MFCCs.begin_col(i) + 1, MFCCs.end_col(i));

    out.push_back(features);
  }

  //for (int i = 0; i < out.size(); i++) {
  //  for (double j: out[i]) {
  //    fout << setw(10) << j << " ";
  //  }
  //  fout << endl;
  //}
  //
  //fout.close();
  return out;
}