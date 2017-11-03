//
// Created by brad on 3/9/16.
//

#include "ZREval.h"
#include <iostream>
#include <fstream>
#include <algorithm> //for sorting matches
#include "cudacompare/include/cudacompare.h"
#include "Timer.h"

using namespace std;

ZREval::ZREval(int lang) :
  preprocessor(lang) {
  overlapRange = 1600; //in frames
}

ZREval::~ZREval() {

}


void ZREval::runEvalFromCfg(std::string path) {
//TODO: check if segments cover the complete length of the input file or if wee need to remember the file size
  Timer preproctimer;
  prepocessed = preprocessor.preProcessDirectory(path);
  cout<<"Preprocessing finished. Took " << preproctimer.elapsed() << " seconds. Segments found:" << prepocessed.size() << "\n";
  
  printSegments(prepocessed);
  cout<<"segments printed \n";

  vector<vector<vector<Match>>> subsequences = findSubSequences(prepocessed);
  cout<<"Comparison complete.\n";

  printMatches(subsequences);
  cout<<"Subsequences printed\n";

  //cluster subsequences
  generateClusters(subsequences);
  cout<<"Clusters generated \n";

  //refine clusters

  //complete transcription
  generateTranscription();
  cout<<"Transcription written \n";
}

vector<vector<vector<Match>>> ZREval::findSubSequences(std::vector<segment> segments) {

  fprintf(stdout, "REV: STARTING matches (GPU)\n");
  vector<vector<vector<Match>>> gpudiscovered;

  bool usecosine = true;
  bool usesmooth = true;
  int maxmatchesper = 10;

  double btmin = 0.1;
  double qThresh = 19;
  double btminFact = 0.05;
  int minmatchlen = 30;
  double penalty = -1;
  double bonus = 1.0;
  Timer gputimer;
  gpudiscovered = findSubSequencesCUDA(segments, usesmooth, usecosine, maxmatchesper, btmin, qThresh,
                                       btminFact, minmatchlen, penalty, bonus);
  fprintf(stdout, "FINISHED GPU. GPU total time: [%lf] seconds\n", gputimer.elapsed());
  

  //REV: re-constructing DPNgram here...so that we can make sure we use same params...
  DPNgramComparison = Comparison( qThresh, btminFact, maxmatchesper, minmatchlen, bonus, penalty, btmin );
  
  //#define COMPARECPU
#ifdef COMPARECPU
  vector<vector<vector<Match>>> cpudiscovered;


  fprintf(stdout, "REV: Starting discovery of matches (CPU)\n");

  int cpucompcount = 0;

  //REV: Are we doing upper or lower triangle?
  //If first index is row, second is column, then let's say, zero-th row should be full, but for itself.
  //And zero, zero should be left blank (I guess). So, (0, 1) is index 1. Fine, let's include the diagnal then?
  //We need to decide how we do it...
  //Anyway, I made it so that [0][0] in your matrix contains "same" guys.
  
  
  Timer cputimer;
  for (int i = 0; i < segments.size(); i++) {

    int segSize = segments.size() - i;
    vector<vector<Match>> current(segSize);
#pragma omp parallel for
    //REV: is there a reason you re-added diagonal computation? Indices are mixed up now... if j=i, should be j-i
    //and you only had segsize - 1...so it would be out of range...
    for (int j = i; j < segments.size(); j++) {
      if (i != j) {
        current[j - i] = DPNgramComparison.compare(segments[i], segments[j], usesmooth);
        //cpucompcount++; //REV: This may cause issues b/c might read old value then iter and lose values, i.e.
	//race condition
      }
      //else {
      //  //TODO:implement comparison with self
      //}
    }

    cpudiscovered.push_back(current);
    //cout<<i<<" ";
  }

  //DPNgramComparison.compare(segments[22], segments[37], true);

  fprintf(stdout, "CPU Finished. CPU time elapsed [%lf] seconds\n", cputimer.elapsed());

  //cout<<"cpucompcount: "<<cpucompcount<<endl;
  
  int gpucount = 0;
  int cpucount = 0;

  if (gpudiscovered.size() != cpudiscovered.size()) {
    cout<<"first dimension not the same: "<<gpudiscovered.size()<<" "<<cpudiscovered.size()<<endl;
  }
  for (size_t i = 0; i < gpudiscovered.size(); i++) {
    if (gpudiscovered[i].size() != cpudiscovered[i].size()) {
      cout<<"second dimension not the same: "<<i<<" "<<gpudiscovered[i].size()<<" "<<
      cpudiscovered[i].size()<<endl;
    }
    //fprintf(stdout, "[%ld]-th row of gpudiscovered has [%ld] columns (should have row-1)\n", i, gpudiscovered[i].size() );
    for (size_t j = 0; j < gpudiscovered[i].size(); j++) {
      if (gpudiscovered[i][j].size() != cpudiscovered[i][j].size()) {
        cout<<"number of found matches not the same (seg1  seg2  gpufound  cpufound): "<<i<<" "<<j<<" "<<gpudiscovered[i][j].size()<<
        " "<<cpudiscovered[i][j].size()<<endl;
      }
      if (gpudiscovered[i][j].size() > 0 && cpudiscovered[i][j].size() > 0) {
        cout<<i<<" "<<j<<" "<<gpudiscovered[i][j].size()<<" "<<cpudiscovered[i][j].size()<<endl;
      }
      gpucount += gpudiscovered[i][j].size();
      cpucount += cpudiscovered[i][j].size();

      for (size_t m = 0; m < cpudiscovered[i][j].size(); m++) {

        //cout<<gpudiscovered[i][j][m].s1.fname.name<<gpudiscovered[i][j][m].s1.fname.file<<" "<<cpudiscovered[i][j][m].s1.fname.name<<cpudiscovered[i][j][m].s1.fname.file<<endl;
        //cout<<gpudiscovered[i][j][m].s2.fname.name<<gpudiscovered[i][j][m].s2.fname.file<<" "<<cpudiscovered[i][j][m].s2.fname.name<<cpudiscovered[i][j][m].s2.fname.file<<endl;


        if (gpudiscovered[i][j].size() > m) {
          cout<<"gpu: "<<gpudiscovered[i][j][m].s1.start / 16000.0<<" "<<
          gpudiscovered[i][j][m].s1.end / 16000.0<<" "<<
          gpudiscovered[i][j][m].s2.start / 16000.0<<" "<<
          gpudiscovered[i][j][m].s2.end / 16000.0<<endl;
        }

        if (cpudiscovered[i][j].size() > m) {
          cout<<"cpu: "<<cpudiscovered[i][j][m].s1.start / 16000.0<<" "<<
          cpudiscovered[i][j][m].s1.end / 16000.0<<" "<<
          cpudiscovered[i][j][m].s2.start / 16000.0<<" "<<
          cpudiscovered[i][j][m].s2.end / 16000.0<<endl;
        }
      }
    }
  }

  cout<<"gpucount: "<<gpucount<<" cpucount: "<<cpucount<<
  endl;

#endif

  return gpudiscovered;
//return cpudiscovered;
}


/* Clustering ****************************************************************/
void ZREval::generateClusters(vector<vector<vector<Match>>> &pairs) {

  for (size_t i = 0; i < pairs.size(); i++) {

    //generate clusters for all unclassified matches from segment i

    vector<Cluster> newClusters = groupByFinfo1(pairs[i]);
    //if (pairs[i].size() > 0) {
    //  cout<<"pairs[i] size "<<pairs[i].size()<<endl;
    //  cout<<"newClusters size "<<newClusters.size()<<endl;
    //}

    //search the rest of pairs for connected components
    // std::cout<<"nc size: "<<newClusters.size()<<endl;
    for (Cluster &c: newClusters) {
      c.addConnected(pairs); //need to be careful about the indexing here
    }
    //add newly created clusters to list of discovered classes;
    if (newClusters.size() > 0) {
      classes.insert(classes.end(), newClusters.begin(), newClusters.end());
      //cout<<"new clusters: "<<newClusters.size()<<" classes: "<<classes.size()<<endl;

    }
  }
}


bool matchSort(Match m1, Match m2) { return (m1.s1.start < m2.s1.start); }

//creates a clusters for non overlapping s1's in matches
//overlapping s1's are added to the same cluster
vector<Cluster> ZREval::groupByFinfo1(vector<vector<Match>> &matches) {


  //flatten
  vector<Match> flat;

  for (int i = 0; i < matches.size(); i++) {
    for (int j = 0; j < matches[i].size(); j++) {
      flat.push_back(matches[i][j]);
    }
  }

  sort(flat.begin(), flat.end(), matchSort);
  //sort matches by s1.start

  vector<Cluster> clusters(0);
  for (Match &m: flat) {
    //don't do anything to clustered matches
    if (m.classification == -1) {
      bool found = false;
      //see if there is already a cluster that contains s1
      for (Cluster &c: clusters) {
        if (c.s1InRange(m.s1) != -1) {
          c.addMatch(m);
          found = true;
        }
      }
      if (!found) { //make a new cluster with m as the seeed
        Cluster newCluster(m, clusters.size() + classes.size(), overlapRange);
        clusters.push_back(newCluster);
      }
    }
  }

  //cout<<"grouped by f1 info: "<<clusters.size()<<endl;
  return clusters;
}

/* Transcription *************************************************************/

void ZREval::generateTranscription() {
  //TODO: figure out if we should add silence or not

  generateSingletons();

  cout<<"Singletons Generated \n";

  std::ofstream fout;
  fout.open("DiscoveredClasses.txt");

  //go through multi member classes first
  cout<<"number of discovered classes: "<<classes.size()<<endl;
  for (int i = 0; i < classes.size(); i++) {
    fout<<"Class "<<i<<endl;
    //fout<<flat[i].quality<<endl;
    for (auto s: classes[i].refined) {
      fout<<s.second.fname.name<<s.second.fname.file<<" "<<
      (double) s.second.start / 16000.0<<" "<<(double) s.second.end / 16000.0<<endl;
      //fout<<"s"<<flat[i].s2.fname.name<<flat[i].s2.fname.file<<" "<<
      //(double) flat[i].s2.start / 16000.0<<" "<<(double) flat[i].s2.end / 16000.0<<endl;
    }
    fout<<endl;
  }

  //then go through singletons
  for (int i = 0; i < singletons.size(); i++) {
    fout<<"Class "<<i + classes.size()<<endl;
    fout<<singletons[i].fname.name<<singletons[i].fname.file<<" "<<
    (double) singletons[i].start / 16000.0<<" "<<(double) singletons[i].end / 16000.0<<endl;
    fout<<endl;
  }
  fout.close();

}

void ZREval::generateSingletons() {

  vector<TranscriptionUnit> TUs;

  //create transcription units from segments
  for (segment s: prepocessed) {
    TranscriptionUnit next;
    next.segmentID = s.fname.index;
    Class base(s, false);
    next.subSequences.push_back(base);
    TUs.push_back(next);
  }

  //add in discovered classes to generate singleton Classes from un matched regions
  for (Cluster c: classes) {
    for (auto r: c.refined) {
      Class newSS(r.second, true);
      TUs[r.first].insertClass(newSS);
    }
  }

  //go through TUs and add every Class that isn't from a cluster

  for (TranscriptionUnit TU: TUs) {
    for (Class c: TU.subSequences) {
      if (!c.fromCluster) singletons.push_back(c);
    }
  }
}


/* Printing ******************************************************************/

void ZREval::printSegments(vector<segment> segments) {

  std::ofstream fout;

  fout.open("segList.txt");

  for (int i = 0; i < segments.size(); i++) {
    fout<<"Class "<<i<<endl;
    fout<<segments[i].fname.name<<segments[i].fname.file<<" "<<
    (double) segments[i].start / 16000.0<<" "<<
    (double) segments[i].end / 16000.0<<std::endl;
    fout<<endl;
  }
  fout.close();

}

void ZREval::printMatches(vector<vector<vector<Match>>> matches) {

  vector<Match> flat;

  for (int i = 0; i < matches.size(); i++) {
    for (int j = 0; j < matches[i].size(); j++) {
      for (int k = 0; k < matches[i][j].size(); k++) {
        flat.push_back(matches[i][j][k]);
      }
    }
  }
  std::ofstream fout;
  cout<<"flattened matches found: "<<flat.size()<<" ";

  fout.open("matchList.txt");

  for (int i = 0; i < flat.size(); i++) {
    fout<<"Class "<<i<<endl;
    //fout<<flat[i].quality<<endl;
    fout<<flat[i].s1.fname.name<<flat[i].s1.fname.file<<" "<<
    (double) flat[i].s1.start / 16000.0<<" "<<(double) flat[i].s1.end / 16000.0<<endl;
    fout<<flat[i].s2.fname.name<<flat[i].s2.fname.file<<" "<<
    (double) flat[i].s2.start / 16000.0<<" "<<(double) flat[i].s2.end / 16000.0<<endl;
    fout<<endl;
  }

  fout.close();
}
