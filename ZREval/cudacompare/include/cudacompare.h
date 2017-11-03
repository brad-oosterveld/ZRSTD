//REV 17 March 2016
//cudacompare.h: Functions called by host code (ZREval.h) to compute
//DTW of a list of chunks

//Computation for each chunk is:
//Euclid (or cosine dist) between each time point. Each time point is a C-dim vector, C is length of MFCC
//Each chunk is:
//A vector of these time points of length T, both .raw and .smooth
//It also has index in original file (start and end).
//Also, it has finfo fname (which is name (of chunk/word/etc.?), and file, which is filename?)

#pragma once

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <string>

#include "../../DataTypes.h"

//typedef double real_t;
typedef float real_t;


//REV: Need to make "workspace" for each, or call pure cuda_malloc?
//I can only know the number of workspaces if I know # of workers...
//Which I can compute beforehand actually haha, by querying the GPU.
//I can select nworkers based on how much memory I have basically

#define matrix_x( arrayval, xsize )  ((arrayval) % (xsize))
#define matrix_y( arrayval, xsize )  ((arrayval) / (xsize))
#define arrayloc( x, y, xsize ) (((y)*(xsize)) + (x))

#define intdivceil( x, y ) (((x)+(y)-1)/(y))

/*-gencode arch=compute_37,code=sm_37struct finfo
{
  std::string name;
  std::string file;

finfo()
:
  name(""), file("")
  {
    //REV: Nothing
  }
};


struct segment
{
  finfo fname;
  int start;
  int end;
  std::vector<double> PCM;
  std::vector<std::vector<double>> smooth;
  std::vector<std::vector<double>> raw;
};

struct Match
{
  double quality;
  int start1;
  int end1;
  finfo chunk1;
  int start2;
  int end2;
  finfo chunk2;
  int classification;
};
*/
//vector<vector<Match>> ZREval::findSubSequencesCUDA( const std::vector<segment>& segments,
//std::vector<std::vector<Match>> findSubSequencesCUDA( const std::vector< segment >& segments,
/*std::vector<std::vector<Match>> findSubSequencesCUDA( const std::vector< segment >& segments,
						      const bool& useraw,
						      const bool& usecosine,
						      const int& maxnmatches
						      );*/

std::vector<std::vector<std::vector<Match>>> findSubSequencesCUDA( const std::vector<segment>& segments,
								   const bool& usesmooth,
								   const bool& usecosine,
								   const int& maxnmatches,
								   const real_t& _btmin,
								   const real_t& _qThresh,
								   const real_t& _btminFact,
								   const int& _minmatchlen,
								   const real_t& _penalty,
								   const real_t& _bonus,
								   const int& cudadev=0
								   );

