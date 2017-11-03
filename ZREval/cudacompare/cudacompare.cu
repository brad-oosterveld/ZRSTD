#include "include/helper_cuda.h"
#include "include/cudacompare.h"
#include <cuda.h>
#include <cuda_runtime.h>




//REV: I need to ensure that threads_per_block works...fuck
//Need to ensure that i'm doing it per-block basically, for "workers". Heh.


__device__
int
matrixmax(int workeroffset, int threadsperworker, int npts,
          real_t* targmatdata, int* mymaxesptr) {
  int mymaxidx = 0;
  int nperchunk = intdivceil(npts, threadsperworker); //ceil(npts / threadsperworker);

  //REV: Max in each worker thread
  if (workeroffset < npts) {
    int targpt = workeroffset; // * (0 * threadsperworker)
    real_t mymaxval = targmatdata[targpt];
    for (int x = 1; x < nperchunk; ++x) {
      targpt = (threadsperworker * x) + workeroffset;
      if (targpt < npts) {
        if (targmatdata[targpt] > mymaxval) {
          mymaxidx = targpt;
          mymaxval = targmatdata[targpt];
        }
      }
    }
  }
  *mymaxesptr = mymaxidx; //set mymaxes pointer to my max value (this will be read by root thread)

  //Sync threads before we go single-thread to combine the maxes.
  __syncthreads();

  //Compute final max in root thread
  if (workeroffset == 0) {
    int len = threadsperworker;
    if (npts < threadsperworker) {
      len = npts;
    }

    //I can use mymaxesptr[X] to access the Xth thread's mymax (just pointer arithmetic).
    for (int x = 1; x < len; ++x) {
      int tmpmaxidx = mymaxesptr[x];
      if (targmatdata[tmpmaxidx] > targmatdata[mymaxidx]) {
        mymaxidx = tmpmaxidx;
      }
    }

    //Set root threads to actual max, so all other guys know what it is so all threads can compute do/while condition in parallel
    *mymaxesptr = mymaxidx;
  }
  __syncthreads();

  int myrootthread_maxidx = mymaxesptr[-workeroffset];

  //Returns the GLOBAL max idx always!
  return myrootthread_maxidx;
} //end matrixmax


__device__
int
findmaxloc(real_t* arr, int size) {
  if (size == 0) {
    return -1;
  }
  int maxloc = 0;
  for (int x = 1; x < size; ++x) {
    if (arr[x] > arr[maxloc]) {
      maxloc = x;
    }
  }
  return maxloc;
}


//REV: I need to know where each thread's individual data is...
//REV: Could launch sub-kernels or something...
__device__
void
computeMatch(int* start1ptr, int* end1ptr,
             int* start2ptr, int* end2ptr,
             int* matchlenptr,
             int workeroffset,
             int threadsperworker,
             real_t* seg1, int seg1len, int seg1flatlen,
             real_t* seg2, int seg2len, int seg2flatlen,
             int mfccsize,
             real_t* mydistworkspace, real_t* myqualworkspace,
             int workspace_size,
             int* mymaxesptr,
             int maxmatches,
             const real_t btmin,
             const real_t qThresh,
             const real_t btMinFact,
             const int minlen,
             const real_t penalty,
             const real_t bonus
) {

  int npts = seg1len * seg2len;

  //Compute DISTANCE
  for (int t = workeroffset; t < npts; t += threadsperworker) {
    int seg1t = matrix_x(t, seg1len);
    int seg2t = matrix_y(t, seg1len);

    //Compute Distance
    real_t dot = 0;
    real_t mag1 = 0;
    real_t mag2 = 0;

    //This only tells me workeroffset, I need to compute for all timepoints...

    //Compute all MFCCs for that time point of each guy
    for (int x = 0; x < mfccsize; ++x) {
      int seg1offset = seg1t + (x * seg1len);
      int seg2offset = seg2t + (x * seg2len);
      real_t p1 = seg1[seg1offset];
      real_t p2 = seg2[seg2offset];
      dot += p1 * p2;
      mag1 += p1 * p1;
      mag2 += p2 * p2;
    }
    real_t dist = dot / (sqrt(mag1) * sqrt(mag2));
    mydistworkspace[t] = dist;
  }

  __syncthreads();

  int rslide = seg2len + 1;

  //compute QUALITY
  for (int tNy = 0; (tNy * threadsperworker) < seg2len; ++tNy) {
    for (int t1x = 0; t1x < (seg1len + rslide); ++t1x) //need to capture topright
    {
      int xpos = t1x - workeroffset;
      int ypos = (tNy * threadsperworker) + (workeroffset);

      int t = arrayloc(xpos, ypos, seg1len);

      if (
        xpos <
        0    //REV: Cheat, it's a signed INT so we can actually check for <0 without unsigned overflow error
        ||
        (xpos > seg1len) //or if it is outside of x dir
        ||
        (ypos > seg2len) //or if it is outside of y dir (cannot be negative bc we start from 0 )
        ) {
        //illegal, don't compute...
      } else if (xpos == 0 || ypos == 0) //init qual values to zero on first column and row.
      {
        myqualworkspace[t] = 0;
      } else {
        //REV: To remember which direction it came from, store a 0, 1 or 2 in a matrix of same size ;)
        real_t diag = mydistworkspace[arrayloc((xpos - 1), (ypos - 1), seg1len)];
        real_t horiz = mydistworkspace[arrayloc((xpos - 1), (ypos), seg1len)];
        real_t vert = mydistworkspace[arrayloc((xpos), (ypos - 1), seg1len)];

        real_t arr[3] = {diag, horiz, vert};
        int maxloc = findmaxloc(arr, 3);

        if (maxloc == 0) //diag
        {
          real_t diagqual = myqualworkspace[arrayloc((xpos - 1), (ypos - 1), seg1len)];
          myqualworkspace[t] = diag * bonus + diagqual;
        } else if (maxloc == 1) //horiz
        {
          real_t horizqual = myqualworkspace[arrayloc((xpos - 1), (ypos), seg1len)];
          myqualworkspace[t] = ((1 - horiz) * horizqual) * penalty + horizqual;
        } else if (maxloc == 2) //vert
        {
          real_t vertqual = myqualworkspace[arrayloc((xpos), (ypos - 1), seg1len)];
          myqualworkspace[t] = ((1 - vert) * vertqual) * penalty + vertqual;
        } else {
          printf("REV: MASSIVE SUPER ERROR from findmaxloc\n");
        }

        //If quality would be less than some min, just set to zero (REV: Does this affect quality traceback at all????)
        if (myqualworkspace[t] < btmin) {
          myqualworkspace[t] = 0;
        }

      } //end else if legal to compute (in the for loop?)
      __syncthreads(); //Sync threads to ensure that each "chunk of rows" is done before next set of columns is done since it will use vert.

    } //end for all work to be done.
  }

  //Compute recursive paths...


  if (workeroffset == 0) {
    *matchlenptr = 0;
  }

  int rootmaxidx = matrixmax(workeroffset, threadsperworker, npts,
                             myqualworkspace, mymaxesptr);

  real_t curr_rootmax = myqualworkspace[rootmaxidx];

  while ((curr_rootmax > qThresh) &&
         ((*matchlenptr) < maxmatches)) {

    //REV: do the matches TRACEBACK
    if (workeroffset == 0) {
      real_t matchqual = curr_rootmax;
      real_t bestqual = matchqual;

      int matchcurr = rootmaxidx;

      bool cont = true;
      while (cont) {
        cont = false;

        int xpos = matrix_x(matchcurr, (seg1len));
        int ypos = matrix_y(matchcurr, (seg1len));
        real_t diag = -1000;
        if (xpos > 0 && ypos > 0) {
          diag = mydistworkspace[arrayloc((xpos - 1), (ypos - 1), (seg1len))];
        }
        real_t horiz = -1000;
        if (xpos > 0) { horiz = mydistworkspace[arrayloc((xpos - 1), (ypos), (seg1len))]; }
        real_t vert = -1000;
        if (ypos > 0) { vert = mydistworkspace[arrayloc((xpos), (ypos - 1), (seg1len))]; }

        real_t arr[3] = {diag, horiz, vert};
        int maxloc = findmaxloc(arr, 3);

        if (maxloc == 0) //diag
        {
          matchqual = myqualworkspace[arrayloc((xpos - 1), (ypos - 1), (seg1len))]; //diag
          if (matchqual > (bestqual * btMinFact)) {
            cont = true;
            matchcurr = arrayloc((xpos - 1), (ypos - 1), (seg1len));
          }
        } else if (maxloc == 1) //horiz
        {
          matchqual = myqualworkspace[arrayloc((xpos - 1), (ypos), (seg1len))]; //horiz;
          if (matchqual > (bestqual * btMinFact)) {
            cont = true;
            matchcurr = arrayloc((xpos - 1), (ypos), (seg1len));
          }
        } else //vert
        {
          matchqual = myqualworkspace[arrayloc((xpos), (ypos - 1), (seg1len))]; //vert;
          if (matchqual > (bestqual * btMinFact)) {
            cont = true;
            matchcurr = arrayloc((xpos), (ypos - 1), (seg1len));
          }
        }
      } //end for while continue this alignment

      //Set tmp start/end variables
      int s1start = matrix_x(matchcurr, (seg1len));
      int s1end = matrix_x(rootmaxidx, (seg1len));
      int s2start = matrix_y(matchcurr, (seg1len));
      int s2end = matrix_y(rootmaxidx, (seg1len));


      //REV: Regardless of criterion, blank out the unneeded regions.
      //Do this afterwards so we can do multithreaded...
      int midx = *matchlenptr;

      start1ptr[midx] = s1start;
      end1ptr[midx] = s1end;
      start2ptr[midx] = s2start;
      end2ptr[midx] = s2end;

      (*matchlenptr) += 1;

    } //end if WORKER==0 (ROOT)

    __syncthreads();


    //ERASE BOXES AROUND DISCOVERED ALIGNMENTS (regardless of whether we keep them or not)

    int currmatchidx = (*matchlenptr) - 1;
    if (currmatchidx >=
        0) //sanity check in case there are NO matches at all (i.e. first max is still below qTHRESH)
    {
      real_t USEDVAL = -5.0;
      //REV: This is theoretically more efficient...but while loop isn't ending?
      int idx = workeroffset;
      int matchwidx = (end1ptr[currmatchidx] - start1ptr[currmatchidx]) + 1;
      int matchwidy = (end2ptr[currmatchidx] - start2ptr[currmatchidx]) + 1;
      int idx_x = matrix_x(idx, matchwidx);
      int idx_y = matrix_y(idx, matchwidx);

      while (matchwidx > 0 && matchwidy > 0 &&
             idx_y < matchwidy && idx_x < matchwidx &&
             idx_y >= 0 && idx_x >= 0) {
        int truex = start1ptr[currmatchidx] + idx_x;
        int truey = start2ptr[currmatchidx] + idx_y;
        int t = arrayloc(truex, truey, seg1len);
        myqualworkspace[t] = USEDVAL;
        idx += threadsperworker;
        idx_x = matrix_x(idx, matchwidx);
        idx_y = matrix_y(idx, matchwidx);
      }
    }
    __syncthreads();

    //Now after blanking check whether we should keep it (Min length within each semgent)
    if (workeroffset == 0) {
      if (currmatchidx >= 0) {
        int s1start = start1ptr[currmatchidx];
        int s1end = end1ptr[currmatchidx];
        int s2start = start2ptr[currmatchidx];
        int s2end = end2ptr[currmatchidx];

        int s1len = (s1end - s1start);
        int s2len = (s2end - s2start);
        if (!(s1len > minlen && s2len > minlen)) {
          (*matchlenptr) -= 1;
        }
      }
    }
    __syncthreads(); //need to synch again for while() check involving curr_rootmax;

    //REV: Find maxes again.
    rootmaxidx = matrixmax(workeroffset, threadsperworker, npts,
                           myqualworkspace, mymaxesptr);

    curr_rootmax = myqualworkspace[rootmaxidx];

  } //end WHILE rootmax > qThresh etc.

} //end __device__ void computeMatch()


__global__
void
computeMatches(int* ret_matches_start1s, int* ret_matches_end1s,
               int* ret_matches_start2s, int* ret_matches_end2s,
               int* ret_matches_per,
               int matches_size,
               int nsegments,
               real_t** segdata_ptrs, int* seglen_ptr, int* segflatlen_ptr, int mfccsize,
               int* comp1_ptr, int* comp2_ptr, int ncomps,
               int nworkers,
               real_t** distworkspace_ptrs,
               real_t** qualworkspace_ptrs,
               int workspace_size,
               int* tmpmaxesptr,
               const real_t _btmin,
               const real_t _qThresh,
               const real_t _btminFact,
               const int _minmatchlen,
               const real_t _penalty,
               const real_t _bonus
) {

  //This is #threads I can do at the "same time" (may be scheduled in unknown way).
  int nthreads = gridDim.x * blockDim.x;

  int threadsperworker = nthreads / nworkers;

  int myidx = threadIdx.x + (blockIdx.x * blockDim.x);
  int ngens = intdivceil(ncomps, nworkers); //ceil(ncomps/nworkers);
  int myworker = myidx / threadsperworker; //e.g. 1024 / 64 = 16. In other words I am worker
  int workeroffset = myidx % threadsperworker;

  for (int gen = 0; gen < ngens; ++gen) {
    int startcomp = gen * nworkers;

    if (gen * nworkers + myworker >= ncomps) {
      return; //we are done...
    }

    int compidx = startcomp + myworker;

    int segidx1 = comp1_ptr[compidx];
    int segidx2 = comp2_ptr[compidx];

    real_t* seg1 = segdata_ptrs[segidx1];
    int seg1len = seglen_ptr[segidx1];
    int seg1flatlen = segflatlen_ptr[segidx1];

    real_t* seg2 = segdata_ptrs[segidx2];
    int seg2len = seglen_ptr[segidx2];
    int seg2flatlen = segflatlen_ptr[segidx2];

    real_t* mydistworkspace = distworkspace_ptrs[myworker];
    real_t* myqualworkspace = qualworkspace_ptrs[myworker];

    int moffset = compidx * matches_size;
    int* mymatchstart1ptr = &(ret_matches_start1s[moffset]);
    int* mymatchstart2ptr = &(ret_matches_start2s[moffset]);
    int* mymatchend1ptr = &(ret_matches_end1s[moffset]);
    int* mymatchend2ptr = &(ret_matches_end2s[moffset]);
    int* mymatchlenptr = &(ret_matches_per[compidx]);

    int* mymaxesptr = &(tmpmaxesptr[myidx]);

    computeMatch(mymatchstart1ptr, mymatchend1ptr,
                 mymatchstart2ptr, mymatchend2ptr,
                 mymatchlenptr,
                 workeroffset,
                 threadsperworker,
                 seg1, seg1len, seg1flatlen,
                 seg2, seg2len, seg2flatlen,
                 mfccsize,
                 mydistworkspace, myqualworkspace,
                 workspace_size,
                 mymaxesptr,
                 matches_size,
                 _btmin,
                 _qThresh,
                 _btminFact,
                 _minmatchlen,
                 _penalty,
                 _bonus
    );

  }
}

//REV: I could compute these indices on the GPU to determine which segments to compare at each point
//Returns COL and ROW of it.
//See testtriangle.py for explanation of algebra.

//REV: Initially
/*
unsigned int row_index_w_diag( unsigned int i, unsigned int M )
{
  double m = M;
  double row = (-2*m - 1 + sqrt( (4*m*(m+1) - 8*(double)i - 7) )) / -2;
  if( row == (double)(int) row ) { row -= 1; }
  return (unsigned int) row;
}


unsigned int column_index_w_diag( unsigned int i, unsigned int M )
{
  unsigned int row = row_index_w_diag( i, M);
  return  i - M * row + row*(row+1) / 2;
}
*/

unsigned int row_index_w_diag(unsigned int i, unsigned int M) {
  double m = M;
  double ii = (m * (m + 1)) / 2 - 1 - i;
  unsigned int K = floor((sqrt(8 * ii + 1) - 1) / 2);
  //if( row == (double)(int) row ) { row -= 1; }
  return (M - 1 - K);
}


unsigned int column_index_w_diag(unsigned int i, unsigned int M) {
  double m = M;
  double ii = (m * (m + 1)) / 2 - 1 - i;
  unsigned int K = floor((sqrt(8 * ii + 1) - 1) / 2);
  unsigned int jj = (i - M * (M + 1) / 2 + (K + 1) * (K + 2) / 2);
  return jj + row_index_w_diag(i, M);
}


unsigned int row_index(unsigned int i, unsigned int M) {
  double m = M;
  double ii = (m * (m - 1)) / 2 - 1 - i;
  unsigned int K = floor((sqrt(8 * ii + 1) + 1) / 2);
  //if( row == (double)(int) row ) { row -= 1; }
  return (M - 1 - K);
}


unsigned int column_index(unsigned int i, unsigned int M) {
  double m = M;
  double ii = (m * (m - 1)) / 2 - 1 - i;
  unsigned int K = floor((sqrt(8 * ii + 1) + 1) / 2);
  unsigned int jj = (i - M * (M - 1) / 2 + (K + 1) * (K + 2) / 2);
  return jj;
}


std::vector <std::vector<std::vector < Match>>>

findSubSequencesCUDA(const std::vector <segment> &segments,
                     const bool &usesmooth,
                     const bool &usecosine,
                     const int &maxnmatches,
                     const real_t &_btmin,
                     const real_t &_qThresh,
                     const real_t &_btminFact,
                     const int &_minmatchlen,
                     const real_t &_penalty,
                     const real_t &_bonus,
                     const int &cudadev /*=0*/
) {
  cudaSetDevice(cudadev);


  int mfcclen = 39;

  int nworkers = 128;
  int threadsperblock = 128;
  int nthreads = nworkers * threadsperblock;
  int blockspergrid = nthreads / threadsperblock;


  std::vector <std::vector<bool>> donematches;
  std::vector < std::vector < std::vector < Match>> > matches;

  if (segments.size() > 0) {
    if (segments[0].smooth.size() > 0) {
      mfcclen = segments[0].smooth[0].size();
    } else {
      fprintf(stderr,
              "REV: ERROR: first segment is of length zero, can't determine correct MFCC size\n");
      exit(1);
    }
  } else {
    fprintf(stderr, "REV: WARNING: no segments passed for computation\n");
    return matches;
  }

  int nsegments = segments.size();

  matches.resize(nsegments);
  donematches.resize(nsegments);

  for (size_t s = 0; s < nsegments; ++s) {
    matches[s].resize((nsegments - s));
    donematches[s].resize((nsegments - s),
                          false); //-1, because even for 0th guy, we NEVER do vs ourself...
  }


  //REV: 26 Mar 2016: Modify code to compute matches iteratively too, to reduce amount of in-memory at a time
  //This is "comparisons" and also the "output" of each... assume 4 bytes each, times 186 million? Still only 800 million bytes, i.e. 800e6, 800 MB?
  //Is it the segments that are taking up the memory? The MFCCs? Yea, I flatten them I guess.....so doubling the space used...

  //Best we can do is to only do s1starts, etc. iteratively, and derive matches from that
  //Same for comparisons
  //s1matches, etc.
  //We could also delete segments data before we go, but that might cause problems if something else wants to use it.


  //I could delete flatsegments after copying to GPU. I guess that works.

  //I need to shrink_to_fit or something...since .clear() will not free memory. Or force it to go out of scope...


  //Only top triangle with** diagnal
  int ncomps = ((nsegments + 1) * nsegments) / 2;

  fprintf(stdout, "Will perform [%d] comparisons on [%d] segments\n", ncomps, nsegments);

  int longestlength = 0;

  std::vector < real_t * > d_seg_ptrs; //(flatsegments.size());
  real_t** d_seg_ptrs_ptr;
  int* d_len_ptr;
  int* d_flatlen_ptr;

  std::vector < real_t * > d_distworkspace_ptrs(nworkers);
  std::vector < real_t * > d_qualworkspace_ptrs(nworkers);
  real_t** d_distworkspace_ptrs_ptr;
  real_t** d_qualworkspace_ptrs_ptr;

  int workspacesize;// = longestlength * longestlength;
  int workspace_sidesize;// = longestlength;
  int* d_tmpmaxes_ptr; //This value is allocated (once), then I pass the ptr location to everybody...



  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  //REV: Local scope so that flatsegments are deallocated after copying to GPU.
  {

    //REV: Flatten the 2d segments.
    std::vector <std::vector<real_t>> flatsegments(segments.size());
    std::vector<int> seglengths(segments.size());
    std::vector<int> flatlengths(segments.size());

    d_seg_ptrs.resize(flatsegments.size());

    //Find max size among all segments, and flatten the segments into a 1d array, with adjacent time points next to each other (i.e. MFCCsize time-length
    //vectors are stacked.)
    for (int s = 0; s < segments.size(); ++s) {
      int seglen = segments[s].smooth.size();
      if (seglen > longestlength) {
        longestlength = seglen;
      }

      if (seglen > 0) {
        if (segments[s].smooth[0].size() != mfcclen) {
          fprintf(stderr, "ERROR MFCC len not as expected\n");
          exit(1);
        }

        std::vector <real_t> flatseg(seglen * mfcclen);

        for (int tx = 0; tx < mfcclen; ++tx) {
          int tpoint = tx * seglen;
          for (int t = 0; t < seglen; ++t) {
            if (segments[s].smooth[t].size() != mfcclen) {
              fprintf(stderr,
                      "Found a timepoint [%d] of seg [%d] in smooth dataset where MFCC vector is not expected length [%d] (is [%ld] instead)\n",
                      t, s, mfcclen, segments[s].smooth[t].size());
              exit(1);
            }

            flatseg[tpoint + t] = segments[s].smooth[t][tx];
          }
        }

        flatlengths[s] = flatseg.size();
        flatsegments[s] = flatseg;

      } else {
        fprintf(stderr,
                "WARNING: segment length is 0!!! Seg [%d]. Not doing anything (no comparisons should be made with this segment)\n",
                s);
      }

      seglengths[s] = seglen;
    } //end for all segments s



    workspacesize = longestlength * longestlength;
    workspace_sidesize = longestlength;

    //Do allocations in GPU
    for (int s = 0; s < flatsegments.size(); ++s) {
      if (flatsegments[s].size() > 0) {
        checkCudaErrors(
          cudaMalloc(&(d_seg_ptrs[s]), flatsegments[s].size() * sizeof(real_t)));

        checkCudaErrors(cudaMemcpy(d_seg_ptrs[s], flatsegments[s].data(),
                                   flatsegments[s].size() * sizeof(flatsegments[s][0]),
                                   cudaMemcpyHostToDevice));
      }
    }

    //Literally copying (vector of) pointers...
    checkCudaErrors(cudaMalloc(&d_seg_ptrs_ptr, d_seg_ptrs.size() * sizeof(d_seg_ptrs[0])));
    checkCudaErrors(
      cudaMemcpy(d_seg_ptrs_ptr, d_seg_ptrs.data(), d_seg_ptrs.size() * sizeof(d_seg_ptrs[0]),
                 cudaMemcpyHostToDevice));


    checkCudaErrors(cudaMalloc(&d_tmpmaxes_ptr, nthreads * sizeof(int)));


    //Just malloc them
    for (int w = 0; w < nworkers; ++w) {
      checkCudaErrors(cudaMalloc(&(d_distworkspace_ptrs[w]), workspacesize * sizeof(real_t)));
      checkCudaErrors(cudaMalloc(&(d_qualworkspace_ptrs[w]), workspacesize * sizeof(real_t)));
    }

    //Literally copying (vector of) pointers...
    checkCudaErrors(cudaMalloc(&d_distworkspace_ptrs_ptr,
                               d_distworkspace_ptrs.size() * sizeof(d_distworkspace_ptrs[0])));
    checkCudaErrors(cudaMemcpy(d_distworkspace_ptrs_ptr, d_distworkspace_ptrs.data(),
                               d_distworkspace_ptrs.size() * sizeof(d_distworkspace_ptrs[0]),
                               cudaMemcpyHostToDevice));

    //Literally copying (vector of) pointers...
    checkCudaErrors(cudaMalloc(&d_qualworkspace_ptrs_ptr,
                               d_qualworkspace_ptrs.size() * sizeof(d_qualworkspace_ptrs[0])));
    checkCudaErrors(cudaMemcpy(d_qualworkspace_ptrs_ptr, d_qualworkspace_ptrs.data(),
                               d_qualworkspace_ptrs.size() * sizeof(d_qualworkspace_ptrs[0]),
                               cudaMemcpyHostToDevice));


    //copy lengths
    checkCudaErrors(cudaMalloc(&d_len_ptr, seglengths.size() * sizeof(seglengths[0])));
    checkCudaErrors(
      cudaMemcpy(d_len_ptr, seglengths.data(), seglengths.size() * sizeof(seglengths[0]),
                 cudaMemcpyHostToDevice));

    //copy flat lengths
    checkCudaErrors(cudaMalloc(&d_flatlen_ptr, flatlengths.size() * sizeof(flatlengths[0])));
    checkCudaErrors(
      cudaMemcpy(d_flatlen_ptr, flatlengths.data(), flatlengths.size() * sizeof(flatlengths[0]),
                 cudaMemcpyHostToDevice));


  }

  //end local scope for deallocation purposes
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////



  int compsperchunk = nworkers * 1000; //128*100 = 128k

  if (ncomps < compsperchunk) {
    compsperchunk = ncomps;
  }


  int nchunks = intdivceil(ncomps, compsperchunk);

  //These have only size related to compsperchunk
  int* d_comp1_ptr;
  int* d_comp2_ptr;

  //copy comp1 vect and comp2 vect
  checkCudaErrors(cudaMalloc(&d_comp1_ptr, compsperchunk * sizeof(int)));
  checkCudaErrors(cudaMalloc(&d_comp2_ptr, compsperchunk * sizeof(int)));


  //Pointers to hold results
  int* d_start1matches_ptr;
  int* d_end1matches_ptr;
  int* d_start2matches_ptr;
  int* d_end2matches_ptr;
  int* d_nmatches_ptr;

  int matchessizeperchunk = maxnmatches * compsperchunk;
  //int matchessize = maxnmatches * ncomps;

  checkCudaErrors(cudaMalloc(&d_start1matches_ptr, matchessizeperchunk * sizeof(int)));
  checkCudaErrors(cudaMalloc(&d_end1matches_ptr, matchessizeperchunk * sizeof(int)));
  checkCudaErrors(cudaMalloc(&d_start2matches_ptr, matchessizeperchunk * sizeof(int)));
  checkCudaErrors(cudaMalloc(&d_end2matches_ptr, matchessizeperchunk * sizeof(int)));
  checkCudaErrors(cudaMalloc(&d_nmatches_ptr, compsperchunk * sizeof(int)));


  //REV: These are output vectors. These are now smaller because they were taking 4GB each.
  std::vector<int> s1matches(matchessizeperchunk);
  std::vector<int> e1matches(matchessizeperchunk);
  std::vector<int> s2matches(matchessizeperchunk);
  std::vector<int> e2matches(matchessizeperchunk);
  std::vector<int> nmatches(compsperchunk);

  std::vector<int> compseg1(compsperchunk);
  std::vector<int> compseg2(compsperchunk);

  //Main LOOP

  //Do all comparisons in "chunks", to minimize unnecessary matches allocation.
  //while( totalcomps < ncomps )
  for (int chunki = 0; chunki < nchunks; ++chunki) {
    compseg1.resize(compsperchunk);
    compseg2.resize(compsperchunk);


    //How many will we actually do (skip zero-length segment-containing comparisons)
    unsigned int cnt = 0;

    unsigned int startcomp = chunki * compsperchunk;
    unsigned int endcomp = startcomp + compsperchunk;

    //Shorten it
    if (endcomp > ncomps) {
      endcomp = ncomps; //REV: Is this right? We index comp from zero, so yes.
    }

    for (unsigned int x = startcomp; x < endcomp; ++x) {
      //Compute row and column (i.e. s1 and s2 index)
      unsigned int s1 = row_index_w_diag(x, nsegments);
      unsigned int s2 = column_index_w_diag(x, nsegments);
      if (s1 >= nsegments || s2 >= nsegments) {
        fprintf(stderr,
                "REV: Big error in segment index computation...got s1=%u  s2=%u despite only %d segments\n",
                s1, s2, nsegments);
        exit(1);
      }

      //REV: "skip" comparisons with zero-legnth s1 or s2
      //REV: SKIP comparisons of s1==s2
      if (segments[s1].smooth.size() > 0 && segments[s2].smooth.size() > 0 &&
          (s1 != s2)) {
        compseg1[cnt] = s1;
        compseg2[cnt] = s2;
        ++cnt;
      }
    }

    //Shorten vector size in case there are less than compsperchunk
    compseg1.resize(cnt);
    compseg2.resize(cnt);

    int compsize = cnt;

    if (compseg1.size() < compsize || compseg2.size() < compsize) {
      fprintf(stderr, "REV: ERROR compseg1 etc. have weird size (<compsize, which is [%d])!!\n",
              compsize);
      exit(1);
    }

    fprintf(stdout, "REV: Doing a chunk [%d] of [%d] chunks ([%d] comps per chunk)\n", (chunki + 1),
            nchunks, compsperchunk);

    //Copy this "turn's" comparison chunk that will be computed in this kernel call
    checkCudaErrors(cudaMemcpy(d_comp1_ptr, compseg1.data(), compsize * sizeof(compseg1[0]),
                               cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_comp2_ptr, compseg2.data(), compsize * sizeof(compseg2[0]),
                               cudaMemcpyHostToDevice));


    computeMatches<< < threadsperblock, blockspergrid>> > (d_start1matches_ptr, d_end1matches_ptr,
      d_start2matches_ptr, d_end2matches_ptr,
      d_nmatches_ptr,
      maxnmatches,
      nsegments,
      d_seg_ptrs_ptr,
      d_len_ptr,
      d_flatlen_ptr,
      mfcclen,
      d_comp1_ptr,
      d_comp2_ptr,
      compsize,
      nworkers,
      d_distworkspace_ptrs_ptr,
      d_qualworkspace_ptrs_ptr,
      workspace_sidesize,
      d_tmpmaxes_ptr,
      _btmin,
      _qThresh,
      _btminFact,
      _minmatchlen,
      _penalty,
      _bonus
    );

    checkCudaErrors(cudaMemcpy(s1matches.data(),
                               d_start1matches_ptr,
                               (compsize * maxnmatches) * sizeof(s1matches[0]),
                               cudaMemcpyDeviceToHost)
    );

    checkCudaErrors(cudaMemcpy(e1matches.data(),
                               d_end1matches_ptr,
                               (compsize * maxnmatches) * sizeof(e1matches[0]),
                               cudaMemcpyDeviceToHost)
    );

    checkCudaErrors(cudaMemcpy(s2matches.data(),
                               d_start2matches_ptr,
                               (compsize * maxnmatches) * sizeof(s2matches[0]),
                               cudaMemcpyDeviceToHost)
    );

    checkCudaErrors(cudaMemcpy(e2matches.data(),
                               d_end2matches_ptr,
                               (compsize * maxnmatches) * sizeof(e2matches[0]),
                               cudaMemcpyDeviceToHost)
    );

    checkCudaErrors(cudaMemcpy(nmatches.data(),
                               d_nmatches_ptr,
                               (compsize) * sizeof(nmatches[0]),
                               cudaMemcpyDeviceToHost)
    );


    cudaDeviceSynchronize();

    //REV: matches is allocated to nsegments. Thus, just push back to each s1 segment.

    unsigned int cnt2 = 0;
    //Do *all* that would be in this one, specifically startcomp to endcomp
    for (unsigned int comp = startcomp; comp < endcomp; ++comp) {
      unsigned int s1 = row_index_w_diag(comp, nsegments);
      unsigned int s2 = column_index_w_diag(comp, nsegments);

#if DEBUGLEVEL > 10
      if( s1 == 0 )
        {
          if( s2 == 0 )
      {
        fprintf(stdout, "GPU sanity check of row/col idxs: [0][0] is idx [%d] (expect 0)\n", comp);
      }

          if( s2 == 1 )
      {
        fprintf(stdout, "GPU sanity check of row/col idxs: [0][1] is idx [%d] (expect 1)\n", comp);
      }
        }
#endif //DEBUGLEVEL > 10
      if (s1 >= nsegments || s2 >= nsegments) {
        fprintf(stderr,
                "REV: Big error in segment index computation (getting matches)...got s1=%u  s2=%u despite only %d segments\n",
                s1, s2, nsegments);
        exit(1);
      }

      std::vector <Match> mymatches;

      if (segments[s1].smooth.size() > 0 && segments[s2].smooth.size() > 0 &&
          (s1 != s2)) {
        int matchesstart =
          maxnmatches * cnt2; //REV: becuase we only put in GPU the ones with >0 length...
        int mynmatches = nmatches[cnt2];

        for (int x = 0; x < mynmatches; ++x) {
          Match m;

          //BHO: scaled these  values to their position in seconds in the files
          int idx = matchesstart + x;
          m.s1.start = 160 * s1matches[idx] + segments[s1].start;
          m.s2.start = 160 * s2matches[idx] + segments[s2].start;
          m.s1.end = 160 * e1matches[idx] + segments[s1].start;
          m.s2.end = 160 * e2matches[idx] + segments[s2].start;
          m.s1.fname = segments[s1].fname;
          m.s2.fname = segments[s2].fname;
          mymatches.push_back(m);
        }

        ++cnt2;

      } //only iter if seg1>0 && seg2>0

      //size_t curr_s1_match=matches[s1].size();
      //if( curr_s1_match != s2 )
      if (matches[s1].size() <= s2 - s1) {
        //fprintf(stderr, "REV (WARNING): In ordering of S1/S2 matches: S2 is [%d], but we were expecting match [%ld]\n", s2, curr_s1_match);
        fprintf(stderr,
                "REV (WARNING): In ordering of S1/S2 matches: S2 is [%d] (literal %d), but S1=[%d] triangle column only has [%ld] size\n",
                s2, s2 - s1, s1, matches[s1].size());
        exit(1);
      }

      if (donematches[s1][s2 - s1] == true) {
        fprintf(stderr, "REV: ERROR, redoing %d %d\n", s1, s2);
        exit(1);
      }

      matches[s1][s2 -
                  s1] = mymatches; //we still push back an (empty) matches vector even if length was 0, to keep with brad's order.

      donematches[s1][s2 - s1] = true;
    } //end for all comps I did this chunk

  } //end for all chunks to compute

  std::vector <size_t> notdones1;
  std::vector <size_t> notdones2;
  for (size_t s1 = 0; s1 < matches.size(); ++s1) {
    for (size_t s2 = 0; s2 < matches[s1].size(); ++s2) {
      if (donematches[s1][s2] == false) {
        notdones1.push_back(s1);
        notdones2.push_back(s2);
        fprintf(stdout, "NOT DONE!! %ld %ld\n", s1, s2);
      }
    }
  }
  if (notdones1.size() > 0) {
    exit(1);
  }


  fprintf(stdout, "FINISHED KERNEL\n");

  /*
  cnt = 0;
  for (int s1 = 0; s1 < nsegments; ++s1) {
    std::vector <std::vector<Match>> smatches;
    for (int s2 = (s1 + 1); s2 < nsegments; ++s2) //REV: Not doing same guys (for now)
    {
      std::vector <Match> mymatches;

      if (segments[s1].smooth.size() > 0 && segments[s2].smooth.size() > 0) {
        int matchesstart = maxnmatches * cnt;
        int mynmatches = nmatches[cnt];

        //fprintf(stdout, "For comparison [%d] ([%d] to [%d]), found [%d] matches\n", cnt, s1, s2, mynmatches);
        for (int x = 0; x < mynmatches; ++x) {
          Match m;

          //temp.quality = best.q;
          //temp.s1.start = alignment[alignment.size() - 1][0];
          //temp.s1.end = alignment[0][0]; //best.pi;
          //temp.s1.fname = novel.fname;
          //temp.s2.start = alignment[alignment.size() - 1][1];
          //temp.s2.end = alignment[0][1];//best.pj;
          //temp.s2.fname = base.fname;
          //temp.classification = -1;

          //BHO: scaled these  values to their position in seconds in the files
          int idx = matchesstart + x;
          m.s1.start = 160 * s1matches[idx] + segments[s1].start;
          m.s2.start = 160 * s2matches[idx] + segments[s2].start;
          m.s1.end = 160 * e1matches[idx] + segments[s1].start;
          m.s2.end = 160 * e2matches[idx] + segments[s2].start;
          //m.s1.start = s1matches[idx];
          //m.s2.start = s2matches[idx];
          //m.s1.end = e1matches[idx];
          //m.s2.end = e2matches[idx];
          m.s1.fname = segments[s1].fname;
          m.s2.fname = segments[s2].fname;
          //REV: Crap: Need to readout quality...
          mymatches.push_back(m);

          //fprintf(stdout, "Match start1 [%d] end1 [%d]   start2 [%d] end2 [%d]\n", m.start1, m.end1, m.start2, m.end2 );
        }

        ++cnt;
      }
      smatches.push_back(mymatches);
    }
    matches.push_back(smatches);
  }
  */

  //Cleanup
  for (int s = 0; s < segments.size(); ++s) {
    if (segments[s].smooth.size() > 0) {
      checkCudaErrors(cudaFree(d_seg_ptrs[s]));
    }
  }
  checkCudaErrors(cudaFree(d_seg_ptrs_ptr));
  checkCudaErrors(cudaFree(d_tmpmaxes_ptr));
  for (int w = 0; w < nworkers; ++w) {
    checkCudaErrors(cudaFree(d_distworkspace_ptrs[w]));
    checkCudaErrors(cudaFree(d_qualworkspace_ptrs[w]));
  }

  checkCudaErrors(cudaFree(d_distworkspace_ptrs_ptr));
  checkCudaErrors(cudaFree(d_qualworkspace_ptrs_ptr));
  checkCudaErrors(cudaFree(d_len_ptr));
  checkCudaErrors(cudaFree(d_flatlen_ptr));

  checkCudaErrors(cudaFree(d_comp1_ptr));
  checkCudaErrors(cudaFree(d_comp2_ptr));

  checkCudaErrors(cudaFree(d_start1matches_ptr));
  checkCudaErrors(cudaFree(d_end1matches_ptr));
  checkCudaErrors(cudaFree(d_start2matches_ptr));
  checkCudaErrors(cudaFree(d_end2matches_ptr));
  checkCudaErrors(cudaFree(d_nmatches_ptr));

  //fprintf(stdout, "FINISHED READING OUT MATCHES\n");
  return matches;
}
