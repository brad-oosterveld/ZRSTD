## Zero Resource Spoken Term Discovery

This is the code used for the evaluation in our ICASSP 2017 paper [A PARALLELIZED DYNAMIC PROGRAMMING APPROACH TO ZERO RESOURCE SPOKEN TERM DISCOVERY](https://hrilab.tufts.edu/publications/oosterveld2017icassp.pdf)

If you want to cite it here is the Bibtex:
~~~~~
@InProceedings{oosterveld2017icassp,
author={Bradley Oosterveld and Richard Veale and Matthias Scheutz},
title={A Parallelized Dynamic Programming Approach to Zero Resource Spoken Term Discovery},
booktitle={Proceedings of the 42nd IEEE International Conference on Acoustics, Speech, and Signal Processing},
year={2017},
project={osl},
}
~~~~~

If you are interested in our other work it can be found here: https://hrilab.tufts.edu/publications/

### 0. Notes:

    This system has only been tested on Ubuntu 14.04. 

    I do not have enough time to properly maintain this code, so it is provided "as is". That being said if you find a bug and fix it, send me a pull request, and I can merge it in.

### 1. Installation:

   1. Dependencies:
      * cmake https://cmake.org/
        ~~~~
        $ sudo apt-get install cmake
        ~~~~

      * libarmadillo-dev http://arma.sourceforge.net/
        ~~~~
        $ sudo apt-get install libarmadillo-dev
        ~~~~
      
      * cuda
      
        follow instructions at:  https://developer.nvidia.com/cuda-downloads
        
        CUDA 5.0 or greater required

   2. Building
      ~~~~
      $ mkdir build
      $ cd build
      $ cmake ..
      $ make
      ~~~~

### 2. Running

   To run the code on the data as provided by the [Zero Speech 2015](http://sapience.dec.ens.fr/bootphon/2015/index.html)

   ~~~~
   $ ./build/bin/ZREval <path_to_data_directory>
   ~~~~

   There is some functionality built in to have the program perform differently for different languages, but it doesn't do anything currently.
   
   As described in the paper evaluation took 346 minutes to generate the complete transcription of the English corpus which is approximately 634 minutes long, and 67 minutes for 244 minutes of Tsonga data. These runtimes were recorded on a computer using 16 GB RAM, Intel Core i7-3820 CPU with 8 cores at 3.60GHz, and NVidia Titan X GPU with 3072 CUDA cores at 1.0GHz and 12 GB VRAM. The parallelized version of our implementation produced a 3200× speedup over the serial CPU implementation.
   
   Things may be faster slower/for you depending on the machine you're using. The GPU is the main factor that affects performance, but if you run out of RAM and start to use swap things will also slow dow a ton.
