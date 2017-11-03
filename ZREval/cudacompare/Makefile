CC = g++
INCLDIR=./include
IMPLDIR=./impl
CPPDIR=./cpp
CPPFLAG=-std=c++11 -O3 #-g -Wall
INCLFLAG=-I. -I$(INCLDIR) -I$(IMPLDIR)
LIBFLAG=-L. -L/usr/lib/x86_64-linux-gnu/
LIBS=
CUDALIBS=-lcuda -lcudart
SHAREDLIBFLAG=-fPIC
CUDALD=-L/usr/local/cuda-7.5/targets/x86_64-linux/lib
NVCC=/usr/local/cuda/bin/nvcc
#GENCODE=-gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 -gencode arch=compute_37,code=sm_37 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_52,code=sm_52
GENCODE=-gencode arch=compute_37,code=sm_37
NVCCFLAGS=-m64 -g -G -DDEBUG=6 $(GENCODE)

INSTALL_PATH=/usr/local

srcs = $(wildcard $(CPPDIR)/*.cpp)
objs = $(srcs:.cpp=.o)
deps = $(srcs:.cpp=.d)

cudacompare.exe: cudacompare.cu main.cpp
	$(NVCC) $(CUDALD) $(NVCCFLAGS) $(CPPFLAG) $(INCLFLAG) $(LIBFLAG) -c cudacompare.cu -lcuda -lcudart
	$(CC) $(CUDALD) $(CPPFLAG) $(INCLFLAG) $(LIBFLAG) main.cpp cudacompare.o -o cudacompare.exe $(LIBS) $(CUDALIBS)

%.o: %.cpp
	$(MPICC) $(SHAREDLIBFLAG) -MMD -MP -c $(CPPFLAG) $(INCLFLAG) $< -o $@

.PHONY: clean

# $(RM) is rm -f by default
clean:
	$(RM) $(objs) $(deps) test

-include $(deps)

