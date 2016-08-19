#CXX = clang++
CXX = g++
AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4
LIBS = $(DEPS_LIBS)
AM_CPPFLAGS = $(DEPS_CFLAGS)
bin_PROGRAMS = neurondyn
INCLUDEPATH = . src
neurondyn_SOURCES = src/adaptationrunner.cc src/adaptationrunner.h src/circularvector.h src/main.cc src/main.h src/netdyn.cc src/netdyn.h
neurondyn_LDADD = -lgsl -lgslcblas -lconfig++
neurondyn_LDFLAGS = -L${HOME}/lib
neurondyn_CPPFLAGS = -I${HOME}/include -std=c++11

if CUDA_ENABLED
  neurondyn_LDADD += -lcurand -lcudart -lpthread -ldl -L/usr/local/cuda/lib64
  neurondyn_CPPFLAGS += -I/usr/local/cuda/include
endif

if OPENMP_ENABLED
  neurondyn_LDFLAGS += -fopenmp
  neurondyn_CPPFLAGS += -fopenmp
endif
