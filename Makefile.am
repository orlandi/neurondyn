#CXX = clang++
CXX = g++-4.9
AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}
LIBS = $(DEPS_LIBS)
AM_CPPFLAGS = $(DEPS_CFLAGS)
bin_PROGRAMS = neurondyn
INCLUDEPATH = . src
neurondyn_SOURCES = src/adaptationrunner.cc src/adaptationrunner.h src/circularvector.h src/main.cc src/main.h src/netdyn.cc src/netdyn.h
neurondyn_LDADD = -lgsl -lgslcblas -lconfig++
neurondyn_LDFLAGS = -fopenmp -L/Users/orlandi/projects/libconfig/lib
neurondyn_CPPFLAGS = -I/Users/orlandi/projects/libconfig/include
neurondyn_CPPFLAGS += -std=c++11 -fopenmp
