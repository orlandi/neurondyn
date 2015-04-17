/*
 * Copyright (c) 2009-2015 Javier G. Orlandi <javierorlandi@javierorlandi.com>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#pragma once
#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
#include "circularvector.h"
#include "adaptationrunner.h"
#include <libconfig.h++>
#include "gsl/gsl_statistics.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

#ifdef CUDA_ENABLED
#include <cuda.h> 
#include <curand.h>
#define CUDA_CALL(x) do { if((x)!=cudaSuccess) { \
                          printf("Error at %s:%d\n",__FILE__,__LINE__);\
                          return EXIT_FAILURE;}} while(0)
#define CURAND_CALL(x) do { if((x)!=CURAND_STATUS_SUCCESS) { \
                            printf("Error at %s:%d\n",__FILE__,__LINE__);\
                            return EXIT_FAILURE;}} while(0)
#define CUDA_CALL_DEST(x) do { if((x)!=cudaSuccess) { \
                          printf("Error at %s:%d\n",__FILE__,__LINE__);\
                          }} while(0)
#define CURAND_CALL_DEST(x) do { if((x)!=CURAND_STATUS_SUCCESS) { \
                            printf("Error at %s:%d\n",__FILE__,__LINE__);\
                            }} while(0)
#endif

#include <math.h>
#include <iostream>
#include <vector>

enum neurotransmitterType
{
  NT_AMPA = 0x01,
  NT_NMDA = 0x02,
  NT_GABA = 0x04
};
enum neuronClass
{
  NEU_RS = 0x01,
  NEU_LS = 0x02,
  NEU_FS = 0x04
};
enum deltaType
{
    LIN = 0,
    LOG = 1
};

#define STD_DATA_COUNT 1000000

class NetDyn
{
  public:
  NetDyn();
  ~NetDyn();
  void setConnectivity(std::string filename);
  void printConnectivityMap();
  void configureNeuron();
  void configureNeuronTypes(double* cs, double* ds, double* ps, int* ts, int num);
  void setDefaultNeuronTypes();
  void configureSimulation();
  void loadConfigFile(std::string filename);
  template <typename T>
  bool assignConfigValue(const char* entry, T& configVariable, bool critical = true);
  void recordSpike(int i, int st);
  void configureSpikeRecord(int maxsize = STD_DATA_COUNT, std::string filename = "SpikeRecord.txt", double subset = 0.05);
  void configureTraces(int num, std::string filename, double sampTime = 1);
  void processSpikes(int size = -1);
  void recordTraces();
  void processTraces(int size = -1);
  bool seedRng();
  int simulationStart();
  bool simulationRun();
  void loadNeuronType(std::string filename);
  inline bool simulationIsRunning() { return running; }
  inline bool simulationIsActive() { return active; }
  inline int getOriginalRngSeed() { return originalRngSeed; }

  bool simulationStep();

  void simulationFinish();

  void initialize();
  void setConstants();
  void setExternalConnections(std::string filename);
  void setPositions(std::string filename);
  void saveNetworkStructure(std::string filename = "NetworkStructure.txt");
  void initSaveResults(std::string tmpStr, std::string filename = "");
  void saveResults(std::string tmpStr, std::string filename = "");

  bool burstCheck();

  double getNormalRng();
  double getUniformRng();

#ifdef OPENMP_ENABLED
  bool seedParallelRng(int seeds);
#endif

  private:
  std::string savedFileName, configFileName;
  gsl_rng* rng; // RNG structure
  libconfig::Config* configFile;

#ifdef CUDA_ENABLED
  bool cuda;
  curandGenerator_t cudaRng;
  float *devNormalRngData, *hostNormalRngData;
  float *devUniformRngData, *hostUniformRngData;
  size_t nextNormalRngIndex, nextUniformRngIndex, cudaRngChunkSize;
#endif

  bool running, active, multiplicativeMini, depressionMini;
#ifdef OPENMP_ENABLED
  bool parallel;
  gsl_rng** parallelRng;
  int* parallelThreadUsage;
#endif
  double* miniStrength, *miniTime;

  int* connectivityMap;
  int* connectivityNumber;
  int* inputConnectivityNumber, *numberExcitatoryInputs, *numberInhibitoryInputs;
  int* connectivityFirstIndex;
  int nNumber, totalConnections;
  int nType, numTypes;
  // Neuron variables
  double* C, *K, *V_r, *V_t, *V_p, *a, *b, *p;
  double* v, *u, *I, *I_GABA, *I_AMPA, *I_NMDA, *I_WNOISE, *c, *d, *D, *meanD_AMPA, *meanD_GABA, *v_d;
  double* positionX, *positionY;
  int* neurotransmitter, *neuronalClass;
  int* ntTypes, *neuronType, neuronTypes;
  double g_AMPA, g_NMDA, g_GABA, g_mAMPA, g_mGABA, g_WNOISE;
  double tau_NMDA, tau_AMPA, tau_GABA, tau_D_AMPA, tau_D_GABA, beta_AMPA, beta_GABA, tau_MINI;
  int totalSpikes;

  // Simulation parameters
  double dt, t, totalTime;
  int simulationSteps, step;

  bool GABA, NMDA, AMPA, MINI, WNOISE, depression;

  double exp_AMPA, exp_NMDA, exp_GABA, dt_times_a, dt_over_tau_D_AMPA, dt_over_tau_D_GABA, dt_over_C, strength_WNOISE;
  // Post processing
  int* dSpikeNeuron, dSpikeRecord, dSpikeRecordLimit, *dSpikeStep, dSpikeSubsetSize;
  std::string dSpikesFile, traceFile, dSpikesSubsetFile;
  std::string resultsFolder, presetTypesFile;
  bool tracing, presetTypes;
  int* tracedNeuron, numberTraces, traceSamplingTime, traceRecord;
  double** traceV, **traceU, **traceI, **traceI_AMPA, **traceI_NMDA, **traceI_GABA;
  double** traceI_WNOISE, **traceD, *traceTime, **traceG_AMPA, **traceG_GABA;
  char currentPath[FILENAME_MAX];
  int simulationReturnValue;
  int originalRngSeed;

  // Burst detector parameters
  bool burstDetector;
  double burstDetectorBinSize, burstDetectorMinimumSpikesPerNeuron;
  int burstDetectorBinNumber, burstDetectorStepsPerBin;
  int burstDetectorFirstStepAboveThreshold, burstDetectorLastStepAboveThreshold;
  double burstDetectorSpikeCount;
  bool burstDetectorPossibleBurst;
  std::string burstDetectorFile;
  // Burst detector internal variables
  circularVector* burstDetectorMemory;
  std::vector<double> burstDetectorStorage;

  // Burst transition exploration parameters
  bool burstTransitionExploration;
  double burstTransitionExplorationDelta, burstTransitionExplorationMaximumIBI;
  int burstTransitionExplorationRequiredBursts;
  int burstTransitionExplorationDeltaType;

  // Burst adaptation parameters
  bool adaptiveIBI;
  double adaptiveIBItotalTime;
  double adaptiveIBItarget, adaptiveIBIaccuracy, adaptiveIBIbaseMultiplier, adaptiveIBIextrapolationMultiplier;
  double adaptiveIBIweightLowerBound, adaptiveIBIweightUpperBound;
  int adaptiveIBImaxIterations;
  // Burst adaptation internal variables
  adaptationRunner* adaptiveIBIrunner;
  int adaptiveIBIcurrentIteration, adaptiveIBItotalTimeSteps;
  bool adaptiveIBIfinished;

  bool dryRun;
};
