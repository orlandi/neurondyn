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

#include "netdyn.h"
#include <sys/time.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <limits>

#include <cmath>
#include <ctime>

#ifdef OPENMP_ENABLED
#include <omp.h>
#endif

NetDyn::NetDyn()
{
  rng = NULL;
  connectivityMap = NULL;
  connectivityNumber = NULL;
  connectivityFirstIndex = NULL;
  burstDetectorMemory = NULL;

  running = false;
  active = false;
  tracing = false;
#ifdef OPENMP_ENABLED
  parallel = false;
  parallelRng = NULL;
  parallelThreadUsage = NULL;
#endif
#ifdef CUDA_ENABLED
  cuda = false;
#endif
  presetTypes = false;
  burstDetector = false;
  adaptiveIBI = false;
  burstTransitionExploration = false;
  dryRun = false;
  WNOISE = false;
  AMPA = false;
  NMDA = false;
  GABA = false;
  MINI = false;
  depression = false;
  simulationReturnValue = 0;

  resultsFolder = "data";
}

NetDyn::~NetDyn()
{
  gsl_rng_free(rng);

#ifdef OPENMP_ENABLED
  if(parallel)
  {
    for (int i = 0; i < omp_get_max_threads(); i++)
    {
      gsl_rng_free(parallelRng[i]);
    }
    delete [] parallelRng;
  }
#endif

#ifdef CUDA_ENABLED
  if(cuda)
  {
    CURAND_CALL_DEST(curandDestroyGenerator(cudaRng));
    CUDA_CALL_DEST(cudaFree(devNormalRngData));
    CUDA_CALL_DEST(cudaFree(devUniformRngData));
    free(hostNormalRngData);
    free(hostUniformRngData);
  }
#endif

  // TODO: free everything else
}

void NetDyn::setConnectivity(std::string filename)
{
  std::ifstream inputFile;
  std::string line;
  int nFrom, nTo, nNeurons, nConnections, nCurrent;
  int connectionIndex;
  std::stringstream tmpStr;
  nNeurons = 0;
  nConnections = 0;
  inputFile.open(filename.c_str(), std::ifstream::in);
  if (!inputFile.is_open())
  {
    std::cerr << "There was an error opening the file " << filename << "\n";
    exit(1);
    return;
  }

  // First loop through the file to get number of synapses and neurons
  while (std::getline(inputFile, line))
  {
    if (line.at(0) == '%')
      continue;
    tmpStr.clear();
    tmpStr.str(line);
    tmpStr >> nFrom >> nTo;
    //        std:: cout << line << " " << nFrom << "\n";
    if (nTo > nNeurons)
      nNeurons = nTo;
    nConnections++;
  }
  if (nFrom > nNeurons)
    nNeurons = nFrom;
  nNeurons++;
  std::cout << "Neurons: " << nNeurons << " Connections: " << nConnections << "\n";
  totalConnections = nConnections;
  nNumber = nNeurons;

  connectivityMap = new int[totalConnections];
  connectivityNumber = new int[nNumber];
  for (int i = 0; i < nNumber; i++)
    connectivityNumber[i] = 0;
  connectivityFirstIndex = new int[nNumber];

  inputFile.clear();
  inputFile.seekg(0, std::ios::beg);
  nCurrent = 0;
  nConnections = 0;
  connectionIndex = 0;
  connectivityFirstIndex[0] = 0;

  std::cout << "Assigning output connectivity...\n";
  while (std::getline(inputFile, line))
  {
    if (line.at(0) == '%')
      continue;
    tmpStr.clear();
    tmpStr.str(line);
    tmpStr >> nFrom >> nTo;

    connectivityMap[connectionIndex] = nTo;
    if (nFrom != nCurrent)
    {
      connectivityNumber[nCurrent] = nConnections;
      nCurrent = nFrom;
      nConnections = 1;
      connectivityFirstIndex[nCurrent] = connectionIndex;
    }
    else
      nConnections++;
    connectionIndex++;
  }
  // For the last neuron
  connectivityNumber[nCurrent] = nConnections;

  // For the input connections
  std::cout << "Assigning input connectivity...\n";

  inputConnectivityNumber = new int[nNumber];
  for (int i = 0; i < nNumber; i++)
    inputConnectivityNumber[i] = 0;
  for (int i = 0; i < totalConnections; i++)
    inputConnectivityNumber[connectivityMap[i]]++;
}

int NetDyn::simulationStart()
{
  std::cout << "Starting the simulation...\n";

  time_t tstart, tend;
  tstart = time(0);
  bool done = false;
  while(!done)
  { 
      done = simulationStep();
  }
  std::cout << "Simulation finished!\n";
  tend = time(0);
  std::cout << "It took " << difftime(tend, tstart) << " second(s).\n";
#ifdef OPENMP_ENABLED
  if (parallel)
  {
    std::cout << "Global thread usage: ";
    for (int i = 0; i < omp_get_max_threads(); i++)
      std::cout << parallelThreadUsage[i] << " ";
    std::cout << "\n";
  }
#endif
  processSpikes(dSpikeRecord);
  if (tracing)
    processTraces(traceRecord);

  if (chdir(currentPath) == -1)
  {
    std::cerr << "Error: " << strerror(errno) << std::endl;
    std::cerr << currentPath << std::endl;
    exit(1);
  }

  return simulationReturnValue;
}

void NetDyn::printConnectivityMap()
{
  std::cout << "Printing connectivity map...\n";
  for (int i = 0; i < nNumber; i++)
  {
    std::cout << i << " (" << connectivityNumber[i] << ") --> ";
    for (int j = connectivityFirstIndex[i]; j < connectivityFirstIndex[i] + connectivityNumber[i]; j++)
      std::cout << connectivityMap[j] << " ";
    std::cout << "\n";
  }
}

void NetDyn::loadNeuronType(std::string filename)
{
  std::ifstream inputFile;
  inputFile.open(filename.c_str(), std::ifstream::in);
  if (!inputFile.is_open())
  {
    std::cerr << "There was an error opening the neuronType file: " << filename << "\n";
    std::cerr << "Using probabilities instead.\n";
    return;
  }
  std::string line;
  int nnum, ntyp;
  std::istringstream tmpStr;

  while (std::getline(inputFile, line))
  {
    if (line.at(0) == '%')
      continue;
    tmpStr.clear();
    tmpStr.str(line);
    tmpStr >> nnum >> ntyp;
    neuronType[nnum] = ntyp;
  }
  std::cout << "neuronType file loaded.\n";
  //    std::cout << "Types: ";
  //    for(int i = 0;i < nNumber; i++)
  //        std::cout << neuronType[i] << " ";
  //    std::cout << "\n";
}

void NetDyn::setPositions(std::string filename)
{
  std::ifstream inputFile;
  inputFile.open(filename.c_str(), std::ifstream::in);
  if (!inputFile.is_open())
  {
    std::cerr << "There was an error opening the file " << filename << "\n";
    exit(1);
  }
  std::string line, nNum, posX, posY;
  int nnum;
  double posx, posy;
  std::istringstream tmpStr;
  positionX = new double[nNumber];
  positionY = new double[nNumber];

  while (std::getline(inputFile, line))
  {
    tmpStr.clear();
    tmpStr.str(line);
    tmpStr >> nnum >> posx >> posy;
    positionX[nnum] = posx;
    positionY[nnum] = posy;
  }
}

void NetDyn::loadConfigFile(std::string filename)
{
  int* nt;
  std::string ntString, tmpStr, lookupStr;
  configFileName = filename;

  configFile = new libconfig::Config();
  try
  {
    configFile->readFile(filename.c_str());
  }
  catch (libconfig::FileIOException& e)
  {
    std::cerr << "Main error. Config file could not be loaded.\n";
    exit(1);
  }

  configFile->setAutoConvert(true);

  // Start setting the connectivity
  assignConfigValue("connectivity.file", tmpStr);
  setConnectivity(tmpStr);

  // Load the positions
  assignConfigValue("positions.file", tmpStr);
  setPositions(tmpStr);

  assignConfigValue("results.path", resultsFolder, false);

  // Configure neuron models TODO: error checking for this part
  libconfig::Setting& typelist = configFile->lookup("neuron.models");
  neuronTypes = typelist.getLength();
  a = new double[neuronTypes];
  b = new double[neuronTypes];
  C = new double[neuronTypes];
  K = new double[neuronTypes];
  V_r = new double[neuronTypes];
  V_t = new double[neuronTypes];
  V_p = new double[neuronTypes];
  c = new double[neuronTypes];
  d = new double[neuronTypes];
  p = new double[neuronTypes];
  nt = new int[neuronTypes];
  neurotransmitter = new int[neuronTypes];
  neuronalClass = new int[neuronTypes];
  for (int i = 0; i < neuronTypes; i++)
  {
    typelist[i].lookupValue("a", a[i]);
    typelist[i].lookupValue("b", b[i]);
    typelist[i].lookupValue("C", C[i]);
    typelist[i].lookupValue("K", K[i]);
    typelist[i].lookupValue("V_r", V_r[i]);
    typelist[i].lookupValue("V_t", V_t[i]);
    typelist[i].lookupValue("V_p", V_p[i]);
    typelist[i].lookupValue("c", c[i]);
    typelist[i].lookupValue("d", d[i]);

    typelist[i].lookupValue("p", p[i]);
    typelist[i].lookupValue("neurotransmitters", ntString);

    nt[i] = 0;
    if (ntString.find("AMPA") != std::string::npos)
      nt[i] = nt[i] | NT_AMPA;
    if (ntString.find("NMDA") != std::string::npos)
      nt[i] = nt[i] | NT_NMDA;
    if (ntString.find("GABA") != std::string::npos)
      nt[i] = nt[i] | NT_GABA;
    neurotransmitter[i] = nt[i];

    typelist[i].lookupValue("class", tmpStr);
    neuronalClass[i] = NEU_RS;
    if (tmpStr.find("RS") != std::string::npos)
      neuronalClass[i] = NEU_RS;
    else if (tmpStr.find("LS") != std::string::npos)
      neuronalClass[i] = NEU_LS;
    else if (tmpStr.find("FS") != std::string::npos)
      neuronalClass[i] = NEU_FS;
  }

  // Configure remaining parameters
  assignConfigValue("simulation.timestep", dt);
  assignConfigValue("simulation.totalTime", totalTime);
#ifdef OPENMP_ENABLED
  assignConfigValue("simulation.openmp", parallel, false);
#endif
#ifdef CUDA_ENABLED
  assignConfigValue("simulation.cuda", cuda, false);
  if(cuda)
  {
    assignConfigValue("simulation.cuda_rng_chunk_size", cudaRngChunkSize);
  }
#endif
  // Configure soma

  assignConfigValue("neuron.soma.whitenoise.active", WNOISE, false);
  if(WNOISE)
    assignConfigValue("neuron.soma.whitenoise.strength", g_WNOISE);

  // Configure synapse

  assignConfigValue("neuron.synapse.AMPA.active", AMPA, false);
  if(AMPA)
  {
    assignConfigValue("neuron.synapse.AMPA.g", g_AMPA);
    assignConfigValue("neuron.synapse.AMPA.tau", tau_AMPA);
  }

  assignConfigValue("neuron.synapse.NMDA.active", NMDA, false);
  if(NMDA)
  {
    assignConfigValue("neuron.synapse.NMDA.g", g_NMDA);
    assignConfigValue("neuron.synapse.NMDA.tau", tau_NMDA);
  }

  assignConfigValue("neuron.synapse.GABA.active", GABA, false);
  if(GABA)
  {
    assignConfigValue("neuron.synapse.GABA.g", g_GABA);
    assignConfigValue("neuron.synapse.GABA.tau", tau_GABA);
  }

  assignConfigValue("neuron.synapse.minis.active", MINI);
  if(MINI)
  {
    assignConfigValue("neuron.synapse.minis.g_AMPA", g_mAMPA);
    assignConfigValue("neuron.synapse.minis.g_GABA", g_mGABA);
    assignConfigValue("neuron.synapse.minis.tau", tau_MINI);
    assignConfigValue("neuron.synapse.minis.multiplicative", multiplicativeMini);
    assignConfigValue("neuron.synapse.minis.depression", depressionMini);
  }

  assignConfigValue("neuron.synapse.depression.active", depression, false);
  if(depression)
  {
    assignConfigValue("neuron.synapse.depression.beta_AMPA", beta_AMPA);
    assignConfigValue("neuron.synapse.depression.tau_AMPA", tau_D_AMPA);
    assignConfigValue("neuron.synapse.depression.beta_GABA", beta_GABA);
    assignConfigValue("neuron.synapse.depression.tau_GABA", tau_D_GABA);
  }

  // Configure Preset types
  assignConfigValue("neuron.presetTypes.active", presetTypes, false);
  if (presetTypes)
    assignConfigValue("neuron.presetTypes.file", presetTypesFile);

  // Configure the protocols

  // burst_detector
  assignConfigValue("protocols.burst_detector.active", burstDetector, false);
  if (burstDetector)
  {
    assignConfigValue("protocols.burst_detector.bin_size", burstDetectorBinSize);
    assignConfigValue("protocols.burst_detector.bin_number", burstDetectorBinNumber);
    assignConfigValue("protocols.burst_detector.minimum_spikes_per_neuron", burstDetectorMinimumSpikesPerNeuron);
    assignConfigValue("protocols.burst_detector.file", burstDetectorFile);
  }

  // adaptive_IBI
  assignConfigValue("protocols.adaptive_IBI.active", adaptiveIBI, false);
  if (adaptiveIBI)
  {
    if(!burstDetector)
    {
        std::cerr << "Error! burst_detector is required for adaptive_IBI." << std::endl;
        exit(1);
    }
    assignConfigValue("protocols.adaptive_IBI.adaptation_time", adaptiveIBItotalTime);
    assignConfigValue("protocols.adaptive_IBI.target_IBI", adaptiveIBItarget);
    assignConfigValue("protocols.adaptive_IBI.accuracy", adaptiveIBIaccuracy);
    assignConfigValue("protocols.adaptive_IBI.base_multiplier", adaptiveIBIbaseMultiplier);
    assignConfigValue("protocols.adaptive_IBI.extrapolation_multiplier", adaptiveIBIextrapolationMultiplier);
    assignConfigValue("protocols.adaptive_IBI.weight_lower_bound", adaptiveIBIweightLowerBound);
    assignConfigValue("protocols.adaptive_IBI.weight_upper_bound", adaptiveIBIweightUpperBound);
    assignConfigValue("protocols.adaptive_IBI.max_iterations", adaptiveIBImaxIterations);
  }

  assignConfigValue("protocols.burst_transition_exploration.active", burstTransitionExploration, false);
  if (burstTransitionExploration)
  {
    if(!burstDetector)
    {
      std::cerr << "Error! burst_detector is required for burst_transition_exploration." << std::endl;
      exit(1);
    }
    assignConfigValue("protocols.burst_transition_exploration.delta", burstTransitionExplorationDelta);
    assignConfigValue("protocols.burst_transition_exploration.delta_type", tmpStr);
    if (tmpStr.find("LIN") != std::string::npos)
    {
      burstTransitionExplorationDeltaType = LIN;
    }
    else if (tmpStr.find("LOG") != std::string::npos)
    {
      burstTransitionExplorationDeltaType = LOG;
    }
    else
    {
        std::cerr << "Unrecognized delta_type. Only (LIN|LOG) are allowed." << "\n";
        exit(1);
    }
    assignConfigValue("protocols.burst_transition_exploration.required_bursts", burstTransitionExplorationRequiredBursts);
    assignConfigValue("protocols.burst_transition_exploration.maximum_IBI", burstTransitionExplorationMaximumIBI);
  }
  if(adaptiveIBI && burstTransitionExploration)
  {
    std::cerr << "Error! adaptive_IBI and burst_transition_exploration are incompatible." << std::endl;
    exit(1);
  }

  // Initialize the rest
  seedRng();
#ifdef OPENMP_ENABLED
  if(parallel)
  {
    omp_set_num_threads(omp_get_max_threads());
    std::cout << "Number of available threads: " << omp_get_max_threads() << "\n";
    seedParallelRng(omp_get_max_threads());
    parallelThreadUsage = new int[omp_get_max_threads()];
    for (int i = 0; i < omp_get_max_threads(); i++)
      parallelThreadUsage[i] = 0;
  }
  else
    omp_set_num_threads(1);
#endif

  // Initialize variables
  initialize();

  setConstants();

  // Configure recordings
  int maxsize = STD_DATA_COUNT;
  int ntraces;
  double subset = 0.05;
  double sampling;
  tmpStr = "spikeRecord.txt";

  assignConfigValue("logging.spikes.spikespersave", maxsize);
  assignConfigValue("logging.spikes.file", tmpStr);
  assignConfigValue("logging.spikes.subset", subset);
  configureSpikeRecord(maxsize, tmpStr, subset);

  assignConfigValue("logging.traces.number", ntraces);
  assignConfigValue("logging.traces.file", tmpStr);
  assignConfigValue("logging.traces.sampling", sampling);
  configureTraces(ntraces, tmpStr, sampling);
}

template <typename T>
bool NetDyn::assignConfigValue(const char* entry, T& configVariable, bool critical)
{
  if (!configFile->lookupValue(entry, configVariable))
  {
    if (critical)
    {
      std::cerr << "Critical Error! " << entry << " not defined in config file.\n";
      exit(1);
    }
    else
      std::cerr << "Warning! " << entry << " not defined in config file.\n";
    return false;
  }
  else
    return true;
}

void NetDyn::setExternalConnections(std::string filename)
{
  std::ifstream inputFile;
  inputFile.open(filename.c_str(), std::ifstream::in);
  if (!inputFile.is_open())
  {
    std::cerr << "There was an error opening the file " << filename << "\n";
    exit(1);
  }
  std::string line;
  int nnum, econs;
  std::istringstream tmpStr;

  while (std::getline(inputFile, line))
  {
    tmpStr.clear();
    tmpStr.str(line);
    tmpStr >> nnum >> econs;
    numberExcitatoryInputs[nnum] = econs;
  }
}

void NetDyn::initialize()
{
  v = new double[nNumber];
  u = new double[nNumber];
  I = new double[nNumber];
  v_d = new double[nNumber];

  neuronType = new int[nNumber];

  if (AMPA)
  {
    I_AMPA = new double[nNumber];
  }
  if (NMDA)
  {
    I_NMDA = new double[nNumber];
  }
  if (GABA)
  {
    I_GABA = new double[nNumber];
  }
  if (WNOISE)
  {
    I_WNOISE = new double[nNumber];
  }
  if (depression)
  {
    D = new double[nNumber];
    meanD_AMPA = new double[nNumber];
    meanD_GABA = new double[nNumber];
  }

  // Preallocate to define the neuron types
  gsl_ran_discrete_t* randist;
  randist = gsl_ran_discrete_preproc(neuronTypes, p);

  // PRELOAD neurontype in here
  if (presetTypes)
  {
    loadNeuronType(presetTypesFile);
  }

  for (int i = 0; i < nNumber; i++)
  {
    u[i] = 0.;
    I[i] = 0.;

    // Assign the neuron type from the type distribution
    if (!presetTypes)
      neuronType[i] = gsl_ran_discrete(rng, randist);

    v[i] = V_r[neuronType[i]];
    v_d[i] = v[i];

    if (AMPA)
      I_AMPA[i] = 0.;
    if (GABA)
      I_GABA[i] = 0.;
    if (NMDA)
      I_NMDA[i] = 0.;
    if (WNOISE)
      I_WNOISE[i] = 0.;
    D[i] = 1.;
  }

  numberExcitatoryInputs = new int[nNumber];
  numberInhibitoryInputs = new int[nNumber];
  for (int i = 0; i < nNumber; i++)
  {
    numberExcitatoryInputs[i] = 0;
    numberInhibitoryInputs[i] = 0;
  }
  int k, neighIndex, ntype;

  for (int i = 0; i < nNumber; i++)
  {
    ntype = neuronType[i];
    k = connectivityFirstIndex[i];
    //    std::cout << connectivityNumber[i] << "\n";
    for (int j = 0; j < connectivityNumber[i]; j++)
    {
      neighIndex = connectivityMap[k + j];

      if ((neurotransmitter[ntype] & NT_AMPA) || (neurotransmitter[ntype] & NT_NMDA))
        numberExcitatoryInputs[neighIndex]++;
      if ((neurotransmitter[ntype] & NT_GABA))
        numberInhibitoryInputs[neighIndex]++;
    }
  }
  saveNetworkStructure();
  simulationSteps = ceil(totalTime / dt);
  step = 0;
  totalSpikes = 0;

  if (burstDetector)
  {
    burstDetectorFirstStepAboveThreshold = 0;
    burstDetectorLastStepAboveThreshold = 0;
    burstDetectorSpikeCount = 0;
    burstDetectorPossibleBurst = false;
    burstDetectorStepsPerBin = floor(burstDetectorBinSize / dt);
    if (fmod(burstDetectorBinSize / dt, 1.0) != 0)
      std::cerr << "Warning! BinSize is not a multiple of dt\n";
    burstDetectorMemory = new circularVector(burstDetectorBinNumber);
    //burstDetectorStorage = new std::vector<double>;
    burstDetectorStorage.clear();
    initSaveResults("% Storing burst times (in s)\n", burstDetectorFile);
  }

  if (adaptiveIBI)
  {
    adaptiveIBIcurrentIteration = 0;
    adaptiveIBItotalTimeSteps = floor(adaptiveIBItotalTime / dt);
    adaptiveIBIfinished = false;
    dryRun = true;
    std::cout << "dryRun is on. Not storing spikes or traces.\n";
    adaptiveIBIrunner = new adaptationRunner(adaptiveIBIbaseMultiplier, adaptiveIBIextrapolationMultiplier, false, true);
    adaptiveIBIrunner->setLowerBound(adaptiveIBIweightLowerBound);
    adaptiveIBIrunner->setUpperBound(adaptiveIBIweightUpperBound);
  }

  if(burstTransitionExploration)
  {
    totalTime = 2.*burstTransitionExplorationMaximumIBI*burstTransitionExplorationRequiredBursts;
    simulationSteps = ceil(totalTime / dt);
    dryRun = true;
  }

}

void NetDyn::saveNetworkStructure(std::string filename)
{
  std::stringstream tmpStr;
  int nt = 0;

  tmpStr
      << "%-------------------------------------------------------------\n"
      << "% Neurondyn \n"
      << "% Copyright (c) 2009-2015 Javier G. Orlandi <javierorlandi@javierorlandi.com> \n"
      << "%-------------------------------------------------------------\n"
      << "% Neuron idx | Type (0 Ex, 1 In) | k_i | k_i(E) | k_i(I) | k_o\n";

  for (int i = 0; i < nNumber; i++)
  {
    if (neurotransmitter[neuronType[i]] & NT_GABA)
      nt = 1;
    else
      nt = 0;
    tmpStr << i << " " << nt << " "
           << inputConnectivityNumber[i] << " "
           << numberExcitatoryInputs[i] << " "
           << numberInhibitoryInputs[i] << " "
           << connectivityNumber[i] << "\n";
    if (i % 1000 == 0)
    {
      saveResults(tmpStr.str(), filename);
      tmpStr.clear();
      tmpStr.str("");
    }
  }
  saveResults(tmpStr.str(), filename);
}

bool NetDyn::simulationStep()
{
  bool returnValue;
  returnValue = false;

  if (tracing)
    if (step % traceSamplingTime == 0)
      recordTraces();

  // Reset the mean depression in the terminals
  if (depression && depressionMini)
  {
    for (int i = 0; i < nNumber; i++)
    {
      meanD_AMPA[i] = 0;
      meanD_GABA[i] = 0;
    }
  }

  // Check for spikes
  for (int i = 0; i < nNumber; i++)
  {
    int ntype = neuronType[i];
    // If theres a spike, reset and propagate
    if (v[i] > V_p[ntype])
    {
      // Soma reset
      v[i] = c[ntype];
      u[i] += d[ntype];
      // Propagate currents to the neighbors;
      int k = connectivityFirstIndex[i];
      //#pragma omp parallel for private(neighIndex) schedule(static, 100)
      for (int j = 0; j < connectivityNumber[i]; j++)
      {
        int neighIndex = connectivityMap[k + j];

        if ((neurotransmitter[ntype] & NT_AMPA) && AMPA)
        {
          I_AMPA[neighIndex] += g_AMPA * D[i];
        }
        if ((neurotransmitter[ntype] & NT_NMDA) && NMDA)
          I_NMDA[neighIndex] += g_NMDA * D[i];
        if ((neurotransmitter[ntype] & NT_GABA) && GABA)
        {
          I_GABA[neighIndex] += g_GABA * D[i];
        }
      }
      // Neurotransmitter depletion
      if (depression)
      {
        if ((neurotransmitter[ntype] & NT_AMPA) && AMPA)
          D[i] *= beta_AMPA;
        if ((neurotransmitter[ntype] & NT_GABA) && GABA)
          D[i] *= beta_GABA;
      }

      // Now store the spike
      recordSpike(i, step);
      totalSpikes++;
      if (burstDetector)
      {
        burstDetectorMemory->increaseCurrentBin();
      }
    }
    // Update the mean depression: assign to each neuron a mean depression of its input connections
    if (depression && depressionMini)
    {
      int k = connectivityFirstIndex[i];
      for (int j = 0; j < connectivityNumber[i]; j++)
      {
        int neighIndex = connectivityMap[k + j];

        if ((neurotransmitter[ntype] & NT_AMPA) && AMPA)
          meanD_AMPA[neighIndex] += D[i];
        if ((neurotransmitter[ntype] & NT_GABA) && GABA)
          meanD_GABA[neighIndex] += D[i];
      }
    }
  }
  // Correctly assign the mean input depression  
  if (depression && depressionMini)
  {
    for (int i = 0; i < nNumber; i++)
    {
      meanD_AMPA[i] = meanD_AMPA[i] / numberExcitatoryInputs[i];
      meanD_GABA[i] = meanD_GABA[i] / numberInhibitoryInputs[i];
    }
  }
// Now we sum the contribution of the different currents
// if a spike happened in the previous loop, it means that it happened in between
// the previous time step, so currents should be updated accordingly. And they are!
// Also, we can update them in the same step
  // Trying to split the loop for vectorization
  #pragma omp parallel for if(parallel)
  for (int i = 0; i < nNumber; i++)
  {
#ifdef OPENMP_ENABLED
    if (parallel)
      parallelThreadUsage[omp_get_thread_num()]++;
#endif
    I[i] = 0;

    // Now add the noise
    if (WNOISE)
    {
      I_WNOISE[i] = strength_WNOISE * getNormalRng();
    }
    // MINIS HERE - only AMPA and GABA
    if (MINI)
    {
      double tmpMiniConstant, tmpMiniStrength;
      if (AMPA)
      {
        if (multiplicativeMini)
          tmpMiniConstant = double(numberExcitatoryInputs[i]) * dt / tau_MINI;
        else
          tmpMiniConstant = dt / tau_MINI;

        if (depressionMini)
          tmpMiniStrength = meanD_AMPA[i] * g_mAMPA;
        else
          tmpMiniStrength = g_mAMPA;

        if(getUniformRng() < tmpMiniConstant)
            I_AMPA[i] += tmpMiniStrength;
      }
      if (GABA)
      {
        if (multiplicativeMini)
          tmpMiniConstant = double(numberInhibitoryInputs[i]) * dt / tau_MINI;
        else
          tmpMiniConstant = dt / tau_MINI;

        if (depressionMini)
          tmpMiniStrength = meanD_GABA[i] * g_mGABA;
        else
          tmpMiniStrength = g_mGABA;

        if(getUniformRng() < tmpMiniConstant)
            I_GABA[i] += tmpMiniStrength;
      }
    }
  }
  for (int i = 0; i < nNumber; i++)
  {
    // Start to add all the currents
    if (WNOISE)
      I[i] += I_WNOISE[i];
    if (AMPA)
    {
      I[i] += I_AMPA[i];
      I_AMPA[i] *= exp_AMPA;
    }
    if (NMDA)
    {
      I[i] += I_NMDA[i];
      I_NMDA[i] *= exp_NMDA;
    }
    if (GABA)
    {
      I[i] += I_GABA[i];
      I_GABA[i] *= exp_GABA;
    }
  }
  for (int i = 0; i < nNumber; i++)
  {
    int ntype = neuronType[i];
    // Now it is time to update the soma - Euler for now

    // Update the membrane potential
    double vtmp = 0.0;
    if (neuronalClass[ntype] & (NEU_RS | NEU_FS))
      vtmp = v[i] + dt / C[ntype] * (K[ntype] * (v[i] - V_r[ntype]) * (v[i] - V_t[ntype]) - u[i] + I[i]);
    // Here there's a dendtric compartment
    else if (neuronalClass[ntype] & NEU_LS)
    {
      vtmp = v[i] + dt / C[ntype] * (K[ntype] * (v[i] - V_r[ntype]) * (v[i] - V_t[ntype]) + 1.2 * (v_d[i] - v[i]) - u[i] + I[i]);
      v_d[i] = v_d[i] + 0.01 * dt * (v[i] - v_d[i]);
    }
    // Update the slow variable
    if (neuronalClass[ntype] & (NEU_RS | NEU_LS))
      u[i] += dt * a[ntype] * (b[ntype] * (v[i] - V_r[ntype]) - u[i]);
    else if (neuronalClass[ntype] & NEU_FS)
    {
      if (v[i] < V_r[ntype])
        u[i] += -dt * a[ntype] * u[i];
      else
        u[i] += dt * a[ntype] * (b[ntype] * (v[i] - V_r[ntype]) * (v[i] - V_r[ntype]) * (v[i] - V_r[ntype]) - u[i]);
    }
    v[i] = vtmp;
  }
  for (int i = 0; i < nNumber; i++)
  {
    if (depression)
    {
      int ntype = neuronType[i];
      if ((neurotransmitter[ntype] & NT_AMPA) && AMPA)
        D[i] += dt_over_tau_D_AMPA * (1. - D[i]);
      if ((neurotransmitter[ntype] & NT_GABA) && GABA)
        D[i] += dt_over_tau_D_GABA * (1. - D[i]);
    }
  }
  // Finished iterating over neurons

  // Now check for bursts
  if (burstDetector)
  {
    if (step > 1 & step % burstDetectorStepsPerBin == 0)
    {
      // If true, there was a burst
      if(burstCheck())
      {
        std::stringstream tmpStr;
        tmpStr.precision(15);

        if(burstTransitionExploration)
        {
          // Store synaptic strength and burst time
          tmpStr << g_AMPA << " " << burstDetectorStorage.back() << "\n";
          saveResults(tmpStr.str(), burstDetectorFile);

          // Check if we had enough bursts
          int numBursts = burstDetectorStorage.size();
          if(numBursts >= burstTransitionExplorationRequiredBursts)
          {
            // Check for end condition
            double meanIBI = (burstDetectorStorage.at(numBursts - 1) - burstDetectorStorage.at(1)) / (double(numBursts) - 2.0);
            if(meanIBI >= burstTransitionExplorationMaximumIBI)
            {
              returnValue = true;
            }
            // Else proceed to next strength
            else
            {
              if(burstTransitionExplorationDeltaType == LIN)
              {
                g_AMPA += burstTransitionExplorationDelta;
              }
              else if(burstTransitionExplorationDeltaType == LOG)
              {
                g_AMPA = pow(10.0,log10(g_AMPA)+burstTransitionExplorationDelta);
              }
              std::cout << "Finished g iteration with an <IBI> of: " << meanIBI << ", new g_AMPA: " << g_AMPA << "\n\n";

              // Reset the system
              step = 0;
              burstDetectorStorage.clear();
              burstDetectorMemory->clear();
            }
          }
        }
        else
        {
          tmpStr << burstDetectorStorage.back() << "\n";
          saveResults(tmpStr.str(), burstDetectorFile);
        }
      }
      if(step > 1)
        burstDetectorMemory->moveToNextBin();
    }
  }
  if(returnValue)
    return returnValue;

  // Move to the next step
  step++;
  if (step % int(10000 / dt) == 0)
    processSpikes(dSpikeRecord);

  if (adaptiveIBI)
  {
    if ((step % adaptiveIBItotalTimeSteps) == 0)
    {
      int numBursts = burstDetectorStorage.size();
      double meanIBI;
      std::cout << "Adaptive iteration finished. Found " << numBursts << " bursts.\n";
      if (numBursts < 3)
      {
        std::cout << "Not enough bursts to do anything. Setting mean IBI to infinity.\n";
        meanIBI = std::numeric_limits<float>::infinity();
      }
      else
      {
        meanIBI = (burstDetectorStorage.at(numBursts - 1) - burstDetectorStorage.at(1)) / (double(numBursts) - 2.0);
        std::cout << "Removing the first burst... Mean IBI: " << meanIBI << " s\n";
      }
      if (fabs(meanIBI - adaptiveIBItarget) <= adaptiveIBIaccuracy)
      {
        std::cout << "Obtained desired IBI at g_AMPA: " << g_AMPA << " . Continuing with the simulation.\n";
        adaptiveIBI = false;
        dryRun = false;
        std::cout << "dryRun off. Storing spikes and/or traces\n";
        step = 0;
      }
      else
      {
        adaptiveIBIrunner->addDataPair(meanIBI, g_AMPA);
        g_AMPA = adaptiveIBIrunner->getPrediction(adaptiveIBItarget);
        std::cout << "New g_AMPA: " << g_AMPA << "\n";
        // Reset the burst counter
        burstDetectorStorage.clear();
        burstDetectorMemory->clear();
        step = 0;
      }
    }
  }

  if (step >= simulationSteps)
    returnValue = true;

  return returnValue;
}

void NetDyn::setConstants()
{
  if (AMPA)
    exp_AMPA = exp(-dt / tau_AMPA);
  if (NMDA)
    exp_NMDA = exp(-dt / tau_NMDA);
  if (GABA)
    exp_GABA = exp(-dt / tau_GABA);
  if (depression)
  {
    dt_over_tau_D_AMPA = dt / tau_D_AMPA;
    dt_over_tau_D_GABA = dt / tau_D_GABA;
  }
  if (WNOISE)
    strength_WNOISE = sqrt(2. * g_WNOISE / dt);
}

bool NetDyn::seedRng()
{
  std::stringstream tmpStr, seedStr;
  struct timeval tv;
  gettimeofday(&tv, NULL);
  struct tm* tm = localtime(&tv.tv_sec);
  int seed = abs(int(tv.tv_usec / 10 + tv.tv_sec * 100000)); // Creates the seed based on actual time
  originalRngSeed = seed;
  rng = gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(rng, seed); // Seeds the previously created RNG
  seedStr << seed << ".txt";
  savedFileName = seedStr.str();

  tmpStr << "% Date: " << tm->tm_mday << "/" << tm->tm_mon + 1 << "/" << tm->tm_year + 1900 << ", "
         << "Time: " << tm->tm_hour << ":" << tm->tm_min << ":" << tm->tm_sec << ", "
         << "RNG Seed: " << seed << "\n";

#ifdef CUDA_ENABLED
  tmpStr << "% CUDA enabled at compile time.\n";
#else
  tmpStr << "% CUDA not enabled at compile time.\n";
#endif

#ifdef OPENMP_ENABLED
  tmpStr << "% OpenMP enabled at compile time.\n";
#else
  tmpStr << "% OpenMP not enabled at compile time.\n";
#endif         

  std::cout << tmpStr.str();

#ifdef CUDA_ENABLED
  if(cuda)
  {
    // Create the RNG and allocate memory
    /* Allocate n floats on host */
    hostNormalRngData = (float*)calloc(cudaRngChunkSize, sizeof(float));
    hostUniformRngData = (float*)calloc(cudaRngChunkSize, sizeof(float));
    /* Allocate n floats on device */
    CUDA_CALL(cudaMalloc((void**)&devNormalRngData, cudaRngChunkSize * sizeof(float)));
    CUDA_CALL(cudaMalloc((void**)&devUniformRngData, cudaRngChunkSize * sizeof(float)));
    /* Create pseudo-random number generator */
    CURAND_CALL(curandCreateGenerator(&cudaRng, CURAND_RNG_PSEUDO_MTGP32));
    /* Set seed */
    CURAND_CALL(curandSetPseudoRandomGeneratorSeed(cudaRng, seed));

    // Generate the first batch of rngs
    /* Generate n floats on device sampled from a normal distribution */
    CURAND_CALL(curandGenerateNormal(cudaRng, devNormalRngData, cudaRngChunkSize, 0.0, 1.0));
    /* Copy device memory to host */
    CUDA_CALL(cudaMemcpy(hostNormalRngData, devNormalRngData, cudaRngChunkSize * sizeof(float), cudaMemcpyDeviceToHost));
    /* Generate n floats on device sampled from an uniform distribution */
    CURAND_CALL(curandGenerateUniform(cudaRng, devUniformRngData, cudaRngChunkSize));
    /* Copy device memory to host */
    CUDA_CALL(cudaMemcpy(hostUniformRngData, devUniformRngData, cudaRngChunkSize * sizeof(float), cudaMemcpyDeviceToHost));
    // Initialize the RNG index
    nextNormalRngIndex = 0;
    nextUniformRngIndex = 0;
  }
#endif
  // EEW
  char newPath[FILENAME_MAX];
  char seedChar[FILENAME_MAX];
  char tmpChar[FILENAME_MAX];

  // Create the new folder
  sprintf(seedChar, "%d", seed);
  strcpy(newPath, resultsFolder.c_str());
  strcat(newPath, "/");
  strcat(newPath, seedChar);
  strcpy(tmpChar, "mkdir ");
  strcat(tmpChar, newPath);
  if (system(tmpChar) == -1)
  {
    std::cerr << "Error: " << strerror(errno) << std::endl;
    std::cerr << tmpChar << std::endl;
    exit(1);
  }

  // Copy the config file
  strcat(newPath, "/");
  if (getcwd(currentPath, FILENAME_MAX) == NULL)
  {
    std::cerr << "Error: " << strerror(errno) << std::endl;
    std::cerr << currentPath << std::endl;
    exit(1);
  }

  strcpy(tmpChar, "cp ");
  strcat(tmpChar, currentPath);
  strcat(tmpChar, "/");
  strcat(tmpChar, configFileName.c_str());
  strcat(tmpChar, " ");
  strcat(tmpChar, newPath);
  strcat(tmpChar, configFileName.c_str());
  if (system(tmpChar) == -1)
  {
    std::cerr << "Error: " << strerror(errno) << std::endl;
    std::cerr << tmpChar << std::endl;
    exit(1);
  }

  // Change to the new working directory
  strcpy(newPath, resultsFolder.c_str());
  strcat(newPath, "/");
  strcat(newPath, seedChar);
  if (chdir(newPath) == -1)
  {
    std::cerr << "Error: " << strerror(errno) << std::endl;
    std::cerr << newPath << std::endl;
    exit(1);
  }
  initSaveResults(tmpStr.str());
  return true;
}

// Wrapper to allow CUDA and GSL RNGs using the same implementation
double NetDyn::getNormalRng()
{
  double returnValue;

#ifdef CUDA_ENABLED
  if(cuda)
  {
    returnValue = hostNormalRngData[nextNormalRngIndex];
    nextNormalRngIndex++;
    // If we have reached the last available RNG, generate a new chunk
    if(nextNormalRngIndex == cudaRngChunkSize)
    {
      //std::cout << "Reached end of RNG chunk. Creating next one at step: " << step << "\n";
      /* Generate n floats on device sampled from a normal distribution */
      CURAND_CALL(curandGenerateNormal(cudaRng, devNormalRngData, cudaRngChunkSize, 0.0, 1.0));
      /* Copy device memory to host */
      CUDA_CALL(cudaMemcpy(hostNormalRngData, devNormalRngData, cudaRngChunkSize * sizeof(float), cudaMemcpyDeviceToHost));
      nextNormalRngIndex = 0;
    }
  }
  else
#endif
  {
#ifdef OPENMP_ENABLED
    if (parallel)
      returnValue = gsl_ran_gaussian_ziggurat(parallelRng[omp_get_thread_num()], 1.);
    else
#endif
      returnValue = gsl_ran_gaussian_ziggurat(rng, 1.);
  }
  return returnValue;
}

// Wrapper to allow CUDA and GSL RNGs using the same implementation
double NetDyn::getUniformRng()
{
  double returnValue;

#ifdef CUDA_ENABLED
  if(cuda)
  {
    returnValue = hostUniformRngData[nextUniformRngIndex];
    nextUniformRngIndex++;
    // If we have reached the last available RNG, generate a new chunk
    if(nextUniformRngIndex == cudaRngChunkSize)
    {
      //std::cout << "Reached end of RNG chunk. Creating next one at step: " << step << "\n";
      /* Generate n floats on device sampled from a normal distribution */
      CURAND_CALL(curandGenerateUniform(cudaRng, devUniformRngData, cudaRngChunkSize));
      /* Copy device memory to host */
      CUDA_CALL(cudaMemcpy(hostUniformRngData, devUniformRngData, cudaRngChunkSize * sizeof(float), cudaMemcpyDeviceToHost));
      nextUniformRngIndex = 0;
    }
  }
  else
#endif
  {
#ifdef OPENMP_ENABLED
    if (parallel)
      returnValue = gsl_rng_uniform(parallelRng[omp_get_thread_num()]);
    else
#endif
      returnValue = gsl_rng_uniform(rng);
  }
  return returnValue;
}

#ifdef OPENMP_ENABLED
bool NetDyn::seedParallelRng(int seeds)
{
  parallelRng = new gsl_rng* [seeds];
  int newseed;
  std::cout << "Initializing parallel RNGs with seeds: ";
  for (int i = 0; i < seeds; i++)
  {
    parallelRng[i] = gsl_rng_alloc(gsl_rng_taus2);
    // There's a reason for that constant (but I forgot)
    newseed = int(gsl_rng_uniform(rng) * 2147483647);
    std::cout << newseed << " ";
    gsl_rng_set(parallelRng[i], newseed);
  }
  std::cout << "\n";
  return true;
}
#endif

void NetDyn::initSaveResults(std::string tmpStr, std::string filename)
{
  if (filename == "")
    filename = savedFileName;
  std::ofstream savedFile(filename.c_str());
  if (!savedFile.is_open())
  {
    std::cerr << "There was an error opening file " << savedFileName;
    exit(1);
  }

  // Logo - file will be used in MATLAB, so use % for non data lines
  savedFile
      << "%-------------------------------------------------------------\n"
      << "% Neurondyn \n"
      << "% Copyright (c) 2009-2015 Javier G. Orlandi <javierorlandi@javierorlandi.com> \n"
      << "%-------------------------------------------------------------\n";

  savedFile << tmpStr;
  savedFile.close();
}

void NetDyn::saveResults(std::string tmpStr, std::string filename)
{
  if (filename == "")
    filename = savedFileName;
  std::ofstream savedFile(filename.c_str(), std::ios::app);
  if (!savedFile.is_open())
  {
    std::cerr << "There was an error opening the file " << savedFileName;
    exit(1);
  }
  savedFile << tmpStr;
  savedFile.close();
  return;
}

void NetDyn::configureSpikeRecord(int maxsize, std::string filename, double subset)
{
  std::ofstream savedFile(filename.c_str(), std::ios::out);
  if (!savedFile.is_open())
  {
    std::cerr << "There was an error opening the file " << savedFileName;
    exit(1);
  }
  savedFile.close();

  dSpikeRecord = 0;
  dSpikeRecordLimit = maxsize;
  dSpikesFile = filename;
  int pos = filename.rfind(".txt");
  dSpikesSubsetFile = filename;
  dSpikesSubsetFile.replace(pos, 4, "Subset.txt");
  if (subset < 1.)
    dSpikeSubsetSize = floor(nNumber * subset);
  else
    dSpikeSubsetSize = subset;

  savedFile.open(dSpikesSubsetFile.c_str());
  if (!savedFile.is_open())
  {
    std::cerr << "There was an error opening the file " << savedFileName;
    exit(1);
  }
  savedFile.close();
  dSpikeNeuron = new int[dSpikeRecordLimit];
  dSpikeStep = new int[dSpikeRecordLimit];
  for (int i = 0; i < dSpikeRecordLimit; i++)
  {
    dSpikeNeuron[i] = 0;
    dSpikeStep[i] = 0;
  }
}

void NetDyn::recordSpike(int i, int st)
{
  dSpikeNeuron[dSpikeRecord] = i;
  dSpikeStep[dSpikeRecord] = st;
  dSpikeRecord++;
  if (dSpikeRecord == dSpikeRecordLimit)
  {
    processSpikes();
  }
}

void NetDyn::processSpikes(int size)
{
  if (size == -1)
    size = dSpikeRecordLimit;
  std::stringstream tmpStr, tmpStrSubset;
  tmpStr.precision(15);
  tmpStrSubset.precision(15);
  for (int i = 0; i < size; i++)
  {
    tmpStr << dSpikeNeuron[i] << " " << dSpikeStep[i] * dt << "\n";
    if (dSpikeNeuron[i] < dSpikeSubsetSize)
      tmpStrSubset << dSpikeNeuron[i] << " " << dSpikeStep[i] * dt << "\n";
    dSpikeNeuron[i] = 0;
    dSpikeStep[i] = 0;
  }
  dSpikeRecord = 0;
  if (!dryRun)
  {
    saveResults(tmpStr.str(), dSpikesFile);
    saveResults(tmpStrSubset.str(), dSpikesSubsetFile);
  }

  std::cout << "Batch of " << size << " spikes processed. Time:" << step* dt * 1e-3 << " s.\n";
  if (tracing)
    processTraces(traceRecord);
}

void NetDyn::recordTraces()
{
  traceTime[traceRecord] = step * dt;
  for (int i = 0; i < numberTraces; i++)
  {
    traceV[i][traceRecord] = v[tracedNeuron[i]];
    traceU[i][traceRecord] = u[tracedNeuron[i]];
    //       traceI[i][traceRecord] = I[tracedNeuron[i]];
    if (depression)
      traceD[i][traceRecord] = D[tracedNeuron[i]];
    if (AMPA)
      traceI_AMPA[i][traceRecord] = I_AMPA[tracedNeuron[i]];
    if (NMDA)
      traceI_NMDA[i][traceRecord] = I_NMDA[tracedNeuron[i]];
    if (GABA)
      traceI_GABA[i][traceRecord] = I_GABA[tracedNeuron[i]];
    if (WNOISE)
    {
      traceI_WNOISE[i][traceRecord] = I_WNOISE[tracedNeuron[i]] * sqrt(dt);
      traceI[i][traceRecord] = I[tracedNeuron[i]] - I_WNOISE[tracedNeuron[i]];
      //                +I_WNOISE[tracedNeuron[i]]*sqrt(dt);
    }
    else
      traceI[i][traceRecord] = I[tracedNeuron[i]];
  }
  traceRecord++;
  if (traceRecord == STD_DATA_COUNT)
    processTraces();
}

void NetDyn::processTraces(int size)
{
  if (size == -1)
    size = STD_DATA_COUNT;
  std::stringstream tmpStr, completefile;
  tmpStr.precision(15);
  for (int i = 0; i < numberTraces; i++)
  {
    tmpStr.clear();
    tmpStr.str("");
    for (int j = 0; j < size; j++)
    {
      tmpStr << traceTime[j] << " " << traceV[i][j] << " "
             << traceU[i][j] << " " << traceI[i][j] << " ";
      if (depression)
        tmpStr << traceD[i][j] << " ";
      if (AMPA)
        tmpStr << traceI_AMPA[i][j] << " ";
      if (NMDA)
        tmpStr << traceI_NMDA[i][j] << " ";
      if (GABA)
        tmpStr << traceI_GABA[i][j] << " ";
      if (WNOISE)
        tmpStr << traceI_WNOISE[i][j] << " ";
      tmpStr << "\n";
    }
    completefile.clear();
    completefile.str("");
    completefile << traceFile << "_" << i << ".txt";
    if (!dryRun)
      saveResults(tmpStr.str(), completefile.str());
  }
  traceRecord = 0;
}

void NetDyn::configureTraces(int num, std::string filename, double sampTime)
{
  tracing = true;

  if (num == 0)
  {
    tracing = false;
    return;
  }

  int* neuronPool;
  neuronPool = new int[nNumber];
  numberTraces = num;
  tracedNeuron = new int[numberTraces];
  traceFile = filename;

  // Simulation steps between consecutive samples
  traceSamplingTime = floor(sampTime / dt);

  for (int i = 0; i < nNumber; i++)
    neuronPool[i] = i;

  // For now select neurons at random
  gsl_ran_choose(rng, tracedNeuron, numberTraces, neuronPool, nNumber, sizeof(int));

  std::cout << "Selected neurons to trace: ";
  for (int i = 0; i < numberTraces; i++)
    std::cout << tracedNeuron[i] << " ";
  std::cout << "\n";

  // Allocate memory

  traceTime = new double[STD_DATA_COUNT];
  traceV = new double* [numberTraces];
  traceI = new double* [numberTraces];
  traceU = new double* [numberTraces];
  if (depression)
    traceD = new double* [numberTraces];
  if (AMPA)
    traceI_AMPA = new double* [numberTraces];
  if (NMDA)
    traceI_NMDA = new double* [numberTraces];
  if (GABA)
    traceI_GABA = new double* [numberTraces];
  if (WNOISE)
    traceI_WNOISE = new double* [numberTraces];
  for (int i = 0; i < numberTraces; i++)
  {
    traceV[i] = new double[STD_DATA_COUNT];
    traceI[i] = new double[STD_DATA_COUNT];
    traceU[i] = new double[STD_DATA_COUNT];
    if (depression)
      traceD[i] = new double[STD_DATA_COUNT];
    if (AMPA)
      traceI_AMPA[i] = new double[STD_DATA_COUNT];
    if (NMDA)
      traceI_NMDA[i] = new double[STD_DATA_COUNT];
    if (GABA)
      traceI_GABA[i] = new double[STD_DATA_COUNT];
    if (WNOISE)
      traceI_WNOISE[i] = new double[STD_DATA_COUNT];
  }
  traceRecord = 0;

  // Create the file headers
  std::stringstream completefile;
  std::stringstream tmpStr;
  for (int i = 0; i < numberTraces; i++)
  {
    completefile.clear();
    completefile.str("");
    completefile << traceFile << "_" << i << ".txt";
    tmpStr.clear();
    tmpStr.str("");
    tmpStr
        << "%-------------------------------------------------------------\n"
        << "% Neurondyn \n"
        << "% Copyright (c) 2009-2015 Javier G. Orlandi <javierorlandi@javierorlandi.com> \n"
        << "%-------------------------------------------------------------\n"
        << "% Trace for neuron " << tracedNeuron[i] << "."
        << " Neuron Type: c = " << neuronType[tracedNeuron[i]] << "\n"
        << "% Available neurotransmitters:";
    if ((neurotransmitter[neuronType[tracedNeuron[i]]] & NT_AMPA) && AMPA)
      tmpStr << " AMPA";
    if ((neurotransmitter[neuronType[tracedNeuron[i]]] & NT_NMDA) && NMDA)
      tmpStr << " NMDA";
    if ((neurotransmitter[neuronType[tracedNeuron[i]]] & NT_GABA) && GABA)
      tmpStr << " GABA";
    tmpStr << "\n"
           << "%-------------------------------------------------------------\n"
           << "% Traces: \n"
           << "% t | v | u | I";
    if (depression)
      tmpStr << " | D";
    if (AMPA)
      tmpStr << " | I_AMPA";
    if (NMDA)
      tmpStr << " | I_NMDA";
    if (GABA)
      tmpStr << " | I_GABA";
    if (WNOISE)
      tmpStr << " | I_WNOISE";
    tmpStr << "\n"
           << "%-------------------------------------------------------------\n";
    saveResults(tmpStr.str(), completefile.str());
  }
}

bool NetDyn::burstCheck()
{
  bool burstFound = false;
  // If we are outside a burst, check for burst
  if (burstDetectorPossibleBurst == false)
  {
    if (burstDetectorMemory->getCurrentValue() > burstDetectorMemory->getMean() + 2 * burstDetectorMemory->getStdDeviation())
    {
      burstDetectorPossibleBurst = true;
      burstDetectorSpikeCount += burstDetectorMemory->getCurrentValue();
      burstDetectorFirstStepAboveThreshold = step;
      burstDetectorLastStepAboveThreshold = step;
    }
  }
  // Else we are in a burst and check for continuation
  else
  {
    // The burst continues
    if (burstDetectorMemory->getCurrentValue() > burstDetectorMemory->getMean() + 2 * burstDetectorMemory->getStdDeviation())
    {
      burstDetectorSpikeCount += burstDetectorMemory->getCurrentValue();
      burstDetectorLastStepAboveThreshold = step;
    }
    // The burst finished in the previous step
    else
    {
      // Check its size for a valid burst
      if (burstDetectorSpikeCount >= burstDetectorMinimumSpikesPerNeuron * nNumber)
      {
        std::cout << "Burst found! At time " << (double(burstDetectorLastStepAboveThreshold - burstDetectorFirstStepAboveThreshold) / 2. + double(burstDetectorFirstStepAboveThreshold)) * dt * 1e-3 << " s, with a duration of " << double(burstDetectorLastStepAboveThreshold - burstDetectorFirstStepAboveThreshold) * dt * 1e-3 << "s and " << burstDetectorSpikeCount << " spikes.\n";
        burstDetectorStorage.push_back((double(burstDetectorLastStepAboveThreshold - burstDetectorFirstStepAboveThreshold) / 2. + double(burstDetectorFirstStepAboveThreshold)) * dt * 1e-3);
        burstFound = true;
      }
      // The burst was too small
      else
      {
      }
      // Reset burst variables
      burstDetectorSpikeCount = 0;
      burstDetectorPossibleBurst = false;
      burstDetectorFirstStepAboveThreshold = 0;
      burstDetectorLastStepAboveThreshold = 0;
    }
  }
  return burstFound;
}
