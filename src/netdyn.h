/*
 * Copyright (c) 2009-2010 Javier G. Orlandi <orlandi@dherkova.com>
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

#ifndef _NETDYN_H_
#define _NETDYN_H_

#include <libconfig.h++>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include <iostream>


enum neurotransmitterType
{
    NT_AMPA = 0x01,
    NT_NMDA = 0x02,
    NT_GABA = 0x04
};
#define STD_DATA_COUNT 1000000

class NetDyn
{
	public:
        NetDyn();
        void setConnectivity(std::string filename);
        void printConnectivityMap();
        void configureNeuron();
        void configureNeuronTypes(double* cs, double* ds, double* ps, int* ts, int num);
        void setDefaultNeuronTypes();
        void configureSimulation();
        void loadConfigFile(std::string filename, int param = 0);
        void recordSpike(int i, int st);
        void configureSpikeRecord(int maxsize = STD_DATA_COUNT, std::string filename = "SpikeRecord.txt", double subset = 0.05);
        void configureTraces(int num, std::string filename, double sampTime = 1);
        void processSpikes(int size = -1);
        void recordTraces();
        void processTraces(int size = -1);
        bool seedRng();
        bool seedParallelRng(int seeds);
        int simulationStart();
        bool simulationRun();

        void setCUX(std::string filename);
        void setMiniExploration(std::string strengthFile, std::string timeFile);

        inline bool simulationIsRunning()
            {return running;}
        inline bool simulationIsActive()
            {return active;}
        inline int getOriginalRngSeed()
            {return originalRngSeed;}
        int simulationStep();

        void simulationFinish();

        void initialize();
        void setConstants();
        void setExternalConnections(std::string filename);
        void setPositions(std::string filename);
        void saveNetworkStructure(std::string filename = "NetworkStructure.txt");
        void saveResults(std::string tmpStr, std::string filename = "");
        void initSaveResults(std::string tmpStr);

    private:
        std::string savedFileName;
        gsl_rng* rng;			// RNG structure
        gsl_rng** parallelRng;
        libconfig::Config* configFile;

        bool running, active, parallel, multiplicativeMini, depressionMini;
        bool CUXactive, *CUXoverexpressed, depressionReset, miniExplorationActive;
        int CUXmodel;
        std::string CUXfile, miniExplorationStrengthFile, miniExplorationTimeFile;
        double *miniStrength, *miniTime;

        int *connectivityMap;
        int *connectivityNumber;
        int *inputConnectivityNumber, *numberExcitatoryInputs, *numberInhibitoryInputs;
        int *connectivityFirstIndex;
        int nNumber, totalConnections;
        int nType, numTypes;
        // Neuron variables
        double *C, *K, *V_r, *V_t, *V_p, *a, *b, *p;
        double *v, *u, *I, *I_GABA, *I_AMPA, *I_NMDA, *I_WNOISE, *c, *d, *D, *meanD_AMPA, *meanD_GABA;
        double *positionX, *positionY;
        int *neurotransmitter;
        int *ntTypes, *neuronType, neuronTypes;
        double g_AMPA, g_NMDA, g_GABA, g_mAMPA, g_mGABA, g_WNOISE;
        double tau_NMDA, tau_AMPA, tau_GABA, tau_D, beta, tau_MINI, depressionResetInterval;
        int totalSpikes;

        // Simulation parameters
        double dt, t, totalTime;
        int  simulationSteps, step;

        bool GABA, NMDA, AMPA, MINI, WNOISE, depression;

        double exp_AMPA, exp_NMDA, exp_GABA, dt_times_a, dt_over_tau_D, dt_over_C, strength_WNOISE;
        // Post processing
        int *dSpikeNeuron, dSpikeRecord, dSpikeRecordLimit, *dSpikeStep, dSpikeSubsetSize;
        std::string dSpikesFile, traceFile, dSpikesSubsetFile;
        std::string resultsFolder;
        bool tracing;
        int *tracedNeuron, numberTraces, traceSamplingTime, traceRecord;
        double **traceV, **traceU, **traceI, **traceI_AMPA, **traceI_NMDA, **traceI_GABA;
        double **traceI_WNOISE, **traceD, *traceTime, **traceG_AMPA, **traceG_GABA;
        int *parallelThreadUsage;

        char currentPath[FILENAME_MAX];
        int simulationReturnValue;
        int originalRngSeed;
}; 

#endif
    // _NETDYN_H_

