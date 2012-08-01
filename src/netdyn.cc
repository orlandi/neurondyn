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

#include "netdyn.h"
#include <sys/time.h>
#include <unistd.h> 
#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <sstream>

#include <cmath>
#include <ctime>
#include <omp.h>

NetDyn::NetDyn()
{
    rng = NULL;

    connectivityMap = NULL;
    connectivityNumber = NULL;
    connectivityFirstIndex = NULL;

    running = false;
    active = false;
    tracing = false;
    parallel = true;
    simulationReturnValue = 0;
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
    inputFile.open (filename.c_str(), std::ifstream::in);
    if (!inputFile.is_open())
    {
        std::cout << "There was an error opening the file " << filename << "\n";
        exit(1);
        return;
    }

    // First loop through the file to get number of synapses and neurons
    while(std::getline(inputFile, line))
    {
        if(line.at(0) == '%')
            continue;
        tmpStr.clear();
        tmpStr.str(line);
        tmpStr >> nFrom >> nTo;
//        std:: cout << line << " " << nFrom << "\n";
        if(nTo > nNeurons)
            nNeurons = nTo;
        nConnections++;
    }
    if(nFrom > nNeurons)
        nNeurons = nFrom;
    nNeurons++;
    std::cout << "Neurons: " << nNeurons << " Connections: " << nConnections << "\n";
    totalConnections = nConnections;
    nNumber = nNeurons;

    connectivityMap = new int[totalConnections];
    connectivityNumber = new int[nNumber];
    for(int i = 0; i < nNumber; i++)
        connectivityNumber[i] = 0;
    connectivityFirstIndex = new int[nNumber];

    inputFile.clear();
    inputFile.seekg (0, std::ios::beg);
    nCurrent = 0;
    nConnections = 0;
    connectionIndex = 0;
    connectivityFirstIndex[0] = 0;

    std::cout << "Assigning output connectivity...\n";
    while(std::getline(inputFile, line))
    {
        if(line.at(0) == '%')
            continue;
        tmpStr.clear();
        tmpStr.str(line);
        tmpStr >> nFrom >> nTo;

        connectivityMap[connectionIndex] = nTo;
        if(nFrom != nCurrent)
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
    for(int i = 0; i < nNumber; i++)
        inputConnectivityNumber[i] = 0;
    for(int i = 0; i < totalConnections; i++)
        inputConnectivityNumber[connectivityMap[i]]++;


}

int NetDyn::simulationStart()
{
    std::cout << "Starting the simulation...\n";

    time_t tstart, tend;
    tstart = time(0);

    while(simulationStep())
    {}
    std::cout << "Simulation finished!\n";
    tend = time(0);
    std::cout << "It took " << difftime(tend, tstart) << " second(s).\n";
    if(parallel)
    {
        std::cout << "Global thread usage: ";
        for(int i = 0; i < omp_get_max_threads(); i++)
            std::cout << parallelThreadUsage[i] << " ";
        std::cout << "\n";
    }
    processSpikes(dSpikeRecord);
    if(tracing)
        processTraces(traceRecord);

    chdir(currentPath);
    return simulationReturnValue;
}

       
void NetDyn::printConnectivityMap()
{
    std::cout << "Printing connectivity map...\n";
    for(int i = 0; i < nNumber; i++)
    {
        std::cout << i << " (" << connectivityNumber[i] << ") --> ";
        for(int j = connectivityFirstIndex[i]; j < connectivityFirstIndex[i]+connectivityNumber[i]; j++)
            std::cout << connectivityMap[j] << " ";
        std::cout << "\n";
    }
}

void NetDyn::setCUX(std::string filename)
{
    std::ifstream inputFile;
    inputFile.open (filename.c_str(), std::ifstream::in);
    if (!inputFile.is_open())
    {
        std::cout << "There was an error opening the file " << filename << "\n";
        return;
    }
    std::string line;
    int nnum, cux;
    std::istringstream tmpStr;

    while(std::getline(inputFile, line))
    {
        if(line.at(0) == '%')
            continue;
        tmpStr.clear();
        tmpStr.str(line);
        tmpStr >> nnum >> cux;
        if(cux == 1)
            CUXoverexpressed[nnum] = true;
        else
            CUXoverexpressed[nnum] = false;
    }
    std::cout << "CUX file loaded.\n";
}

void NetDyn::setMiniExploration(std::string strengthFile, std::string timeFile)
{
    std::ifstream inputFile;
    std::string line;
    int nnum;
    double value;
    std::istringstream tmpStr;

    inputFile.open(strengthFile.c_str(), std::ifstream::in);
    if (!inputFile.is_open())
    {
        std::cout << "There was an error opening the file " << strengthFile << "\n";
        return;
    }

    while(std::getline(inputFile, line))
    {
        if(line.at(0) == '%')
            continue;
        tmpStr.clear();
        tmpStr.str(line);
        tmpStr >> nnum >> value;
        miniStrength[nnum] = value;
    }
    std::cout << "Mini strength file loaded.\n";
    inputFile.close();

    inputFile.open(timeFile.c_str(), std::ifstream::in);
    if (!inputFile.is_open())
    {
        std::cout << "There was an error opening the file " << timeFile << "\n";
        return;
    }

    while(std::getline(inputFile, line))
    {
        if(line.at(0) == '%')
            continue;
        tmpStr.clear();
        tmpStr.str(line);
        tmpStr >> nnum >> value;
        miniTime[nnum] = value;
    }
    std::cout << "Mini time file loaded.\n";
    inputFile.close();
}

void NetDyn::setPositions(std::string filename)
{
    std::ifstream inputFile;
    inputFile.open (filename.c_str(), std::ifstream::in);
    if (!inputFile.is_open())
    {
        std::cout << "There was an error opening the file " << filename << "\n";
        return;
    }
    std::string line, nNum, posX, posY;
    int nnum;
    double posx, posy;
    std::istringstream tmpStr;
    positionX = new double[nNumber];
    positionY = new double[nNumber];

    while(std::getline(inputFile, line))
    {
        tmpStr.clear();
        tmpStr.str(line);
        tmpStr >> nnum >> posx >> posy;
        positionX[nnum] = posx;
        positionY[nnum] = posy;
    }
}

void NetDyn::loadConfigFile(std::string filename, int param)
{ 
    int *nt;
    std::string ntString, tmpStr, lookupStr;
    bool tmpbool = false;
    
    configFile = new libconfig::Config();
    configFile->readFile(filename.c_str());
    configFile->setAutoConvert(true);

    // Start setting the connectivity
    if(!configFile->lookupValue("connectivity.file", tmpStr))
    {
        std::cout << "Main error. Connectivity file not set\n";
        exit(1);
    }
    setConnectivity(tmpStr);

    // Load the positions
    if(!configFile->lookupValue("positions.file", tmpStr))
    {
        std::cout << "Main error. Positions file not set\n";
        exit(1);
    }
    setPositions(tmpStr);
    
    if(!configFile->lookupValue("results.path", tmpStr))
    {
        std::cout << "Warning. There's no path. Using default one\n";
        tmpStr = "data";
    }
    resultsFolder = tmpStr;

    // Configure neuron models
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
    for(int i = 0; i < neuronTypes; i ++)
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
        if(ntString.find("AMPA") != std::string::npos)
            nt[i] = nt[i] | NT_AMPA;    
        if(ntString.find("NMDA") != std::string::npos)
            nt[i] = nt[i] | NT_NMDA;
        if(ntString.find("GABA") != std::string::npos)
            nt[i] = nt[i] | NT_GABA;
//            std::cout << C[i] << " " << c[i] << " " << K[i] << " " << nt[i] << "\n";
        neurotransmitter[i] = nt[i];
    }
        // Configure remaining parameters
        if(!configFile->lookupValue("simulation.timestep", dt))
            std::cout << "Warning! Missing timestep\n";
        if(!configFile->lookupValue("simulation.totalTime", totalTime))
            std::cout << "Warning! Missing totalTime\n";
        if(!configFile->lookupValue("simulation.parallel", parallel))
            std::cout << "Warning! Missing parallel\n";
        
        // Configure soma
        if(!configFile->lookupValue("neuron.soma.whitenoise.active", WNOISE))
            std::cout << "Warning! White noise not defined\n";
        else if((WNOISE))
            configFile->lookupValue("neuron.soma.whitenoise.strength", g_WNOISE);
        
        // Configure synapse
        if(!configFile->lookupValue("neuron.synapse.AMPA.active", AMPA))
            std::cout << "Warning! AMPA not defined\n";
        else if((AMPA))
        {
            if(!configFile->lookupValue("neuron.synapse.AMPA.g", g_AMPA))
                std::cout << "Warning! g_AMPA not defined\n";
            if(!configFile->lookupValue("neuron.synapse.AMPA.tau", tau_AMPA))
                std::cout << "Warning! tau_AMPA not defined\n";
        }

        if(!configFile->lookupValue("neuron.synapse.NMDA.active", NMDA))
            std::cout << "Warning! NMDA not defined\n";
        else if((NMDA))
        {
            if(!configFile->lookupValue("neuron.synapse.NMDA.g", g_NMDA))
                std::cout << "Warning! g_NMDA not defined\n";
            if(!configFile->lookupValue("neuron.synapse.NMDA.tau", tau_NMDA))
                std::cout << "Warning! tau_NMDA not defined\n";
        }
        if(!configFile->lookupValue("neuron.synapse.GABA.active", GABA))
            std::cout << "Warning! GABA not defined\n";
        else if((GABA))
        {
            if(!configFile->lookupValue("neuron.synapse.GABA.g", g_GABA))
                std::cout << "Warning! g_GABA not defined\n";
            if(!configFile->lookupValue("neuron.synapse.GABA.tau", tau_GABA))
                std::cout << "Warning! tau_GABA not defined\n";
        }
        if(!configFile->lookupValue("neuron.synapse.minis.active", MINI))
            std::cout << "Warning! minis not defined\n";
        else if((MINI))
        {
            if(!configFile->lookupValue("neuron.synapse.minis.g_AMPA", g_mAMPA))
                std::cout << "Warning! g_mAMPA not defined\n";
            if(!configFile->lookupValue("neuron.synapse.minis.g_GABA", g_mGABA))
                std::cout << "Warning! g_mGABA not defined\n";
            if(!configFile->lookupValue("neuron.synapse.minis.tau", tau_MINI))
                std::cout << "Warning! tau_MINI not defined\n";
            if(!configFile->lookupValue("neuron.synapse.minis.multiplicative", multiplicativeMini))
                std::cout << "Warning! multiplicativeMini not defined\n";
            if(!configFile->lookupValue("neuron.synapse.minis.depression", depressionMini))
                std::cout << "Warning! depressionMini not defined\n";
        }
        if(!configFile->lookupValue("neuron.synapse.depression.active", depression))
            std::cout << "Warning! depression not defined\n";
        else if((depression))
        {
            if(!configFile->lookupValue("neuron.synapse.depression.beta", beta))
                std::cout << "Warning! beta not defined\n";
            if(!configFile->lookupValue("neuron.synapse.depression.tau", tau_D))
                std::cout << "Warning! tau_D not defined\n";
        }
     
        // Configure CUX
        if(!configFile->lookupValue("neuron.CUX.active", CUXactive))
        {
            std::cout << "Warning! neuron.CUX.active not defined\n";
            CUXactive = false;
        }
        if(CUXactive)
        {
            if(!configFile->lookupValue("neuron.CUX.model", CUXmodel))
                std::cout << "Warning! neuron.CUX.model not defined\n";
            if(!configFile->lookupValue("neuron.CUX.file", CUXfile))
                std::cout << "Warning! neuron.CUX.file not defined\n";
        }

        // Configure the protocols
        // Depression_reset
        if(!configFile->lookupValue("protocols.depression_reset.active", depressionReset))
        {
            std::cout << "Warning! protocols.depression_reset.active not defined\n";
            depressionReset = false;
        }
        if(depressionReset)
        {
            if(!configFile->lookupValue("protocols.depression_reset.interval", depressionResetInterval))
                std::cout << "Warning! protocols.depression_reset.interval not defined\n";
        }
        // mini_exploration
        if(!configFile->lookupValue("protocols.mini_exploration.active", miniExplorationActive))
        {
            std::cout << "Warning! protocols.mini_exploration.active not defined\n";
            miniExplorationActive = false;
        }
        if(miniExplorationActive)
        {
            if(!configFile->lookupValue("protocols.mini_exploration.strength_file", miniExplorationStrengthFile))
                std::cout << "Warning! protocols.mini_exploration.strength_file\n";
            if(!configFile->lookupValue("protocols.mini_exploration.time_file", miniExplorationTimeFile))
                std::cout << "Warning! protocols.mini_exploration.time_file not defined\n";
        }


        // Initialize the rest
        seedRng();
        if(parallel)
        {
            std::cout << "Number of available threads: " << omp_get_max_threads() << "\n";
            seedParallelRng(omp_get_max_threads());
            parallelThreadUsage = new int[omp_get_max_threads()];
            for(int i = 0; i < omp_get_max_threads(); i++)
                parallelThreadUsage[i] = 0;
        }


        // Initialize variables
        initialize();

        setConstants();

        // Configure recordings
        int maxsize = STD_DATA_COUNT;
        int ntraces;
        double subset = 0.05;
        double sampling;
        tmpStr = "spikeRecord.txt";

        if(!configFile->lookupValue("logging.spikes.spikespersave", maxsize))
            std::cout << "Warning! Missing spikes per save\n";
        if(!configFile->lookupValue("logging.spikes.file", tmpStr))
            std::cout << "Warning! Missing spike record file name\n";
        if(!configFile->lookupValue("logging.spikes.subset", subset))
            std::cout << "Warning! Missing subset size\n";
        configureSpikeRecord(maxsize, tmpStr, subset);

        if(!configFile->lookupValue("logging.traces.number", ntraces))
            std::cout << "Warning! Missing number traces\n";
        if(!configFile->lookupValue("logging.traces.file", tmpStr))
            std::cout << "Warning! Missing traces file name\n";
        if(!configFile->lookupValue("logging.traces.sampling", sampling))
            std::cout << "Warning! Missing sampling time\n";

        configureTraces(ntraces, tmpStr, sampling);

    }

    void NetDyn::setExternalConnections(std::string filename)
    {
        std::ifstream inputFile;
        inputFile.open (filename.c_str(), std::ifstream::in);
        if (!inputFile.is_open())
        {
            std::cout << "There was an error opening the file " << filename << "\n";
            return;
        }
        std::string line;
        int nnum, econs;
        std::istringstream tmpStr;

        while(std::getline(inputFile, line))
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
            neuronType = new int[nNumber];

            if(AMPA)
            {
                I_AMPA = new double[nNumber];
            }
            if(NMDA)
            {
                I_NMDA = new double[nNumber];
            }
            if(GABA)
            {
                I_GABA = new double[nNumber];
            }
            if(WNOISE)
            {
                I_WNOISE = new double[nNumber];
            }
            if(depression)
            {
                D = new double[nNumber];
                meanD_AMPA = new double[nNumber];
                meanD_GABA = new double[nNumber];
            }
            if(CUXactive)
            {
                CUXoverexpressed = new bool[nNumber];
                setCUX(CUXfile);
            }
            if(miniExplorationActive)
            {
                miniStrength = new double[nNumber];
                miniTime = new double[nNumber];
                setMiniExploration(miniExplorationStrengthFile,miniExplorationTimeFile);
            }

            // Preallocate to define the neuron types
            gsl_ran_discrete_t *randist;
            randist =  gsl_ran_discrete_preproc (neuronTypes, p);

            for(int i = 0; i < nNumber; i++)
            {
                u[i] = 0.;
                I[i] = 0.;

                // Assign the neuron type from the type distribution 
                neuronType[i] = gsl_ran_discrete (rng, randist);

                // Override neuronal type if CUX is overexpressed
                if(CUXactive)
                    if(CUXoverexpressed[i] == true)
                        neuronType[i] = CUXmodel;


                v[i] = V_r[neuronType[i]];


                if(AMPA)
                    I_AMPA[i] = 0.;
                if(GABA)
                    I_GABA[i] = 0.;
                if(NMDA)
                    I_NMDA[i] = 0.;
                if(WNOISE)
                    I_WNOISE[i] = 0.;
                D[i] = 1.;
            }

            numberExcitatoryInputs = new int[nNumber];
            numberInhibitoryInputs = new int[nNumber];
            for(int i = 0; i < nNumber; i++)
            {
                numberExcitatoryInputs[i] = 0;
                numberInhibitoryInputs[i] = 0;
            }
            int k, neighIndex, ntype;

            for(int i = 0; i < nNumber; i++)
            {
                ntype = neuronType[i];
                k = connectivityFirstIndex[i];
            //    std::cout << connectivityNumber[i] << "\n";
                for(int j = 0; j < connectivityNumber[i]; j++)
                {
                    neighIndex = connectivityMap[k+j];

                    if((neurotransmitter[ntype] & NT_AMPA) || (neurotransmitter[ntype] & NT_NMDA))
                        numberExcitatoryInputs[neighIndex]++;
                    if((neurotransmitter[ntype] & NT_GABA))
                        numberInhibitoryInputs[neighIndex]++;
                }
            }
            saveNetworkStructure();
            simulationSteps = ceil(totalTime/dt);
            step = 0;
            totalSpikes = 0;
    }

    void NetDyn::saveNetworkStructure(std::string filename)
    {
        std::stringstream tmpStr;
        int nt = 0;

        tmpStr
            << "%-------------------------------------------------------------\n"
            << "% Neuron11 \n"
            << "% Copyright (c) 2009-10 Javier G. Orlandi <dherkova@gmail.com> \n"
            << "%-------------------------------------------------------------\n"
            << "% Neuron idx | Type (0 Ex, 1 In) | k_i | k_i(E) | k_i(I) | k_o\n";

        for(int i = 0; i < nNumber; i++)
        {
            if(neurotransmitter[neuronType[i]] & NT_GABA)
                nt = 1;
            else
                nt = 0;
            tmpStr << i << " " << nt << " "
                   << inputConnectivityNumber[i] << " "
                   << numberExcitatoryInputs[i] << " "
                   << numberInhibitoryInputs[i] << " "
                   << connectivityNumber[i] << "\n";
            if(i % 1000 == 0)
            {
                saveResults(tmpStr.str(), filename);
                tmpStr.clear();
                tmpStr.str("");
            }
        }
        saveResults(tmpStr.str(), filename);
    }

    int NetDyn::simulationStep()
    {
        int i, j, k, neighIndex;
        int ntype;
        double vtmp, tmpMiniStrength, tmpMiniConstant;

        if(tracing)
            if(step % traceSamplingTime == 0)
                recordTraces();

        // Reset the mean depression in the terminals
        if(depression && depressionMini)
        {
            for(i = 0; i < nNumber; i++)
            {
                meanD_AMPA[i] = 0;
                meanD_GABA[i] = 0;
            }
        }

        // Check for spikes
        for(i = 0; i < nNumber; i++)
        {
            ntype = neuronType[i];
            // If there's a spike, reset and propagate
            if(v[i] > V_p[ntype])
            {
                // Soma reset
                v[i] = c[ntype];
                u[i] += d[ntype];
                // Propagate currents to the neighbors;
                k = connectivityFirstIndex[i];
                for(j = 0; j < connectivityNumber[i]; j++)
                {
                    neighIndex = connectivityMap[k+j];

                    if((neurotransmitter[ntype] & NT_AMPA) && AMPA)
                    {
                        I_AMPA[neighIndex] += g_AMPA*D[i];
                    }
                    if((neurotransmitter[ntype] & NT_NMDA) && NMDA)
                        I_NMDA[neighIndex] += g_NMDA*D[i];
                    if((neurotransmitter[ntype] & NT_GABA) && GABA)
                    {
                        I_GABA[neighIndex] += g_GABA*D[i];
                    }
                }
                // Neurotransmitter depletion
                if(depression)
                    D[i] *= beta;
    //            std::cout << D[i]<< "\n";
                // Add the spike to the list
                recordSpike(i, step);
                totalSpikes++;
            }
            // Update the mean depression: assign to each neuron a mean depression of its input connections
            if(depression && depressionMini)
            {
                k = connectivityFirstIndex[i];
                for(j = 0; j < connectivityNumber[i]; j++)
                {
                    neighIndex = connectivityMap[k+j];

                    if((neurotransmitter[ntype] & NT_AMPA) && AMPA)
                        meanD_AMPA[neighIndex] += D[i];
                    if((neurotransmitter[ntype] & NT_GABA) && GABA)
                        meanD_GABA[neighIndex] += D[i];
                }
     
            } 
        }
        // Correctly assign the mean input depression
        if(depression && depressionMini)
        {
            for(i = 0; i < nNumber; i++)
            {
                meanD_AMPA[i] = meanD_AMPA[i]/numberExcitatoryInputs[i];
                meanD_GABA[i] = meanD_GABA[i]/numberInhibitoryInputs[i];
            }
        }

        // Now we sum the contribution of the different currents
        // if a spike happened in the previous loop, it means that it happened in between
        // the previous time step, so currents should be updated accordingly. And they are!
        // Also, we can update them in the same step
        #pragma omp parallel for private(vtmp, ntype, tmpMiniStrength, tmpMiniConstant)
        for(i = 0; i < nNumber; i++)
        {
            ntype = neuronType[i];
            if(parallel)
                parallelThreadUsage[omp_get_thread_num()]++;
            I[i] = 0;

            // Now add the noise
            if(WNOISE)
            {
                if(parallel)
                    I_WNOISE[i] = strength_WNOISE*
                        gsl_ran_gaussian_ziggurat(parallelRng[omp_get_thread_num()], 1.);
                else
                    I_WNOISE[i] = strength_WNOISE*gsl_ran_gaussian_ziggurat(rng, 1.);
            }
            // MINIS HERE - only AMPA and GABA
            if(MINI)
            {
                if(AMPA)
                {
                    if(multiplicativeMini)
                        tmpMiniConstant = double(numberExcitatoryInputs[i])*dt/tau_MINI;
                    else
                        tmpMiniConstant = dt/tau_MINI;

                    if(depressionMini)
                        tmpMiniStrength = meanD_AMPA[i]*g_mAMPA;
                    else
                        tmpMiniStrength = g_mAMPA;

                    // The protocol, overrides everything done before
                    if(miniExplorationActive)
                    {
                        tmpMiniStrength = miniStrength[i];
                        tmpMiniConstant = dt/miniTime[i];
                    }

                    if(parallel)
                    {
                        if(gsl_rng_uniform(parallelRng[omp_get_thread_num()]) < tmpMiniConstant)
                            I_AMPA[i] += tmpMiniStrength;
                    }
                    else
                    {
                        if(gsl_rng_uniform(rng) < tmpMiniConstant)
                            I_AMPA[i] += tmpMiniStrength;
                    }
                }
                if(GABA)
                {
                    if(multiplicativeMini)
                        tmpMiniConstant = double(numberInhibitoryInputs[i])*dt/tau_MINI;
                    else
                        tmpMiniConstant = dt/tau_MINI;

                    if(depressionMini)
                        tmpMiniStrength = meanD_GABA[i]*g_mGABA;
                    else
                        tmpMiniStrength = g_mGABA;

                    if(parallel)
                    {
                        if(gsl_rng_uniform(parallelRng[omp_get_thread_num()]) < tmpMiniConstant)
                            I_GABA[i] += tmpMiniStrength;
                    }
                    else
                    {
                        if(gsl_rng_uniform(rng) < tmpMiniConstant)
                            I_GABA[i] += tmpMiniStrength;
                    }

                }
            }

        // Start to add all the currents
        if(WNOISE)
            I[i] += I_WNOISE[i];
        if(AMPA)
        {
            I[i] += I_AMPA[i];
            I_AMPA[i] *= exp_AMPA;
        }
        if(NMDA)
        {
            I[i] += I_NMDA[i];
            I_NMDA[i] *= exp_NMDA;
        }
        if(GABA)
        {
            I[i] += I_GABA[i];
            I_GABA[i] *= exp_GABA;
        }
        // Now it is time to update the soma - Euler for now
        vtmp = v[i]+dt/C[ntype]*(K[ntype]*(v[i]-V_r[ntype])*(v[i]-V_t[ntype])-u[i]+I[i]);
        u[i] += dt*a[ntype]*(b[ntype]*(v[i]-V_r[ntype])-u[i]);
        v[i] = vtmp;

        if(depression)
            D[i] += dt_over_tau_D*(1.-D[i]);
    }

        // DepressionReset protocol
        if(depressionReset)
        {
            if(step % int(depressionResetInterval/dt) == 0)
            {
                #pragma omp parallel
                for(i = 0; i < nNumber; i++)
                    D[i] = 1.;
            }
        }

    step++;
    if(step % int(10000/dt) == 0)
        processSpikes(dSpikeRecord);
    if(step >= simulationSteps)
        return 0;
    else
        return 1;
}
 
void NetDyn::setConstants()
{
    if(AMPA)
        exp_AMPA = exp(-dt/tau_AMPA);
    if(NMDA)
        exp_NMDA = exp(-dt/tau_NMDA);
    if(GABA)
        exp_GABA = exp(-dt/tau_GABA);
    if(depression)
        dt_over_tau_D = dt/tau_D;
    if(WNOISE)
        strength_WNOISE = sqrt(2.*g_WNOISE/dt);
}

bool NetDyn::seedRng()
{
    std::stringstream tmpStr, seedStr;
	struct timeval tv;
	gettimeofday(&tv,NULL);
    struct tm *tm = localtime(&tv.tv_sec);
	int seed = abs(int(tv.tv_usec/10+tv.tv_sec*100000));	// Creates the seed based on actual time
    originalRngSeed = seed;
	rng = gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set(rng,seed);			// Seeds the previously created RNG
    seedStr << seed << ".txt";
    savedFileName = seedStr.str();

    tmpStr << "% Date: " << tm->tm_mday << "/" << tm->tm_mon +1 << "/" << tm->tm_year + 1900 << ", "
           << "Time: " << tm->tm_hour << ":" << tm->tm_min << ":" << tm->tm_sec << ", "
           << "RNG Seed: " << seed << "\n";
    std::cout << tmpStr.str();

    // EEW
    char newPath[FILENAME_MAX];
    char seedChar [FILENAME_MAX];
    char tmpChar [FILENAME_MAX];

/*    sprintf(seedChar,"%d",seed);
    getcwd(currentPath, FILENAME_MAX);
    strcpy(tmpChar,"data");
    strcat(tmpChar,"/");
    strcat(tmpChar,seedChar);

    mkdir("data",0777);
    
    if(mkdir(tmpChar,0777)==-1) //creating a directory
    {
        std::cerr << "Warning:  " << strerror(errno) << std::endl;
    }
    strcat(currentPath,"/");
    strcat(currentPath, tmpChar);
    chdir(currentPath);
    system("cp ../../config.cfg .");*/

    // Create the new folder
    sprintf(seedChar,"%d",seed);
    strcpy(newPath,resultsFolder.c_str());
    strcat(newPath,"/");
    strcat(newPath,seedChar);
    strcpy(tmpChar, "mkdir ");
    strcat(tmpChar, newPath);
    system(tmpChar);

    // Copy the config file
    strcat(newPath,"/");
    getcwd(currentPath, FILENAME_MAX);
    strcpy(tmpChar, "cp ");
    strcat(tmpChar, currentPath);
    strcat(tmpChar, "/config.cfg ");
    strcat(tmpChar, newPath);
    strcat(tmpChar, "config.cfg");
    system(tmpChar);

    // Change to the new working directory
    strcpy(newPath,resultsFolder.c_str());
    strcat(newPath,"/");
    strcat(newPath,seedChar);
    if(chdir(newPath) == -1)
        std::cerr << "warning: " << strerror(errno) << std::endl;
    initSaveResults(tmpStr.str());
    return true;
}

bool NetDyn::seedParallelRng(int seeds)
{
    parallelRng = new gsl_rng*[seeds];
    int newseed;
    std::cout << "Initializing parallel RNGs with seeds: ";
    for(int i = 0; i < seeds; i++)
    {
        parallelRng[i] = gsl_rng_alloc(gsl_rng_taus2);
        newseed = int(gsl_rng_uniform(rng)*2147483647);
        std::cout << newseed << " ";
    	gsl_rng_set(parallelRng[i], newseed);
    }
    std::cout << "\n";
    return true;
}

void NetDyn::initSaveResults(std::string tmpStr)
{
    std::ofstream savedFile(savedFileName.c_str());
    if (!savedFile.is_open())
    {
        std::cout << "There was an error opening file " << savedFileName;
        return;
    }

    // Logo - file will be used in MATLAB, so use % for non data lines
    savedFile
        << "%-------------------------------------------------------------\n"
        << "% Neuron11 \n"
        << "% Copyright (c) 2009-10 Javier G. Orlandi <dherkova@gmail.com> \n"
        << "%-------------------------------------------------------------\n";

    savedFile << tmpStr;
//    std::cout << tmpStr;
//    savedFile << "\n";
    savedFile.close();
}

void NetDyn::saveResults(std::string tmpStr, std::string filename)
{
    if(filename == "")
        filename = savedFileName;
	std::ofstream savedFile(filename.c_str(), std::ios::app);
    if (!savedFile.is_open())
    {
        std::cout << "There was an error opening the file " << savedFileName;
        return;
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
        std::cout << "There was an error opening the file " << savedFileName;
        return;
    }
    savedFile.close();

    dSpikeRecord = 0;
    dSpikeRecordLimit = maxsize;
    dSpikesFile = filename;
    int pos = filename.rfind(".txt");
    dSpikesSubsetFile = filename;
    dSpikesSubsetFile.replace (pos, 4, "Subset.txt");
    if(subset < 1.)
        dSpikeSubsetSize = floor(nNumber*subset);
    else
        dSpikeSubsetSize = subset;

	savedFile.open(dSpikesSubsetFile.c_str());
    if (!savedFile.is_open())
    {
        std::cout << "There was an error opening the file " << savedFileName;
        return;
    }
    savedFile.close();
    dSpikeNeuron = new int[dSpikeRecordLimit];
    dSpikeStep = new int[dSpikeRecordLimit];
    for(int i = 0; i < dSpikeRecordLimit; i++)
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
    if(dSpikeRecord == dSpikeRecordLimit)
    {
        processSpikes();
    }
}

void NetDyn::processSpikes(int size)
{
    if(size == -1)
        size = dSpikeRecordLimit;
    std::stringstream tmpStr, tmpStrSubset;
    tmpStr.precision(15);
    tmpStrSubset.precision(15);
    for(int i = 0; i < size; i ++)
    {
        tmpStr << dSpikeNeuron[i] << " " << dSpikeStep[i]*dt << "\n";
        if(dSpikeNeuron[i] < dSpikeSubsetSize)
            tmpStrSubset << dSpikeNeuron[i] << " " << dSpikeStep[i]*dt << "\n";
        dSpikeNeuron[i] = 0;
        dSpikeStep[i] = 0;
    }
    dSpikeRecord = 0;
    saveResults(tmpStr.str(), dSpikesFile);
    saveResults(tmpStrSubset.str(), dSpikesSubsetFile);
    
    std::cout << "Batch of " << size << " spikes processed. Time:" << step*dt*1e-3 << " s.\n";
    if(tracing)
        processTraces(traceRecord);
}

void NetDyn::recordTraces()
{
    traceTime[traceRecord] = step*dt;
    for(int i = 0; i < numberTraces; i++)
    {
       traceV[i][traceRecord] = v[tracedNeuron[i]];
       traceU[i][traceRecord] = u[tracedNeuron[i]];
//       traceI[i][traceRecord] = I[tracedNeuron[i]];
       if(depression)
            traceD[i][traceRecord] = D[tracedNeuron[i]];
       if(AMPA)
            traceI_AMPA[i][traceRecord] = I_AMPA[tracedNeuron[i]];
       if(NMDA)
            traceI_NMDA[i][traceRecord] = I_NMDA[tracedNeuron[i]];
       if(GABA)
            traceI_GABA[i][traceRecord] = I_GABA[tracedNeuron[i]];
       if(WNOISE)
       {
            traceI_WNOISE[i][traceRecord] = I_WNOISE[tracedNeuron[i]]*sqrt(dt);
            traceI[i][traceRecord] = I[tracedNeuron[i]]-I_WNOISE[tracedNeuron[i]];
//                +I_WNOISE[tracedNeuron[i]]*sqrt(dt);
       }
       else
            traceI[i][traceRecord] = I[tracedNeuron[i]];
    }
    traceRecord++;
    if(traceRecord == STD_DATA_COUNT)
        processTraces();
}

void NetDyn::processTraces(int size)
{
    if(size == -1)
        size = STD_DATA_COUNT;
    std::stringstream tmpStr, completefile;
    tmpStr.precision(15);
    for(int i = 0; i < numberTraces; i++)
    {
        tmpStr.clear();
        tmpStr.str("");
        for(int j = 0; j < size; j ++)
        {
            tmpStr << traceTime[j] << " " << traceV[i][j] << " " 
                   << traceU[i][j] << " " << traceI[i][j] << " ";
            if(depression)
                tmpStr << traceD[i][j] << " "; 
            if(AMPA)
                tmpStr << traceI_AMPA[i][j] << " ";
            if(NMDA)
                tmpStr << traceI_NMDA[i][j] << " ";
            if(GABA)
                tmpStr << traceI_GABA[i][j] << " ";
            if(WNOISE)
                tmpStr << traceI_WNOISE[i][j] << " ";
            tmpStr << "\n";
        }
        completefile.clear();
        completefile.str("");
        completefile << traceFile << "_" << i << ".txt"; 
        saveResults(tmpStr.str(), completefile.str());
    }
    traceRecord = 0;
}
        
void NetDyn::configureTraces(int num, std::string filename, double sampTime)
{
    tracing = true;

    if(num == 0)
    {
        tracing = false;
        return;
    }

    int *neuronPool;
    neuronPool = new int[nNumber];
    numberTraces = num;
    tracedNeuron = new int[numberTraces];
    traceFile = filename;

    // Simulation steps between consecutive samples
    traceSamplingTime = floor(sampTime/dt);

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
    traceV = new double*[numberTraces];
    traceI = new double*[numberTraces];
    traceU = new double*[numberTraces];
    if(depression)
        traceD = new double*[numberTraces];
    if(AMPA)
        traceI_AMPA = new double*[numberTraces];
    if(NMDA)
        traceI_NMDA = new double*[numberTraces];
    if(GABA)
        traceI_GABA = new double*[numberTraces];
    if(WNOISE)
        traceI_WNOISE = new double*[numberTraces];
    for(int i = 0; i < numberTraces; i++)
    {
        traceV[i] = new double[STD_DATA_COUNT];
        traceI[i] = new double[STD_DATA_COUNT];
        traceU[i] = new double[STD_DATA_COUNT];
        if(depression)
            traceD[i] = new double[STD_DATA_COUNT];
        if(AMPA)
            traceI_AMPA[i] = new double[STD_DATA_COUNT];
        if(NMDA)
            traceI_NMDA[i] = new double[STD_DATA_COUNT];
        if(GABA)
            traceI_GABA[i] = new double[STD_DATA_COUNT];
        if(WNOISE)
            traceI_WNOISE[i] = new double[STD_DATA_COUNT];
    }
    traceRecord = 0;

    // Create the file headers
    std::stringstream completefile;
    std::stringstream tmpStr;
    for(int i = 0; i < numberTraces; i++)
    {
        completefile.clear();
        completefile.str("");
        completefile << traceFile << "_" << i << ".txt";
        tmpStr.clear();
        tmpStr.str("");
        tmpStr
            << "%-------------------------------------------------------------\n"
            << "% Neuron11 \n"
            << "% Copyright (c) 2009-10 Javier G. Orlandi <dherkova@gmail.com> \n"
            << "%-------------------------------------------------------------\n"
            << "% Trace for neuron " << tracedNeuron[i] << "."
            << " Neuron Type: c = " << neuronType[tracedNeuron[i]] << "\n"
            << "% Available neurotransmitters:";
        if((neurotransmitter[neuronType[tracedNeuron[i]]] & NT_AMPA) && AMPA)
            tmpStr << " AMPA";
        if((neurotransmitter[neuronType[tracedNeuron[i]]] & NT_NMDA) && NMDA)
            tmpStr << " NMDA";
        if((neurotransmitter[neuronType[tracedNeuron[i]]] & NT_GABA) && GABA)
            tmpStr << " GABA";
        tmpStr << "\n"
            << "%-------------------------------------------------------------\n"
            << "% Traces: \n"
            << "% t | v | u | I";
        if(depression)
            tmpStr << " | D";
        if(AMPA)
            tmpStr << " | I_AMPA";
        if(NMDA)
            tmpStr << " | I_NMDA";
        if(GABA)
            tmpStr << " | I_GABA";
        if(WNOISE)
            tmpStr << " | I_WNOISE";
        tmpStr << "\n"
            << "%-------------------------------------------------------------\n";
        saveResults(tmpStr.str(), completefile.str());
    }
}

