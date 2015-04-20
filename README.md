# neurondyn

Simulates neuronal culture dynamics using a modified Izhikevich model.

A previous version of this code was used in the following paper, and its algorithmic details are explained in its supplementary information:
[Noise focusing and the emergence of coherent activity in neuronal cultures](http://www.nature.com/nphys/journal/v9/n9/full/nphys2686.html)  
Nature Physics **9**, 582–590 (2013)

# Dependencies
[libconfig++](http://www.hyperrealm.com/libconfig/)
[gnu/gsl](http://www.gnu.org/software/gsl/)

# Additional dependencies
[OpenMP](http://openmp.org/wp/)
[CUDA](http://www.nvidia.com/object/cuda_home_new.html)

## Step by step installation guide

### Ubuntu-based systems

If you are on a fresh install of ubuntu server you might need some very basic stuff, like make, gcc, g++ and git, so install them:
  
    sudo apt-get install make gcc g++ git

Now let's start with the real dependencies libconfig++ and gsl

    sudo apt-get install libgsl0-dev libconfig++-dev 

If you want to use OpenMP or CUDA you should already know how to install them.

We should be set, now let's clone this repository:

    cd ~
    git clone https://github.com/orlandi/neurondyn

Now let's configure:

    ./configure

OpenMP and CUDA are enabled by default, you can disable them by running:

    ./configure --enable-openmp=no --enable-cuda=no
    
Now let's compile with make:

    make

We are set. If everything went ok you should have the `neurondyn` executable.

## Additional info
OpenMP is used to parallelize the main loop, although it's only recommended for extremely large networks (more than 100k neurons).

CUDA is used instead of gsl to generate most random numbers. As long as you have a modern mid-range NVIDIA card you should try and enable it. It will run much faster.

You should be able to run OpenMP and CUDA simultaneously (although its behavior is untested).

Make sure that the results path you set in the config file exist!

The input connectivity file is a txt file with two columns (I J). Each row identifies an existing connection of the form (I->J). Neuron index starts at 0

The positions file is a txt file with two columns (X Y). Each row indicates the (X,Y) coordiantes of a given neuron (starting at 0). In the current version the positions are not used for anything.

## Usage

The program reads the file config.cfg located in the same folder and
runs the dynamics. Check the comments on that file to know the available options.
Although you can use any kind of network (from its adjacency matrix in sparse format), this program is bets used in tandem with [neurongen](https://github.com/orlandi/neurongen)


The provided config.cfg file looks something like this:
```
# Default configuration file 

version = 1.0;

connectivity:
{
    file = "data/networks/circle5x5r300/cons1_a066.txt";
};

positions:
{
    file = "data/networks/circle5x5r300/map1.txt";
};

results:
{
    path = "data";
};


logging:
{
    traces:
    {
        number = 0;             # Number of neurons to trace
        file = "neuronTrace";   # Beginning of the file name
        sampling = 1;           # Take a sample each XX ms
    };
    spikes:
    {
        spikespersave = 100000;    # Number of spikes stored in memory before writing to disk
        file = "SpikeRecord.txt";   # File where the spikes are stored
        subset = 0.1;              # Create another spike record for a subset of neurons
    };
};

simulation:
{
    openmp = false;             # If the system uses OpenMP
    cuda = true;                # If the system uses cuda to generate the most costly rngs (remember that CUDA needs to be enabled at comile time with the DCUDA_ENABLED flag)
    #cuda_rng_chunk_size = 81920000; # Make it a multiple of 16384 (5000 of it ~ 300MB)
    cuda_rng_chunk_size = 32768000; # Make it a multiple of 16384 (2000 of it ~ 125MB)
    timestep = 0.1;             # Simulation step (in ms)
    totalTime = 50e3;           # Total simulation time (in ms)
    algorithm = "Euler";        # Algorithm
};

neuron:
{
    soma:
    {
        whitenoise:
        {
            active = true;
            strength = 3e2;
        };
    };
    models = ({
                C = 50;
                K = 0.5;
                V_r = -60;
                V_t = -45;
                V_p = 35;
                a = 0.02;
                b = 0.5;
                c = -50;
                d =  50;
                p =  1.0;
                neurotransmitters = "AMPA NMDA";
              },
              {
                C = 100;
                K = 1;
                V_r = -56;
                V_t = -42;
                V_p = 35;
                a = 0.03;
                b = 8;
                c = -53;
                d =  20;
                p =  0.;
                neurotransmitters = "GABA";
              });
    synapse:
    {
        AMPA:
        {
            active = true;
            g = 9.04640144754;
            tau = 10;
        };
        NMDA:
        {
            active = false;
            g = 1.3;
            tau = 150;
        };
        GABA:
        {
            active = false;
            g = -42;
            tau = 20;
        };
       minis:
        {
            active = true;
            multiplicative = false;
            depression = false;
            g_AMPA = 10.0;
            g_GABA = -0;
            tau = 6.0;
        };
        depression:
        {
            active = true;
            beta_AMPA = 0.8;
            tau_AMPA = 3000;
            beta_GABA = 0.95;
            tau_GABA = 100;
            #beta = 0.01;
            #tau = 20;
        };};
};

protocols:
{
    # Automatically detects bursts and stores their mean time
    burst_detector:
    {
        active = true;
        bin_size = 20; # in ms
        bin_number = 50; # number of bins to use for burst detection
        minimum_spikes_per_neuron = 2; # number of bins to use for burst detection
        file = "bursts.txt";
    };

    # Finds the correct synaptic strength for the desired IBI (burst detector
    # is required)
    adaptive_IBI:
    {
        active = false;
        adaptation_time = 25e3; # in ms
        target_IBI = 5.0; # in s
        accuracy = 0.2; # in s
        base_multiplier = 1.1;
        extrapolation_multiplier = 1.0;
        weight_lower_bound = 0.0;
        weight_upper_bound = 1000.0;
        max_iterations = 20;
    };

    # Computes the IBI for a range of synaptic strength while
    # computing enough bursts for good statistics until it reaches
    # the desired maximum IBI or the final time (the final time
    # is overriden by 2*required_bursts*maximum_IBI)
    # Note that this protocol does not store spike data and you need the burst detector
    burst_transition_exploration:
    {
        active = true;
        delta = -0.1; # Initial g variation. Initial g is defined in the AMPA section
        delta_type = "LOG"; # If the increment is performed in "lin" or "log" scale
        required_bursts = 100; # Number of bursts used to compute the mean IBI
        maximum_IBI = 100; # in s
    };
};
```

