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

Now let's run autogen configure:

    ./autogen.sh
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
    file = "networks/periodic5x5r400/cons100_var.txt";
};

positions:
{
    file = "networks/periodic5x5r400/cons100_var_pos.txt";
};

results:
{
    path = "data";
};

# To determine the size and time of nuclei igntion
subnetwork:
{
    active = false;
    files = 81;
    # Reduced matrix
    connectivity:
    {
        file = "../neuron11_data/networks/periodic10x10r50v4/subnetworkRadial.txt";
    };
    # Real number of connections of each neuron in the new matrix
    externalConnections:
    {
        active = false;
        file = "../networks/periodic10x10r50v4/extconRadial.txt";
    };
};

logging:
{
    traces:
    {
        number = 1;             # Number of neurons to trace
        file = "neuronTrace";   # Beginning of the file name
        sampling = 1;           # Take a sample each XX ms
    };
    spikes:
    {
        spikespersave = 200000;    # Number of spikes stored in memory before writing to disk
        file = "SpikeRecord.txt";   # File where the spikes are stored
        subset = 0.1;              # Create another spike record for a subset of neurons
    };
};

simulation:
{
    parallel = true;            # If the system uses OpenMP
    timestep = 0.1;             # Simulation step (in ms)
    totalTime = 2000e3;          # Total simulation time (in ms)
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
                p =  1.;
                neurotransmitters = "AMPA NMDA";
              });
    synapse:
    {
        AMPA:
        {
            active = true;
            g = 170.612;
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
            g = -17;
            tau = 20;
        };
        minis:
        {
            active = true;
	        multiplicative = false;
            depression = false;
            g_AMPA = 18;
            g_GABA = -1;
            tau = 30;
            tau_AMPA = 10;
            tau_GABA = 20;
        };
        depression:
        {
            active = true;
            beta_AMPA = 0.8;
            tau_AMPA = 7500;
            beta_GABA = 0.95;
            tau_GABA = 100;
	        #beta = 0.01;
	        #tau = 20;
        };
        # Activity Dependent Plasticity
        adp:
        {
            active = false;
            AMPA:
            {
                active = false;
                logging_dt = 1000;
                time_window = 10000;
                desired_period = 10000;
                Delta_g = 5;
                D_jump = 0.2;
                store_g = true;
                file = "gAMPA.txt";
            };
            GABA:
            {
                active = false;
                D_minimum = 0.4;
                Delta_g = -5;
                store_g = true;
                file = "gGABA.txt";
            };
        };
    };
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
        adaptation_time = 50e3; # in ms
        target_IBI = 3.38; # in s
        accuracy = 0.05; # in s
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
        active = false;
        delta = -0.1; # Initial g variation. Initial g is defined in the AMPA section
        delta_type = "LOG"; # If the increment is performed in "lin" or "log" scale
        required_bursts = 100; # Number of bursts used to compute the mean IBI
        maximum_IBI = 100; # in s
    };

    # Simulates an external stimulation protocol by applying a bipolar square pulse
    # to the culture every X seconds with varying amplitude and number of repetitions
    external_stimulation:
    {
        active = false;
        pulse_duration = 10; # in ms (this is the total duration of half wave)
        pulse_period = 5e3; # in ms
        pulse_initial_amplitude = 1; # current amplitude A (the pulse goes from +A to -A)
        pulse_delta_amplitude = 1;    # increments in amplitude
        pulse_final_amplitude = 50; #
        pulse_repetitions = 1; # times it repeats a given pulse
        resistance_type = "Gaussian"; # type of distribution to generate the resistance of each neuron to the pulse
        resistance_mean = 10;
        resistance_std = 2;
        inhibit_minis = true;  # sets the mini str to 0 so spontaneous bursting does not interfere with the protocol
        speed_up = true;       # stops after the full activation is reached
        fine_transition = true; # as soon as the system reaches full activation it goes back one delta and continues with a finer grain
        fine_transition_delta_divider = 10;
        fine_transition_repetitions_multiplier = 10;
        fine_activation_threshold = 0.8;
        full_activation_threshold = 0.99;
        file = "stimulation.txt";
        store_spikes = true;
    }
};
```

