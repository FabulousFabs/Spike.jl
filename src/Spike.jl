"""
Spike.jl module. Contains all functionality and manages exports.
"""

module Spike

include("Magic.jl");
include("Model.jl");
include("Neurons.jl");
include("Synapses.jl");
include("Operations.jl");
include("Monitors.jl");
include("Build.jl");
include("Solvers.jl");

export SpikeObject, cast_magic;
export Model, run;
export NeuronGroup, LIF;
export Synapses;
export Operation;
export StateMonitor, EventMonitor;
export rk2, euler, heun;

end