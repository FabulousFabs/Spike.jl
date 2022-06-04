module Spike

include("Magic.jl");
include("Model.jl");
include("Neurons.jl");
include("Synapses.jl");
include("Operations.jl");
include("Build.jl");
include("Solvers.jl");

export SpikeObject, cast_magic;
export Model;
export NeuronGroup, LIF;
export Synapses;
export Operation;
export rk2, euler, heun;

end