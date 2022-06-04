module Spike

include("Magic.jl");
include("Model.jl");
include("Neurons.jl");
include("Synapses.jl");
include("Build.jl");
include("Solvers.jl");

export SpikeObject, cast_magic;
export Model;
export NeuronGroup, LIF;
export rk2;

end