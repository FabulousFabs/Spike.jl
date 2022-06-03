module Spike

include("Magic.jl");
include("Model.jl");
include("Neurons.jl");
include("Build.jl");

export SpikeObject, cast_magic;
export Model;
export NeuronGroup, LIF;


end