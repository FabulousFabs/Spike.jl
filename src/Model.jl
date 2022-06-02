using Parameters;

include("Neurons.jl");
include("Synapses.jl");

@with_kw mutable struct Model
    Neurons::Array{NeuronGroup} = NeuronGroup[];
    Synapses::Array{Synapses} = Synapses[];
end