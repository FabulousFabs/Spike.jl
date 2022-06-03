"""

"""

using Parameters;

include("Neurons.jl");
include("Synapses.jl");

@with_kw mutable struct Model
    Neurons::Vector{NeuronGroup} = NeuronGroup[];
    Synapses::Vector{Synapses} = Synapses[];
end