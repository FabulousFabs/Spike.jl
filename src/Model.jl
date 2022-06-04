"""

"""

using Parameters;

include("Neurons.jl");
include("Synapses.jl");
include("Operations.jl");

@with_kw mutable struct Model
    Neurons::Dict{Symbol, NeuronGroup} = Dict();
    Synapses::Dict{Symbol, Synapses} = Dict();
    Operations::Dict{Symbol, Operation} = Dict();
end
