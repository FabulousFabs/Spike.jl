using Parameters;

include("Neurons.jl");

@with_kw mutable struct Model
    Neurons::Array{Neuron}
end