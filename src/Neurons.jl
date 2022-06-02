using Parameters;

include("Magic.jl");

@with_kw mutable struct NeuronGroup <: SpikeObject
    N::Int
    eq::Expr
    method::String
end