using Parameters;

include("Magic.jl");

@with_kw mutable struct Synapses <: SpikeObject
    pre::NeuronGroup
    post::NeuronGroup
    eq::Expr

    on_pre::Expr = :()
    on_post::Expr = :()
    cond::Expr = :()
    prob::Float64 = 0.0
    M::Array{Int, 2} = Int[Int[]];
end