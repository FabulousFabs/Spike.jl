"""
Provides the main structure for NeuronGroups, standard operations that will be executed
at runtime, and presets of typically used neurons.
"""

using Parameters;

include("Magic.jl");

@with_kw mutable struct NeuronGroup <: SpikeObject
    """
    Main structure for creating neurons. Note that, for convenience, you can often simply use
    one of the auxiliary functions provided in this file (e.g., LIF(N = 10)) rather than re-
    writing expressions yourself (unless that is explicitly required).

    INPUTS:
        N::Int                          -   Number of neurons in group
        eq::Expr                        -   Equation determining the behaviour of neurons in this group.
        thr::Expr                       -   Condition for determining events. (default = :(v > v_th))
        reset::Expr                     -   Reset rule after events. (default = :(v = v_reset))
        method::String                  -   Method used for differentiating. (default = "rk2")
        parameters::Dict{Symbol, Any}   -   All parameters for all neurons. (default = Dict())
    """

    N::Int
    eq::Expr
    thr::Expr = :(v > v_th)
    reset::Expr = :(v = v_reset)
    method::String = "rk2"
    parameters::Dict{Symbol, Any} = Dict()

    __built::Bool = false
    __normeqs::Dict{Symbol, Expr} = Dict()
    __diffeqs::Dict{Symbol, Expr} = Dict()
end

function NeuronGroup_update(; neurons::NeuronGroup, dt::Float64)
    """

    """


end

function LIF(; N::Int = 1, normalised = false)::NeuronGroup
    """
    
    """

    if normalised == false
        return NeuronGroup(N = N, 
                           eq = :(dv_dt = (.-(v .- E_L) .+ R .* I) ./ 𝜏_m;
                                  I_t = .-(v .- E_L)), 
                           parameters = Dict(:v => -70.5 * ones(N), 
                                             :E_L => -70.5 * ones(N), 
                                             :R => 1.0 * ones(N),
                                             :I => 0.0 * ones(N),
                                             :𝜏_m => 10.0 * ones(N)));
    end

    return NeuronGroup()
end
