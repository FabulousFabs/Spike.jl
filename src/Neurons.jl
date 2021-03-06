"""
Provides the main structure for NeuronGroups, standard operations that will be executed
at runtime, and presets of typically used neurons.
"""

using Parameters;

include("Magic.jl");
include("Solvers.jl");

"""
    NeuronGroup <: SpikeObject

Main structure for creating neurons. Note that, for convenience, you can often simply use one of the auxiliary functions provided in this file (e.g., [`LIF`](@ref)) rather than rewriting expressions yourself (unless that is explicitly required).

    INPUTS:
        N::Int                                          -   Number of neurons in group
        eq::Expr                                        -   Equation determining the behaviour of neurons in this group.
        method::Function                                -   Function to use for solving differential equations. See Spike::Solvers. (default = euler)
        parameters::Dict{Symbol, Any}                   -   All parameters for all neurons. (default = Dict())
        events::Dict{Symbol, Tuple{Expr, Expr}}         -   Event specifications in the form of events = Dict(:name => (:(condition), :(effect;))). (default = Dict())
        __built::Bool                                   -   (Internal) Has this model been built? (default = false)
        __normeqs::Dict{Symbol, Expr}                   -   (Internal) Built equations. (default = Dict())
        __diffeqs::Dict{Symbol, Expr}                   -   (Internal) Built differential equations. (default = Dict())
        __eventeqs::Dict{Symbol, Dict{Symbol, Expr}}    -   (Internal) Built event equations. (defualt = Dict())
        __eventlog::Dict{Symbol, Vector{Bool}}          -   (Internal) Events at runtime. (default = Dict())
"""
@with_kw mutable struct NeuronGroup <: SpikeObject
    N::Int
    eq::Expr
    method::Function = euler
    parameters::Dict{Symbol, Any} = Dict()
    events::Dict{Symbol, Tuple{Expr, Expr}} = Dict()

    __built::Bool = false
    __normeqs::Dict{Symbol, Expr} = Dict()
    __diffeqs::Dict{Symbol, Expr} = Dict()
    __eventeqs::Dict{Symbol, Dict{Symbol, Expr}} = Dict()
    __eventlog::Dict{Symbol, Vector{Bool}} = Dict()
end

"""
    step(neurons::NeuronGroup; dt::Float64,
                               t::Float64)::NeuronGroup

Performs one time step for all the equations specified for a [`NeuronGroup`](@ref). Note that this includes both state updates as well as event checks and logging.

    INPUTS:
        neuron::NeuronGroup     -   NeuronGroup to perform a step on.
        dt::Float64             -   Time step size.
        t::Float64              -   Current time.
    
    OUTPUTS:
        neuron::NeuronGroup     -   Self
"""
@fastmath function step(neurons::NeuronGroup; dt::Float64, t::Float64)::NeuronGroup
    # update general parameters
    neurons.parameters[:N] = neurons.N;
    neurons.parameters[:t] = t;
    neurons.parameters[:dt] = dt;
    
    # solve normal equations
    for eq::Pair{Symbol, Expr} ??? neurons.__normeqs
        neurons.parameters[eq[1]] = eval(interpolate_from_dict(eq[2], neurons.parameters));
    end
    
    # solve differential equations
    for eq::Pair{Symbol, Expr} ??? neurons.__diffeqs
        neurons.parameters[eq[1]] = neurons.method(sym = eq[1], eq = eq[2], par = neurons.parameters, dt = dt, t = t);
    end

    # evaluate any events and their effects
    for event::Pair{Symbol, Dict{Symbol, Expr}} ??? neurons.__eventeqs
        neurons.__eventlog[event[1]] = eval(interpolate_from_dict(neurons.events[event[1]][1], neurons.parameters));

        if any(neurons.__eventlog[event[1]])
            ps::Dict{Symbol, Any} = Dict()

            for p::Pair{Symbol, Any} ??? neurons.parameters
                if size(p[2], 1) > 1
                    ps[p[1]] = p[2][neurons.__eventlog[event[1]]];
                else
                    ps[p[1]] = p[2];
                end
            end

            for eq::Pair{Symbol, Expr} ??? event[2]
                ps[eq[1]] = neurons.parameters[eq[1]][neurons.__eventlog[event[1]]] = eval(interpolate_from_dict(eq[2], ps));
            end
        end
    end

    neurons;
end

"""
    LIF(; N::Int = 1, 
          normalised::Bool = false)::NeuronGroup

Shorthand to create standard Leaky Integrate-and-Fire neurons.

    INPUTS:
        N::Int                  -   Number of neurons (default = 1)
        normalised::Bool        -   Should values be normalised rather than biological (i.e., spanning [0, 1])? (default = false)
    
    OUTPUTS:
        neurons::NeuronGroup    -   LIF-group of neurons
"""
function LIF(; N::Int = 1, normalised::Bool = false)::NeuronGroup
    if normalised == false
        return NeuronGroup(N = N, 
                           eq = :(dv_dt = (.-(v .- E_L) .+ R .* I) ./ ????_m;
                                  dI_dt = -I),
                           method = rk2,
                           events = Dict(:spike => (:(v .> v_th), :(v = v_reset;))),
                           parameters = Dict(:v => -70.5 * ones(N), 
                                             :v_th => -55.0 * ones(N),
                                             :v_reset => -75.5 * ones(N),
                                             :E_L => -70.5 * ones(N), 
                                             :R => 1000.0 * ones(N),
                                             :I => 0.0 * ones(N),
                                             :????_m => 10.0 * ones(N)));
    end

    return NeuronGroup()
end