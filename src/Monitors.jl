"""
Provides the wrapper for monitors that can be added to the network.
"""

using Parameters;

include("Magic.jl");
include("Neurons.jl");
include("Synapses.jl");

"""
    StateMonitor <: SpikeObject

Main structure for monitoring states of neurons or synapses.

    INPUTS:
        obj::Any                                    -   Object to be monitored.
        vars::Vector{Symbol}                        -   Parameters to be monitored in object.
        every::Float64                              -   Delay between calls. (default = 1e-3)
        t::Vector{Float64}                          -   (Internal) Time vector. (default = Float64[])
        states::Dict{Symbol, Array{Float64, 2}}     -   (Internal) State vectors. (default = Dict())
        __built::Bool                               -   (Internal) Has this monitor been built? (default = false)
        __last::Float64                             -   (Internal) Last check. (default = -Inf)
"""
@with_kw mutable struct StateMonitor <: SpikeObject
    obj::Any
    vars::Vector{Symbol}
    every::Float64 = 1e-3

    t::Vector{Float64} = Float64[];
    states::Dict{Symbol, Array{Float64, 2}} = Dict()

    __built::Bool = false
    __last::Float64 = -Inf
end

"""
    step(monitor::StateMonitor; dt::Float64, 
                                t::Float64)::StateMonitor
    
Performs one time step of the monitor. This is an internal function and should not be called manually.

    INPUTS:
        monitor::StateMonitor   -   StateMonitor to perform time step on.
        dt::Float64             -   Time step size.
        t::Float64              -   Current time.

    OUTPUTS:
        monitor::StateMonitor   -   Self
"""
@fastmath function step(monitor::StateMonitor; dt::Float64, t::Float64)::StateMonitor
    # run step if in timing
    if monitor.__last + monitor.every <= t
        push!(monitor.t, t);

        for var::Symbol ∈ monitor.vars
            if typeof(monitor.obj) == NeuronGroup
                monitor.states[var] = cat(monitor.states[var], monitor.obj.parameters[var]; dims = 2);
            elseif typeof(monitor.obj) == Synapses
                monitor.states[var] = cat(monitor.states[var], monitor.obj.__parameters[var]; dims = 2);
            end
        end
    end

    monitor;
end

"""
    EventMonitor <: SpikeObject

Main structure for monitoring events.

    INPUTS:
        obj::NeuronGroup                    -   Object to be monitored.
        event::Symbol                       -   Event to be tracked. (default = :spike)
        t::Vector{Float64}                  -   (Internal) Time vector. (default = Float64[])
        i::Vector{Int}                      -   (Internal) State vectors. (default = Int[])
        __built::Bool                       -   (Internal) Has this monitor been built? (default = false)
"""
@with_kw mutable struct EventMonitor <: SpikeObject
    obj::NeuronGroup
    event::Symbol = :spike

    t::Vector{Float64} = Float64[];
    i::Vector{Int} = Int[];

    __built::Bool = false
end

"""
    step(monitor::EventMonitor; dt::Float64,
                                t::Float64)::EventMonitor

Performs one time step of the monitor. This is an internal fuction and should not be called manually.
"""
@fastmath function step(monitor::EventMonitor; dt::Float64, t::Float64)::EventMonitor
    # collect events and log
    ids::Vector{Int} = collect(1:length(monitor.obj.__eventlog[monitor.event]))[monitor.obj.__eventlog[monitor.event]];
    append!(monitor.t, t .* ones(size(ids, 1)));
    append!(monitor.i, ids);

    monitor;
end

