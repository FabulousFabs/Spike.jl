"""
Provides the wrapper for monitors that can be added to the network.
"""

using Parameters;

include("Magic.jl");
include("Neurons.jl");
include("Synapses.jl");

@with_kw mutable struct StateMonitor <: SpikeObject
    """
    Main structure for creating a StateMonitor.

    INPUTS:
        obj::Any                                    -   Object to be monitored.
        vars::Vector{Symbol}                        -   Parameters to be monitored in object.
        every::Float64                              -   Delay between calls. (default = 1e-3)
        t::Vector{Float64}                          -   (Internal) Time vector. (default = Float64[])
        states::Dict{Symbol, Array{Float64, 2}}     -   (Internal) State vectors. (default = Dict())
        __built::Bool                               -   (Internal) Has this monitor been built? (default = false)
        __last::Float64                             -   (Internal) Last check. (default = -Inf)
    """

    obj::Any
    vars::Vector{Symbol}
    every::Float64 = 1e-3

    t::Vector{Float64} = Float64[];
    states::Dict{Symbol, Array{Float64, 2}} = Dict()

    __built::Bool = false
    __last::Float64 = -Inf
end

@fastmath function step(monitor::StateMonitor; dt::Float64, t::Float64)::StateMonitor
    """
    Perform a time step of the monitor.

    INPUTS:
        monitor::StateMonitor   -   StateMonitor to perform time step on.
        dt::Float64             -   Time step size.
        t::Float64              -   Current time.
    
    OUTPUTS:
        monitor::StateMonitor   -   Self
    """

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

@with_kw mutable struct EventMonitor <: SpikeObject
    """
    Main structure for creating an EventMonitor.

    INPUTS:
        obj::NeuronGroup                    -   Object to be monitored.
        event::Symbol                       -   Event to be tracked. (default = :spike)
        t::Vector{Float64}                  -   (Internal) Time vector. (default = Float64[])
        i::Vector{Int}                      -   (Internal) State vectors. (default = Int[])
        __built::Bool                       -   (Internal) Has this monitor been built? (default = false)
    """

    obj::NeuronGroup
    event::Symbol = :spike

    t::Vector{Float64} = Float64[];
    i::Vector{Int} = Int[];

    __built::Bool = false
end

@fastmath function step(monitor::EventMonitor; dt::Float64, t::Float64)::EventMonitor
    """
    Perform a time step of the monitor.

    INPUTS:
        monitor::EventMonitor   -   EventMonitor to perform time step on.
        dt::Float64             -   Time step size.
        t::Float64              -   Current time.
    
    OUTPUTS:
        monitor::EventMonitor   -   Self
    """

    # collect events and log
    ids::Vector{Int} = collect(1:length(monitor.obj.__eventlog[monitor.event]))[monitor.obj.__eventlog[monitor.event]];
    append!(monitor.t, t .* ones(size(ids, 1)));
    append!(monitor.i, ids);

    monitor;
end

