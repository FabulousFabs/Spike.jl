"""
Provides the main structure for Synapses as well as standard operations to
be executed at runtime.
"""

using Parameters;

include("Magic.jl");
include("Solvers.jl");

"""
    Synapses <: SpikeObject

Main structure for creating synapses between NeuronGroups. Note that both `cond::Expr` and `prob::Float64` may be used complemantarily (or left as is to default to full connections).

    INPUTS:
        pre::NeuronGroup                                -   Presynaptic neuron group.
        post::NeuronGroup                               -   Postsynaptic neuron group.
        cond::Expr                                      -   Conditional expression for building connectivity matrix. (default = :())
        prob::Float64                                   -   Probability of realising a potential entry in connectivity matrix. (default = 1.0)
        eq::Expr                                        -   Equation determining the behaviour of the synapses. (default = :())
        method::Function                                -   Function to use for solving differential equations. See Spike::Solvers. (default = euler)
        parameters::Dict{Symbol, Any}                   -   Parameter pre-specifications; either func where func(N) is valid or ::Number. (default = Dict())
        on_pre::Dict{Symbol, Expr}                      -   Hooks into presynaptic events and their effects (e.g., Dict(:spike => :(post_I = post_I .+ w;))). (default = Dict())
        on_post::Dict{Symbol, Expr}                     -   Hooks into postsynaptic events and their effects (e.g., Dict(:spike => :(pre_f_t_f = t;))). (default = Dict())
        __built::Bool                                   -   (Internal) Have these synapses been built? (default = false)
        __normeqs::Dict{Symbol, Expr}                   -   (Internal) Built equations. (default = Dict())
        __diffeqs::Dict{Symbol, Expr}                   -   (Internal) Built differential equations. (default = Dict())
        __preeqs::Dict{Symbol, Dict{Symbol, Expr}}      -   (Internal) Built presynaptic event equations. (default = Dict())
        __posteqs::Dict{Symbol, Dict{Symbol, Expr}}     -   (Internal) Built postsynaptic event equations. (default = Dict())
        M_pre::Vector{Vector{Int}}                      -   (Internal) Built forwards connectivity matrix. (default = Vector[Int[]])
        M_post::Vector{Vector{Int}}                     -   (Internal) Built backwards connectivity matrix. (default = Vector[Int[]])
        N::Int                                          -   (Internal) Built number of synapses. (default = 0)
        i::Vector{Int}                                  -   (Internal) Built i-th index of every synapse, indicating presynaptic neuron. (default = Int[])
        j::Vector{Int}                                  -   (Internal) Built j-th index of every synapse, indicating postsynaptic neuron. (default = Int[])
"""
@with_kw mutable struct Synapses <: SpikeObject
    pre::NeuronGroup
    post::NeuronGroup
    cond::Expr = :()
    prob::Float64 = 1.0
    eq::Expr = :()
    method::Function = euler
    parameters::Dict{Symbol, Any} = Dict()
    on_pre::Dict{Symbol, Expr} = Dict()
    on_post::Dict{Symbol, Expr} = Dict()
    __built::Bool = false
    __normeqs::Dict{Symbol, Expr} = Dict()
    __diffeqs::Dict{Symbol, Expr} = Dict()
    __preeqs::Dict{Symbol, Dict{Symbol, Expr}} = Dict()
    __posteqs::Dict{Symbol, Dict{Symbol, Expr}} = Dict()
    M_pre::Vector{Vector{Int}} = Vector[Int[]]
    M_post::Vector{Vector{Int}} = Vector[Int[]]
    N::Int = 0
    i::Vector{Int} = Int[]
    j::Vector{Int} = Int[]
end

"""
    step(synapses::Synapses; dt::Float64,
                             t::Float64)::Synapses

Performs one time step for all equations specified for a group of synapses. Note that this also includes state updates as well as event hooks and effects.

    INPUTS:
        synapses::Synapses      -   Synapses to perform step on.
        dt::Float64             -   Time step size.
        t::Float64              -   Current time.
    
    OUTPUTS:
        synapses::Synapses      -   Self
"""
@fastmath function step(synapses::Synapses; dt::Float64, t::Float64)::Synapses
    # update time parameters
    synapses.parameters[:t] .= t;
    synapses.parameters[:dt] .= dt;
    
    # update general parameters
    synapses.parameters[:N] = synapses.N;
    synapses.parameters[:pre_N] = synapses.pre.N;
    synapses.parameters[:post_N] = synapses.post.N;

    # make renamed presynaptic parameters available
    for par_pre::Pair{Symbol, Any} ∈ synapses.pre.parameters
        alias::Symbol = Meta.parse("pre_" * string(par_pre[1]));
        if isa(par_pre[2], Vector)
            synapses.parameters[alias] = par_pre[2][synapses.i];
        elseif isa(par_pre[2], Number)
            synapses.parameters[alias] = par_pre[2] .* ones(synapses.N);
        else
            @assert false "\nSpike::Synapses::step():\nCould not broadcast parameter `" * string(alias) * "` of unsupported type `" * string(typeof(par_pre[2])) * "`. Allowed = [Vector, Number].";
        end
    end

    # make renamed postsynaptic parameters available
    for par_post::Pair{Symbol, Any} ∈ synapses.post.parameters
        alias::Symbol = Meta.parse("post_" * string(par_post[1]));
        if isa(par_post[2], Vector)
            synapses.parameters[alias] = par_post[2][synapses.j];
        elseif isa(par_post[2], Number)
            synapses.parameters[alias] = par_post[2] .* ones(synapses.N);
        else
            @assert false "\nSpike::Synapses::step():\nCould not broadcast parameter `" * string(alias) * "` of unsupported type `" * string(typeof(par_post[2])) * "`. Allowed = [Vector, Number].";
        end
    end

    # solve normal equations
    for eq::Pair{Symbol, Expr} ∈ synapses.__normeqs
        synapses.parameters[eq[1]] = eval(interpolate_from_dict(eq[2], synapses.parameters));
    end

    # solve differential equations
    for eq::Pair{Symbol, Expr} ∈ synapses.__diffeqs
        synapses.parameters[eq[1]] = synapses.method(sym = eq[1], eq = eq[2], par = synapses.parameters, dt = dt, t = t);
    end

    # evaluate pre-synaptic events
    for preeq::Pair{Symbol, Dict{Symbol, Expr}} ∈ synapses.__preeqs
        if !any(synapses.pre.__eventlog[preeq[1]])
            continue;
        end

        presyn_mask::Vector{Bool} = synapses.pre.__eventlog[preeq[1]];
        presyn_indx::Vector{Int} = collect(1:synapses.pre.N)[presyn_mask];

        synapse_indx::Vector{Int} = findall(∈(presyn_indx), synapses.i);
        synapse_mask::Vector{Bool} = falses(synapses.N);
        synapse_mask[synapse_indx] .= true;

        possyn_indx::Vector{Int} = synapses.j[synapse_indx];
        possyn_mask::Vector{Bool} = falses(synapses.post.N);
        possyn_mask[possyn_indx] .= true;

        for eq::Pair{Symbol, Expr} ∈ preeq[2]
            synapse_update::Vector{Float64} = synapse_mask .* (eval(interpolate_from_dict(eq[2], synapses.parameters)) .- synapses.parameters[eq[1]]);
            target_str::String = String(eq[1]);
            internal_sym::Symbol = Symbol();

            if length(target_str) > 5 && target_str[1:5] == "post_"
                internal_sym = Symbol(target_str[6:end]);
                synapses.post.parameters[internal_sym][unique(possyn_indx)] .+= [sum(synapse_update[synapses.j .== post_indx]) for post_indx ∈ unique(possyn_indx)];
                synapses.parameters[eq[1]][synapse_mask] .= synapses.post.parameters[internal_sym][synapses.j[synapse_mask]];
            elseif length(target_str) > 4 && target_str[1:4] == "pre_"
                internal_sym = Symbol(target_str[5:end]);
                synapses.pre.parameters[internal_sym][unique(presyn_indx)] .+= [sum(synapses_update[synapses.j .== post_indx]) for post_indx ∈ unique(possyn_indx)];
                synapses.parameters[eq[1]][synapse_mask] .= synapses.pre.parameters[internal_sym][synapses.i[synapse_mask]];
            else
                internal_sym = eq[1];
                synapses.parameters[eq[1]][synapse_mask] .+= synapse_update;
            end
        end
    end

    # evaluate post-synaptic events
    for posteq::Pair{Symbol, Dict{Symbol, Expr}} ∈ synapses.__posteqs
        if !any(synapses.post.__eventlog[posteq[1]])
            continue;
        end

        possyn_mask::Vector{Bool} = synapses.post.__eventlog[posteq[1]];
        possyn_indx::Vector{Int} = collect(1:synapses.post.N)[possyn_mask];

        synapse_indx::Vector{Int} = findall(∈(possyn_indx), synapses.j);
        synapse_mask::Vector{Bool} = falses(synapses.N);
        synapse_mask[synapse_indx] .= true;

        presyn_indx::Vector{Int} = synapses.i[synapse_indx];
        presyn_mask::Vector{Bool} = falses(synapses.pre.N);
        presyn_mask[presyn_indx] .= true;
        
        for eq::Pair{Symbol, Expr} ∈ posteq[2]
            synapse_update::Vector{Float64} = synapse_mask .* (eval(interpolate_from_dict(eq[2], synapses.parameters)) .- synapses.parameters[eq[1]]);
            target_str::String = String(eq[1]);
            internal_sym::Symbol = Symbol();

            if length(target_str) > 5 && target_str[1:5] == "post_"
                internal_sym = Symbol(target_str[6:end]);
                synapses.post.parameters[internal_sym][unique!(possyn_indx)] .+= [sum(synapse_update[synapses.j .== post_indx]) for post_indx ∈ unique!(possyn_indx)];
                synapses.parameters[eq[1]][synapse_mask] .= synapses.post.parameters[internal_sym][synapses.j[synapse_mask]];
            elseif length(target_str) > 4 && target_str[1:4] == "pre_"
                internal_sym = Symbol(target_str[5:end]);
                synapses.pre.parameters[internal_sym][unique!(presyn_indx)] .+= [sum(synapses_update[synapses.j .== post_indx]) for post_indx ∈ unique!(possyn_indx)];
                synapses.parameters[eq[1]][synapse_mask] .= synapses.pre.parameters[internal_sym][synapses.i[synapse_mask]];
            else
                internal_sym = eq[1];
                synapses.parameters[eq[1]][synapse_mask] .+= synapse_update;
            end
        end
    end


    synapses;
end