"""

"""

using Parameters;

include("Magic.jl");
include("Solvers.jl");

@with_kw mutable struct Synapses <: SpikeObject
    """
    
    """

    pre::NeuronGroup
    post::NeuronGroup
    cond::Expr = :()
    prob::Float64 = 0.0
    eq::Expr = :()
    method::Function = rk2
    parameters::Dict{Symbol, Any} = Dict()
    on_pre::Dict{Symbol, Expr} = Dict()
    on_post::Dict{Symbol, Expr} = Dict()

    __parameters::Dict{Symbol, Any} = Dict()
    __normeqs::Dict{Symbol, Expr} = Dict()
    __diffeqs::Dict{Symbol, Expr} = Dict()
    __preeqs::Dict{Symbol, Dict{Symbol, Expr}} = Dict()
    __posteqs::Dict{Symbol, Dict{Symbol, Expr}} = Dict()
    __M_pre::Vector{Vector{Int}} = Vector[Int[]];
    __M_post::Vector{Vector{Int}} = Vector[Int[]];
    __N::Int = 0;
    __i::Vector{Int} = Int[];
    __j::Vector{Int} = Int[];
end

macro add_prefix_to(prefix, symbol)
    return :(eval(Meta.parse(":" * $prefix * string($symbol))))
end

function step(synapses::Synapses; dt::Float64, t::Float64)
    """

    """

    # update general parameters
    synapses.__parameters[:N] = synapses.__N;
    synapses.__parameters[:pre_N] = synapses.pre.N;
    synapses.__parameters[:post_N] = synapses.post.N;
    synapses.__parameters[:t] = t;
    synapses.__parameters[:dt] = dt;

    # make renamed presynaptic parameters available
    for par_pre::Pair{Symbol, Any} ∈ synapses.pre.parameters
        alias::Symbol = Meta.parse("pre_" * string(par_pre[1]));
        if isa(par_pre[2], Vector)
            synapses.__parameters[alias] = par_pre[2][synapses.__i];
        elseif isa(par_pre[2], Number)
            synapses.__parameters[alias] = par_pre[2] .* ones(synapses.__N);
        else
            @assert false "Spike::Synapses::step(): Could not broadcast parameter `" * string(alias) * "` of unsuppoted type `" * string(typeof(par_pre[2])) * "`.";
        end
    end

    # make renamed postsynaptic parameters available
    for par_post::Pair{Symbol, Any} ∈ synapses.post.parameters
        alias::Symbol = Meta.parse("post_" * string(par_post[1]));
        if isa(par_post[2], Vector)
            synapses.__parameters[alias] = par_post[2][synapses.__j];
        elseif isa(par_post[2], Number)
            synapses.__parameters[alias] = par_post[2] .* ones(synapses.__N);
        else
            @assert false "Spike::Synapses::step(): Could not broadcast parameter `" * string(alias) * "` of unsupported type `" * string(typeof(par_post[2])) * "`.";
        end
    end

    # solve normal equations
    for eq::Pair{Symbol, Expr} ∈ synapses.__normeqs
        synapses.__parameters[eq[1]] = eval(interpolate_from_dict(eq[2], synapses.__parameters));
    end

    # solve differential equations
    for eq::Pair{Symbol, Expr} ∈ synapses.__diffeqs
        synapses.__parameters[eq[1]] = synapses.method(sym = eq[1], eq = eq[2], par = synapses.__parameters, dt = dt, t = t);
    end

    # evaluate pre-synaptic events that synapses are sensitive to
    for preeq::Pair{Symbol, Dict{Symbol, Expr}} ∈ synapses.__preeqs
        if any(synapses.pre.__eventlog[preeq[1]])
            mask_pre::Vector{Bool} = synapses.pre.__eventlog[preeq[1]];
            
            for eq::Pair{Symbol, Expr} ∈ preeq[2]
                str_sym::String = string(eq[1]);
                sender::Vector{Int} = collect(1:synapses.pre.N)[mask_pre];
                receiver::Vector{Vector{Int}} = synapses.__M_pre[sender];
                
                for single::Tuple{Int, Vector{Int}} ∈ zip(sender, receiver)
                    for rec::Int ∈ single[2]
                        indx::Int = collect(1:synapses.__N)[(single[1] .== synapses.__i) .&& (rec .== synapses.__j)][1];
                        ps::Dict{Symbol, Any} = Dict()

                        for p::Pair{Symbol, Any} ∈ synapses.__parameters
                            if size(p[2], 1) > 1
                                ps[p[1]] = synapses.__parameters[p[1]][indx];
                            elseif size(p[2], 1) == 1
                                ps[p[1]] = p[2];
                            end
                        end

                        in_sym::Symbol = Symbol();
                        if length(str_sym) >= 6 && str_sym[1:5] == "post_"
                            in_sym = Symbol(str_sym[6:end]);
                            synapses.post.parameters[in_sym][rec] = eval(interpolate_from_dict(eq[2], ps));
                            synapses.__parameters[eq[1]][rec .== synapses.__j] .= synapses.post.parameters[in_sym][rec];
                        elseif length(str_sym >= 5) && str_sym[1:4] == "pre_"
                            in_sym = Symbol(str_sym[5:end]);
                            synapses.pre.parameters[in_sym][rec] = eval(interpolate_from_dict(eq[2], ps));
                            synapses.__parameters[eq[1]][rec .== synapses.__i] .= synapses.pre.parameters[in_sym][rec];
                        end
                    end
                end
            end
        end
    end

    # evaluate post-synaptic events that synapses are sensitive to
    for posteq::Pair{Symbol, Dict{Symbol, Expr}} ∈ synapses.__posteqs
        if any(synapses.post.__eventlog[posteq[1]])
            mask_post::Vector{Bool} = synapses.post.__eventlog[posteq[1]];
            
            for eq::Pair{Symbol, Expr} ∈ posteq[2]
                str_sym::String = string(eq[1]);
                sender::Vector{Int} = collect(1:synapses.post.N)[mask_post];
                receiver::Vector{Vector{Int}} = synapses.__M_post[sender];
                
                for single::Tuple{Int, Vector{Int}} ∈ zip(sender, receiver)
                    for rec::Int ∈ single[2]
                        indx::Int = collect(1:synapses.__N)[(single[1] .== synapses.__j) .&& (rec .== synapses.__i)][1];
                        ps::Dict{Symbol, Any} = Dict()

                        for p::Pair{Symbol, Any} ∈ synapses.__parameters
                            if size(p[2], 1) > 1
                                ps[p[1]] = synapses.__parameters[p[1]][indx];
                            elseif size(p[2], 1) == 1
                                ps[p[1]] = p[2];
                            end
                        end

                        in_sym::Symbol = Symbol();
                        if length(str_sym) >= 6 && str_sym[1:5] == "post_"
                            in_sym = Symbol(str_sym[6:end]);
                            synapses.post.parameters[in_sym][rec] = eval(interpolate_from_dict(eq[2], ps));
                            synapses.__parameters[eq[1]][rec .== synapses.__j] .= synapses.post.parameters[in_sym][rec];
                        elseif length(str_sym >= 5) && str_sym[1:4] == "pre_"
                            in_sym = Symbol(str_sym[5:end]);
                            synapses.pre.parameters[in_sym][rec] = eval(interpolate_from_dict(eq[2], ps));
                            synapses.__parameters[eq[1]][rec .== synapses.__i] .= synapses.pre.parameters[in_sym][rec];
                        end
                    end
                end
            end
        end
    end
end