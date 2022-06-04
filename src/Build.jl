"""
Provides functions to build model components before running the model such that
some operations do not have to be computed at runtime. This is particularly
the case for analysing expressions.
"""

using Parameters;

include("Neurons.jl");
include("Synapses.jl");
include("Expressions.jl");

@fastmath function build(target::NeuronGroup)::NeuronGroup
    """
    Builds a group of neurons such that their expressions may be evaluated
    more directly and performs safety checks.

    INPUTS:
        target::NeuronGroup     -   Group of neurons to build.
    
    OUTPUTS:
        target::NeuronGroup     -   Built group of neurons.
    """

    target.__normeqs, target.__diffeqs = categorise_subexpressions(target.eq);

    for event::Pair{Symbol, Tuple{Expr, Expr}} ∈ target.events
        target.__eventeqs[event[1]], _ = categorise_subexpressions(event[2][2]);
    end

    for par::Pair{Symbol, Any} ∈ target.parameters
        @assert size(par[2], 1) ∈ [1, target.N] "Spike::Build::build(::NeuronGroup): Detected parameter `" * string(par[1]) * "` with first dimension not in [1, " * string(target.N) * "].";
    end

    target.__built = true;
    target;
end

@fastmath function build(target::Synapses)::Synapses
    """
    Builds synapses between two groups of neurons such that their equations can be
    evaluated more directly, performs safety checks and evaluates the connectivity
    matrix and, subsequently, internal parametrisation.

    INPUTS:
        target::Synapses    -   Synapses to build.
    
    OUTPUTS:
        target::Synapses    -   Built group of synapses.
    """

    target.__normeqs, target.__diffeqs = categorise_subexpressions(target.eq);

    for event_pre::Pair{Symbol, Expr} ∈ target.on_pre
        target.__preeqs[event_pre[1]], _ = categorise_subexpressions(event_pre[2]);
    end

    for event_post::Pair{Symbol, Expr} ∈ target.on_post
        target.__posteqs[event_post[1]], _ = categorise_subexpressions(event_post[2]);
    end

    pre_N::Int, post_N::Int = (target.pre.N, target.post.N);
    i::Vector{Int}, j::Vector{Int} = (repeat(1:pre_N, inner = post_N), repeat(1:post_N, outer = pre_N));
    mask::Vector{Bool} = ones(size(i, 1));

    if target.cond != :()
        data::Dict{Symbol, Any} = Dict(:pre_N => pre_N, :post_N => post_N, :i => i, :j => j);
        mask = eval(interpolate_from_dict(target.cond, data));
        i, j = (i[mask], j[mask]);
    end

    if target.prob > 0.0 && target.prob < 1.0
        mask = rand(size(i, 1)) .> target.prob;
        i, j = (i[mask], j[mask]);
    end

    target.__i, target.__j = (i, j);

    target.__M_pre = Vector{Vector{Int}}(undef, target.pre.N);
    for pre_i::Int ∈ 1:size(target.__M_pre, 1)
        target.__M_pre[pre_i] = Int[];
    end
    target.__M_post = Vector{Vector{Int}}(undef, target.post.N);
    for post_i::Int ∈ 1:size(target.__M_post, 1)
        target.__M_post[post_i] = Int[];
    end
    
    for ij::Tuple{Int, Int} ∈ zip(i, j)
        push!(target.__M_pre[ij[1]], ij[2]);
        push!(target.__M_post[ij[2]], ij[1]);
    end

    for par::Pair{Symbol, Any} ∈ target.parameters
        if isa(par[2], Function)
            target.__parameters[par[1]] = par[2](size(target.__i, 1))
        elseif isa(par[2], Number)
            target.__parameters[par[1]] = ones(size(target.__i, 1)) .* par[2];
        else
            @assert false "Spike::Build::build(::Synapses): Detected parameter `" * string(par[1]) * "` of unsupported type `" * string(typeof(par[2])) * "`.";
        end
    end

    target.__N = size(target.__i, 1);

    target;
end