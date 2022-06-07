"""
Provides functions to build model components before running the model such that
some operations do not have to be computed at runtime. This is particularly
the case for analysing expressions.
"""

using Parameters;

include("Neurons.jl");
include("Synapses.jl");
include("Operations.jl");
include("Monitors.jl");
include("Expressions.jl");

"""
    collect_and_test!(symbols::Vector{Symbol},
                      dict::Dict{Symbol, Expr})

Traverses a dictionary of symbol => expression pairs, collects symbols into references object and tests equations with dummy parameters.

    INPUTS:
        symbols::Vector{Symbol}     -   Referenced list of symbols (will be updated!).
        dict::Dict{Symbol, Expr}    -   Dictionary of symbols and equations.
"""
function collect_and_test!(symbols::Vector{Symbol}, dict::Dict{Symbol, Expr})
    for eq::Pair{Symbol, Expr} ∈ dict
        push!(symbols, eq[1]);
        eq[2] |> push_symbols!(symbols);
        dummies::Dict{Symbol, Any} = dummy_parameters(unique!(symbols));
        if haskey(dummies, :N)
            dummies[:N] = 2;
        end
        test_expression(; expr = eq[2], parameters = dummies);
    end
end

"""
    defaults_for_missing(required::Vector{Symbol}, 
                         parameters::Dict{Symbol, Any}; 
                         default::Float64 = 0.0)::Dict{Symbol, Any}

In-paints missing required symbols in a dictionary with a default value.

    INPUTS:
        required::Vector{Symbol}        -   Required symbols (keys).
        parameters::Dict{Symbol, Any}   -   Parameter dictionary.
        N::Int                          -   Number of values required.
        default::Float64                -   Value to insert.
    
    OUTPUTS:
        parameters::Dict{Symbol, Any}   -   Updated parameter dictionary.
"""
function defaults_for_missing(required::Vector{Symbol}, parameters::Dict{Symbol, Any}, N::Int; default::Float64 = 0.0)::Dict{Symbol, Any}
    for required_symbol::Symbol ∈ required
        if !haskey(parameters, required_symbol)
            parameters[required_symbol] = ones(N) .* default;
        end
    end

    parameters;
end

"""
    build(target::NeuronGroup)::NeuronGroup

Builds a group of neurons such that their expressions may be evaluated more directly. Also performs safety checks. This is an internal function and should not be called manually.

    INPUTS:
        target::NeuronGroup     -   Group of neurons to build.

    OUTPUTS:
        target::NeuronGroup     -   Built group of neurons.
"""
@fastmath function build(target::NeuronGroup)::NeuronGroup
    # don't rebuild
    if target.__built == true
        return target;
    end

    # standard parameters
    symbols::Vector{Symbol} = Symbol[:N, :dt, :t];

    # build internal state expressions
    target.__normeqs, target.__diffeqs = categorise_subexpressions(target.eq);
    collect_and_test!(symbols, merge(target.__normeqs, target.__diffeqs));

    # build internal event expressions
    for event::Pair{Symbol, Tuple{Expr, Expr}} ∈ target.events
        target.__eventeqs[event[1]], _ = categorise_subexpressions(event[2][2]);
        collect_and_test!(symbols, Dict(:N => event[2][1]));
        collect_and_test!(symbols, target.__eventeqs[event[1]]);
    end

    # create default parameters if necessary
    target.parameters = defaults_for_missing(symbols, target.parameters, target.N);

    # check parameters
    for par::Pair{Symbol, Any} ∈ target.parameters
        @assert size(par[2], 1) ∈ [1, target.N] "\nSpike::Build::build(::NeuronGroup):\nDetected parameter `" * string(par[1]) * "` with first dimension not in [1, " * string(target.N) * "].";
    end

    # finalise
    target.__built = true;
    target;
end

"""
    build(target::Synapses)::Synapses

Builds synapses between two groups of neurons such that their equations can be evaluated more directly. Also performs safety checks and handles evaluation of connectivity matrices as well as internal parametrisation. This is an internal function and should not be called manually.

    INPUTS:
        target::Synapses    -   Synapses to build.

    OUTPUTS:
        target::Synapses    -   Built group of synapses.
"""
@fastmath function build(target::Synapses)::Synapses
    # don't rebuild
    if target.__built == true
        return target;
    end

    # standard parameters
    symbols::Vector{Symbol} = Symbol[:N, :N_pre, :N_post, :dt, :t];

    # build internal state equations
    target.__normeqs, target.__diffeqs = categorise_subexpressions(target.eq);
    collect_and_test!(symbols, merge(target.__normeqs, target.__diffeqs));

    # build internal presynaptic event hooks and equations
    for event_pre::Pair{Symbol, Expr} ∈ target.on_pre
        target.__preeqs[event_pre[1]], _ = categorise_subexpressions(event_pre[2]);
        collect_and_test!(symbols, target.__preeqs[event_pre[1]]);
    end

    # build internal postsynaptic event hooks and equations
    for event_post::Pair{Symbol, Expr} ∈ target.on_post
        target.__posteqs[event_post[1]], _ = categorise_subexpressions(event_post[2]);
        collect_and_test!(symbols, target.__posteqs[event_post[1]]);
    end

    # create full connectivity vectors
    pre_N::Int, post_N::Int = (target.pre.N, target.post.N);
    i::Vector{Int}, j::Vector{Int} = (repeat(1:pre_N, inner = post_N), repeat(1:post_N, outer = pre_N));
    mask::Vector{Bool} = ones(size(i, 1));

    # create conditional mask from expression and apply to vectors
    if target.cond != :()
        data::Dict{Symbol, Any} = Dict(:pre_N => pre_N, :post_N => post_N, :i => i, :j => j);
        mask = eval(interpolate_from_dict(target.cond, data));
        i, j = (i[mask], j[mask]);
    end

    # apply probabilistic mask to vectors
    if target.prob > 0.0 && target.prob < 1.0
        mask = rand(size(i, 1)) .> (1 - target.prob);
        i, j = (i[mask], j[mask]);
    end

    # set final connectivity vectors
    target.i, target.j = (i, j);

    # create empty forwards connectivity vectors
    target.M_pre = Vector{Vector{Int}}(undef, target.pre.N);
    for pre_i::Int ∈ 1:size(target.M_pre, 1)
        target.M_pre[pre_i] = Int[];
    end

    # create empty backwards connectivity vectors
    target.M_post = Vector{Vector{Int}}(undef, target.post.N);
    for post_i::Int ∈ 1:size(target.M_post, 1)
        target.M_post[post_i] = Int[];
    end
    
    # fill matrices
    for ij::Tuple{Int, Int} ∈ zip(i, j)
        push!(target.M_pre[ij[1]], ij[2]);
        push!(target.M_post[ij[2]], ij[1]);
    end

    # setup new parameters
    parameters::Dict{Symbol, Any} = Dict();

    # insert pre-specified parameters
    for par::Pair{Symbol, Any} ∈ target.parameters
        if isa(par[2], Function)
            parameters[par[1]] = par[2](size(target.i, 1))
        elseif isa(par[2], Number)
            parameters[par[1]] = ones(size(target.i, 1)) .* par[2];
        else
            @assert false "\nSpike::Build::build(::Synapses):\nDetected parameter `" * string(par[1]) * "` of unsupported type `" * string(typeof(par[2])) * "`.";
        end
    end

    # setup constant internal parameters
    target.N = size(target.i, 1);
    parameters[:N] = target.N;
    parameters[:pre_N] = target.pre.N;
    parameters[:post_N] = target.post.N;

    # make renamed presynaptic parameters available
    for par_pre::Pair{Symbol, Any} ∈ target.pre.parameters
        alias::Symbol = Meta.parse("pre_" * string(par_pre[1]));
        if isa(par_pre[2], Vector)
            parameters[alias] = par_pre[2][target.i];
        elseif isa(par_pre[2], Number)
            parameters[alias] = par_pre[2] .* ones(target.N);
        else
            @assert false "\nSpike::Build::build(::Synapses):\nCould not broadcast parameter `" * string(alias) * "` of unsuppoted type `" * string(typeof(par_pre[2])) * "`. Allowed = [Vector, Number].";
        end
    end

    # make renamed postsynaptic parameters available
    for par_post::Pair{Symbol, Any} ∈ target.post.parameters
        alias::Symbol = Meta.parse("post_" * string(par_post[1]));
        if isa(par_post[2], Vector)
            parameters[alias] = par_post[2][target.j];
        elseif isa(par_post[2], Number)
            parameters[alias] = par_post[2] .* ones(target.N);
        else
            @assert false "\nSpike::Build::build(::Synapses):\nCould not broadcast parameter `" * string(alias) * "` of unsupported type `" * string(typeof(par_post[2])) * "`. Allowed = [Vector, Number].";
        end
    end

    # create default parameters if necessary
    target.parameters = defaults_for_missing(symbols, parameters, target.N);

    # finalise
    target.__built = true;
    target;
end

"""
    build(target::Operation)::Operation

Builds an operation into the model, which will be executed at prespecified intervals. This is an internal function and should not be called manually.

    INPUTS:
        target::Operation       -   Operation to build.

    OUTPUTS:
        target::Operation       -   Built operation.
"""
function build(target::Operation)::Operation
    # don't rebuild
    if target.__built == true
        return target;
    end

    # check cycle
    @assert target.cycle ∈ ["pre", "post"] "\nSpike::Build::build(::Operation):\nReceived unsupported value for cycle = `" * target.cycle * ". Allowed = [pre, post].`";

    # check timing
    @assert target.every >= 1e-3 "\nSpike::Build::build(::Operation):\nReceived unsupported value for every = `" * string(target.every) * "`. Allowed >= 1e-3.";
    
    # finalise
    target.__built = true;
    target;
end

"""
    build(target::StateMonitor)::StateMonitor

Builds a state monitor into the model that will be updated periodically. This is an internal function and should not be called manually.

    INPUTS:
        target::StateMonitor    -   StateMonitor to build.
    
    OUTPUTS:
        target::StateMonitor    -   Built monitor.
"""
function build(target::StateMonitor)::StateMonitor
    # don't rebuild
    if target.__built == true
        return target;
    end
    
    # check object type
    @assert typeof(target.obj) ∈ [NeuronGroup, Synapses] "\nSpike::Build::build(::StateMonitor):\nReceived unsupported obj of type `" * string(typeof(target.obj)) * "`. Allowed [NeuronGroup, Synapses].";

    # check timing
    @assert target.every >= 1e-3 "\nSpike::Build::build(::StateMonitor):\nReceived unsupported value for every = `" * string(target.every) * "`. Allowed >= 1e-3.";

    # initialise
    for var::Symbol ∈ target.vars
        N::Int = target.obj.N;
        target.states[var] = Array{Float64}(undef, N, 0);
    end

    # finalise
    target.__built = true;
    target;
end

"""
    build(target::EventMonitor)::EventMonitor

Builds an event monitor into the model that will be updated as events occur. This is an internal function and should not be called manually.

    INPUTS:
        target::EventMonitor    -   EventMonitor to build.

    OUTPUTS:
        target::EventMonitor    -   Built monitor.
"""
function build(target::EventMonitor)::EventMonitor
    # don't rebuild
    if target.__built == true
        return target;
    end

    # finalise
    target.__built = true;
    target;
end