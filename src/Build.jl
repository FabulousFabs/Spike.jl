"""
Provides functions to build model components before running the model such that
some operations do not have to be computed at runtime. This is particularly
the case for analysing expressions.
"""

using Parameters;

include("Neurons.jl");
include("Synapses.jl");
include("Expressions.jl");

function build(target::NeuronGroup)::NeuronGroup
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