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
    more directly.

    INPUTS:
        target::NeuronGroup     -   Group of neurons to build.
    
    OUTPUTS:
        target::NeuronGroup     -   Built group of neurons.
    """

    target.__normeqs, target.__diffeqs = categorise_subexpressions(target.eq);
    target.__built = true;
    target;
end