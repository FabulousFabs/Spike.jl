"""
Provides the wrapper for custom operations to be added to a model as well
as basic timing control.
"""

using Parameters;

include("Magic.jl");

"""
    Operation <: SpikeObject

Wrapper structure for custom operations to be run during simulation.

    INPUTS:
        op::Function        -   External function to be called.
        every::Float64      -   Delay between calls. (default = 1e-3)
        cycle::String       -   Timing within cycle. (default = "pre")
        __built::Bool       -   (Internal) Has this operation been built? (default = false)
        __last::Float64     -   (Internal) Last operation. (default = -Inf)
"""
@with_kw mutable struct Operation <: SpikeObject
    op::Function
    every::Float64 = 1e-3
    cycle::String = "pre"

    __built::Bool = false
    __last::Float64 = -Inf
end

"""
    step(operation::Operation; dt::Float64, 
                               t::Float64, 
                               cycle::String)::Operation

Performs one time step of an operation. This is an internal function and should not be called manually.

    INPUTS:
        operation::Operation    -   Operation to perform time step on.
        dt::Float64             -   Time step size.
        t::Float64              -   Current time.
        cycle::String           -   Current cycle.
    
    OUTPUTS:
        operation::Operation    -   Self
"""
function step(operation::Operation; dt::Float64, t::Float64, cycle::String)::Operation
    # run operation if in cycle and timing
    if cycle == operation.cycle
        if operation.__last + operation.every <= t
            operation.op()
            operation.__last = t;
        end
    end

    operation;
end

