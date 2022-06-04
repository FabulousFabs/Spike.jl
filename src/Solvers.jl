"""
Provides functions for solving differential equations.
"""

function rk2(; sym::Symbol, eq::Expr, par::Dict{Symbol, Any}, dt::Float64, t::Float64)::Any
    """
    Second order Runge-Kutta time stepper for differential equations.

    INPUTS:
        sym::Symbol                 -   Symbol of the parameters to be differentiated.
        eq::Expr                    -   Expression of the differential equation.
        par::Dict{Symbol, Any}      -   Paramters for the differential equation.
        dt::Float64                 -   Step size.
        t::Float64                  -   Time at step zero.
    
    OUTPUTS:
        phi_n_plus_1::Any           -   Approximate value at step one.
    """
    
    par[sym] .+= (dt .* eval(interpolate_from_dict(eq, par))) ./ 2;
    par[:t] += dt ./ 2;
    par[sym] .+ (dt .* eval(interpolate_from_dict(eq, par)));
end

function euler(; sym::Symbol, eq::Expr, par::Dict{Symbol, Any}, dt::Float64, t::Float64)::Any
    """
    Euler's method for differential equations.

    INPUTS:
        sym::Symbol                 -   Symbol of the parameters to be differentiated.
        eq::Expr                    -   Expression of the differential equation.
        par::Dict{Symbol, Any}      -   Paramters for the differential equation.
        dt::Float64                 -   Step size.
        t::Float64                  -   Time at step zero.
    
    OUTPUTS:
        phi_n_plus_1::Any           -   Approximate value at step one.
    """

    par[sym] .+ (dt .* eval(interpolate_from_dict(eq, par)));
end

function heun(; sym::Symbol, eq::Expr, par::Dict{Symbol, Any}, dt::Float64, t::Float64)::Any
    """
    Heun's method for differential equations.

    INPUTS:
        sym::Symbol                 -   Symbol of the parameters to be differentiated.
        eq::Expr                    -   Expression of the differential equation.
        par::Dict{Symbol, Any}      -   Paramters for the differential equation.
        dt::Float64                 -   Step size.
        t::Float64                  -   Time at step zero.
    
    OUTPUTS:
        phi_n_plus_1::Any           -   Approximate value at step one.
    """

    par_dy::Dict{Symbol, Any} = par;
    par_dy[:t] += dt;
    par_dy[sym] .+= dt .* eval(interpolate_from_dict(eq, par));
    par[sym] .+ (dt ./ 2) .* (par_dy[sym] .- par[sym] .+ eval(interpolate_from_dict(eq, par_dy)));
end