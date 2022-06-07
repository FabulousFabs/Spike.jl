"""
Provides expression functionality to ease their use as equations in model specification.
"""

"""
    interpolate_from_dict(expr::Expr, 
                          dict::Dict{Symbol, Any})::Any

Set of functions to allow for data to be supplied easily for evaluating expressions.
"""
interpolate_from_dict(expr::Expr, dict::Dict{Symbol, Any}) = Expr(expr.head, interpolate_from_dict.(expr.args, Ref(dict))...);
interpolate_from_dict(expr::Symbol, dict::Dict{Symbol, Any}) = get(dict, expr, expr);
interpolate_from_dict(expr::Any, dict::Dict{Symbol, Any}) = expr;

"""
    __regex_subexpression_eq_norm

Regular expression that defines the syntax of normal equations in [`categorise_subexpressions`](@ref). This is an internal variable and should not be changed.
"""
__regex_subexpression_eq_norm = r"([a-zA-Z0-9\_]+)\_t";

"""
    __regex_subexpression_eq_diff

Regular expression that defines the syntax of differential equations in [`categorise_subexpressions`](@ref). This is an internal variable and should not be changed.
"""
__regex_subexpression_eq_diff = r"d([a-zA-Z0-9\_]+)\_dt";

"""
    __regex_subexpression_eq_assi

Regular expression that defines the syntax of assignments in [`categorise_subexpressions`](@ref). This is an internal variable and should not be changed.
"""
__regex_subexpression_eq_assi = r"([a-zA-Z0-9\_]+)";

"""
    categorise_subexpressions(expr::Expr)::Tuple{Dict{Symbol, Expr}, Dict{Symbol, Expr}}

Takes a set of equations formulated as an expression and categorises by whether or not they need to be differentiated at runtime. Currently, all equations are implicitly interpreted with respect to time. Typical syntax of equations should be:

```julia
:(I_t = A ./ 2 .+ A ./ 2 .* sin(2 * π * t);
  dv_dt = (.-(v .- v_rest) .+ I) ./ 𝜏;
  A = rand(N);)
```

which would create a regular function I(t), a differential equation dv/dt and an assignment of A for the model at runtime.

    INPUTS:
        expr::Expr                                                              -   Expression to evaluate
    
    OUTPUTS:
        (eqs_norm, eqs_diff)::Tuple{Dict{Symbol, Expr}, Dict{Symbol, Expr}}     -   Tuple of normal and differential equations.
"""
function categorise_subexpressions(expr::Expr)::Tuple{Dict{Symbol, Expr}, Dict{Symbol, Expr}}
    # get regexes
    global __regex_subexpression_eq_norm, __regex_subexpression_eq_diff, __regex_subexpression_eq_assi;

    # setup dictionaries
    eqs_norm::Dict{Symbol, Expr} = Dict()
    eqs_diff::Dict{Symbol, Expr} = Dict() 

    # find subexpressions (requires ; termination!)
    for sub_expr::Any ∈ expr.args
        if typeof(sub_expr) != Expr
            continue;
        end

        in_eq::Bool = false;
        in_sy::Symbol = :Placeholder;

        for i::Int ∈ 1:size(sub_expr.args, 1)
            if in_eq == false && typeof(sub_expr.args[i]) == Symbol
                in_eq = true;
                in_sy = sub_expr.args[i];
            elseif in_eq == true
                str_sy::String = string(in_sy);

                if occursin(__regex_subexpression_eq_diff, str_sy) == true
                    eqs_diff[Meta.parse(match(__regex_subexpression_eq_diff, str_sy)[1])] = sub_expr.args[i];
                elseif occursin(__regex_subexpression_eq_norm, str_sy) == true
                    eqs_norm[Meta.parse(match(__regex_subexpression_eq_norm, str_sy)[1])] = sub_expr.args[i];
                elseif occursin(__regex_subexpression_eq_assi, str_sy) == true
                    x::Any = sub_expr.args[i];
                    eqs_norm[Meta.parse(match(__regex_subexpression_eq_assi, str_sy)[1])] = :($x * 1);
                else
                    @assert false "\nSpike::Expressions::categorise_subexpressions():\nCould not categorise subexpression:\n" * str_sy * "\nAre you missing a semicolon anywhere in `" * string(expr) * "` to terminate an expression?";
                end

                in_eq = false;
                in_sy = :Placeholder;
            end
        end
    end

    (eqs_norm, eqs_diff);
end

"""
    push_symbols!(symbols::Vector{Symbols})

Recursively walks an expression to find all its symbols. Call using pipe operator, e.g.: 
```Julia
symbols::Vector{Symbol} = Symbol[];
:(a + b + c) |> push_symbols!(symbols);
```

    INPUTS:
        symbols::Vector{Symbol}     -   Reference of where to store symbols.
"""
push_symbols!(symbols) = expr -> begin
    expr isa Symbol && push!(symbols, expr);
    expr isa Expr && expr.head == :call && map(push_symbols!(symbols), expr.args[2:end]);
end

"""
    dummy_parameters(symbols::Vector{Symbol}; N::Int)::Dict{Symbol, Any}

Creates a dictionary of dummy parameters for all symbols.

    INPUTS:
        symbols::Vector{Symbol}         -   Vector of symbols.
        N::Int                          -   Number of entries per parameter. (default = 2)
    
    OUTPUTS:
        parameters::Dict{Symbol, Any}   -   Dummy dictionary
"""
function dummy_parameters(symbols::Vector{Symbol}; N::Int = 2)::Dict{Symbol, Any}
    parameters::Dict{Symbol, Any} = Dict([x => ones(N) for x::Symbol ∈ symbols]);
end

"""
    test_expression(; expr::Expr, 
                      parameters::Dict{Symbol, Any})::Bool

Evaluates an expression in a try...catch...-block to generate an error (or output passes = true).

    INPUTS:
        expr::Expr                      -   Expression to test
        paramters::Dict{Symbol, Any}    -   Parameters to use for testing.
    
    OUTPUTS:
        passes::Bool                    -   Does the expression pass?
"""
function test_expression(; expr::Expr, parameters::Dict{Symbol, Any})::Bool
    passes::Bool = false;
    
    try 
        outputs::Vector{Any} = eval(interpolate_from_dict(expr, parameters));
        passes = true;
    catch exc
        @assert false "Spike::Expressions::test_expressions(): Detected a faulty expression `" * string(expr) * "`. Cannot proceed. Error:\n" * string(exc);
    end

    passes;
end

"""
    test_expressions(; expr::Vector{Expr}, 
                       parameters::Dict{Symbol, Any})::Bool

Wrapper for a vector of expressions to be evaluated using the same set of parameters. See `test_expressions(; expr::Expr ...)`.

    INPUTS:
        Vector{Expr}                    -   Vector of expressions to be tested.   
        parameters::Dict{Symbol, Any}   -   Parameters to use for testing.
    
    OUTPUTS:
        all_passes::Bool                -   Do all expressions pass?
"""
function test_expressions(; expr::Vector{Expr}, parameters::Dict{Symbol, Any})::Bool
    all([test_expression(; expr = ex, parameters = parameters) for ex::Expr ∈ expr]);
end