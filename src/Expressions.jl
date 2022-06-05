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
    global __regex_subexpression_eq_norm, __regex_subexpression_eq_diff, __regex_subexpression_eq_assi;

    eqs_norm::Dict{Symbol, Expr} = Dict()
    eqs_diff::Dict{Symbol, Expr} = Dict() 

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
                    @assert false "Spike::Expressions::categorise_subexpressions(): Could not categorise subexpression:\n" * str_sy * "\nIs this a parameter instead?";
                end

                in_eq = false;
                in_sy = :Placeholder;
            end
        end
    end

    (eqs_norm, eqs_diff);
end