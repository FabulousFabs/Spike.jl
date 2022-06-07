"""
Provides magic network functionality such that models may be built implicitly in the main file.
More concretely, this allows for model specifications such as:

```
    using Spike;

    cast_magic();

    population = LIF(N = 100);
    ...

    model = run(; T = 1.0, dt = 1e-3);
```

whereby the model is collected from the global Main scope before being built and executed.
"""

"""
    SpikeObject

Abstract type definition that enables look-up and matching of magic objects.
"""
abstract type SpikeObject end;

"""
    __cast_magic::Bool

Internal variable to keep track of whether or not magic has been cast. Do not change this manually. Please use [`cast_magic`](@ref) instead.
"""
__cast_magic = false;

"""
    cast_magic()

Enables creation of a magic network that will be built from the module's global scope. This should always be called when trying to utilise magic networks.
"""
function cast_magic()
    global __cast_magic = true;
end

"""
    collect_magic_objects(; magic_obj::DataType = SpikeObject, 
                            magic_tar::Module = Main)::Vector{Symbol}

Collects all magic objects from the module. This is an internal function and should not be called manually.

    INPUTS:
        magic_obj::DataType     -   DataType that identifies magic objects. (default = SpikeObject)
        magic_tar::Module       -   Module where magic is utilised. (default = Main)
    
    OUTPUTS:
        cols::Vector{Symbol}    -   Vector of symbols of magic objects in module.
"""
function collect_magic_objects(; magic_obj::DataType = SpikeObject, magic_tar::Module = Main)::Vector{Symbol}
    objs::Vector{Symbol} = names(magic_tar, all = true);
    cols::Vector{Symbol} = Symbol[];

    for obj::Symbol ∈ objs
        try 
            if isa(getproperty(magic_tar, obj), magic_obj) && obj != :ans
                push!(cols, obj);
            end
        catch exc
            if !isa(exc, UndefVarError)
                @assert false "\nSpike::Magic::collect_magic_objects():\nEncountered an unexpected error of type `" * string(typeof(exc)) * "`:\n" * string(exc);
            end
        end
    end

    cols;
end

"""
    create_magic_model(; magic_obj::DataType = SpikeObject, 
                         magic_tar::Module = Main)::Model

Creates a model from all available magic objects that have corresponding vectorised forms in the standard model structure. This is an internal function and should not be called manually.

    INPUTS:
        magic_obj::DataType     -   DataType that identifies magic objects. (default = SpikeObject)
        magic_tar::Module       -   Module where magic is utilised. (default = Main)

    OUTPUTS:
        model::Model            -   The magic network.
"""
function create_magic_model(; magic_obj::DataType = SpikeObject, magic_tar::Module = Main)::Model
    global __cast_magic;

    @assert __cast_magic == true "\nSpike::Magic::build_magic():\nReceived call while no initial magic was cast.";

    model::Model = Model();
    magic_toks::Vector{Symbol} = Symbol[];
    magic_typs::Vector{DataType} = DataType[];

    for tok::Symbol ∈ fieldnames(Model)
        if isa(getproperty(model, tok), Dict)
            if tok ∉ magic_toks && typeof(getproperty(model, tok)) ∉ magic_typs
                push!(magic_toks, tok);
                push!(magic_typs, typeof(getproperty(model, tok)));
            end
        end
    end

    magic_objs::Vector{Symbol} = collect_magic_objects(; magic_obj = magic_obj, magic_tar = magic_tar);

    for obj::Symbol ∈ magic_objs
        indx::Vector{Int} = findall(x -> x == Dict{Symbol, typeof(getproperty(magic_tar, obj))}, magic_typs);

        @assert size(indx, 1) > 0 "\nSpike::Magic::build_magic():\nA magic object was found that has no vector-type correspondence with Spike::Model::Model().";
        @assert size(indx, 1) == 1 "\nSpike::Magic::build_magic():\nA magic object was found that has multiple vector-type correspondences with Spike::Model::Model().";

        getproperty(model, magic_toks[indx[1]])[obj] = getproperty(magic_tar, obj);
    end

    model;
end