"""
Provides the main structure for models, standard operations for creating models
magically before running them and, finally, running simulations of them.
"""

using Parameters;

include("Magic.jl");
include("Neurons.jl");
include("Synapses.jl");
include("Operations.jl");
include("Monitors.jl");
include("Build.jl");

@with_kw mutable struct Model
    """
    Main structure for creating a model. Note that, for convenience, this structure does
    not have to be created yourself. If you can provide adeqaute scope and have called
    Spike::Magic::cast_magic() to start the scope, Spike::Model::run(; ...) will magically 
    created the model for you (and return it after simulation).

    INPUTS:
        Neurons::Dict{Symbol, NeuronGroup}          -   All neuron groups in the model. (default = Dict())
        Synapses::Dict{Symbol, Synapses}            -   All synapse groups in the model. (default = Dict())
        Operations::Dict{Symbol, Operation}         -   All operations in the model. (default = Dict())
        StateMonitors::Dict{Symbol, StateMonitor}   -   All state monitors in the model. (default = Dict())
        EventMonitors::Dict{Symbol, EventMonitor}   -   All event monitors in the model. (default = Dict())
        verbose::Bool                               -   Verbose updates? (default = true)
    """

    Neurons::Dict{Symbol, NeuronGroup} = Dict();
    Synapses::Dict{Symbol, Synapses} = Dict();
    Operations::Dict{Symbol, Operation} = Dict();
    StateMonitors::Dict{Symbol, StateMonitor} = Dict();
    EventMonitors::Dict{Symbol, EventMonitor} = Dict();
    verbose::Bool = true;
end

function status(model::Model, status::String; e::String = "\n")
    """
    Logging function at runtime.
    
    INPUTS:
        model::Model        -   The model object.
        status::String      -   The status.
        e::String           -   End-of-line characters.
    """

    if model.verbose == false
        return;
    end

    print(status * e);
end

@fastmath function run(model::Model; T::Float64, dt::Float64 = 1e-3)::Model
    """
    Main simulation entry point.

    INPUTS:
        model::Model    -   The model to run.
        T::Float64      -   Total time to run the model for.
        dt::Float64     -   Time step size. (default = 1e-3)
    
    OUTPUTS:
        model::Model    -   Self
    """

    # build the model
    status(model, "Building model components...");

    for tok::Symbol ∈ fieldnames(Model)
        if isa(getproperty(model, tok), Dict)
            for obj::Pair{Any, Any} ∈ getproperty(model, tok)
                if typeof(obj[1]) != Symbol || !isa(obj[2], SpikeObject)
                    continue;
                end

                getproperty(model, tok)[obj[1]] = build(getproperty(model, tok)[obj[1]]);
            end
        end
    end

    # run a simulation
    for t::Float64 ∈ collect(0.0:dt:T)
        progress::String = repeat("-", Int(round((t / T) * 25))) * repeat(" ", Int(round(25 - (t / T) * 25)));
        status(model, "Simulating...\t[" * progress * "]  \t" * string(round((t / T)*100)) * "%\t", e = "\r");

        for op::Pair{Symbol, Operation} ∈ model.Operations
            step(model.Operations[op[1]]; t = t, dt = dt, cycle = "pre");
        end

        for tok::Symbol ∈ fieldnames(Model)
            if isa(getproperty(model, tok), Dict)
                for obj::Pair{Any, Any} ∈ getproperty(model, tok)
                    if typeof(obj[1]) != Symbol || !isa(obj[2], SpikeObject) || typeof(obj[2]) == Operation
                        continue;
                    end
    
                    step(getproperty(model, tok)[obj[1]]; t = t, dt = dt);
                end
            end
        end

        for op::Pair{Symbol, Operation} ∈ model.Operations
            step(model.Operations[op[1]]; t = t, dt = dt, cycle = "post");
        end
    end

    status(model, "\nDone.");

    model;
end

function run(; T::Float64, dt::Float64 = 1e-3, magic_obj::DataType = SpikeObject, magic_tar::Module = Main, verbose::Bool = true)::Model
    """
    Entry point for magic networks; creates the model and yields it
    to simulation entry point.

    INPUTS:
        T::Float64              -   Total time to run the model for.
        dt::Float64             -   Time step size.
        magic_obj::DataType     -   DataType that identifies magic objects. (default = SpikeObject)
        magic_tar::Module       -   Module where magic is utilised. (default = Main)
        verbose::Bool           -   Verbose updates? (default = true)
    
    OUTPUTS:
        model::Model            -   Self
    """

    # create a magic model
    model::Model = create_magic_model(; magic_obj = magic_obj, magic_tar = magic_tar);
    model.verbose = verbose;

    # yield to simulation
    status(model, "Created magic model.");
    model = run(model; T = T, dt = dt);
    model;
end