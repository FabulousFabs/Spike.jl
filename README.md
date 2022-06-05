# Spike.jl

[![Build Status](https://github.com/FabulousFabs/Spike.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/FabulousFabs/Spike.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package allows simulations of spiking neural networks in Julia. Experimenting with spiking networks typically means working with many changing parts in equations (e.g., most projects will differ slightly in the models used), which can create a lot of unnecessary overhead in (re-)implementing models.

Spike.jl addresses this issue by providing extreme flexibility in model specification by utilising a heavy meta-programming approach that allows models to be specified through expressions, whereupon Spike.jl will handle everything else behind the scenes.


**Note that** this is an alpha version and bugs do in all likelihood exist. **Note also that** Spike.jl requires Julia >= 1.7.

## How do I get it?
Currently, you will have to clone the repository and build it yourself.

```Julia
] add https://github.com/FabulousFabs/Spike.jl.git
```

## How do I use it?
Here's a very simple "Hello World!" example of creating a model in Spike.jl:

```Julia
using Spike;
using Plots;

# enable magic model creation
cast_magic();

# create random spikers
spikers = NeuronGroup(N = 10,
                      eq = :(v_t = rand(N);),
                      method = euler,
                      events = Dict(:spike => (:(v .> 0.95), :(v = v_r;))),
                      parameters = Dict(:v => zeros(10),
                                        :v_r => zeros(10)));

# monitor events
monitor = EventMonitor(obj = spikers, event = :spike);

# run a magic simulation
run(; T = 1.0, dt = 1e-3);

# plot
scatter(monitor.t, monitor.i, markershape = :vline)
```

Okay, but what about using it in a real project? I'll want to embed it in my own simulation code to make it control an agent, for example. Here's a sketch:

```Julia
using Spike;
using Parameters;

@with_kw mutable struct MySimulation
  env::Environment
  agent::Agent
  ...
  model::Model = Model()
end

mysim = MySimulation(...)
mysim.model.Neurons[:input_layer] = NeuronGroup(...)
mysim.model.Neurons[:decision_layer] = NeuronGroup(...)
mysim.model.Synapses[:feedforward] = Synapses(...)

function my_callback()
  println("I can do some magic here to control the agent and environment!");
end

mysim.model.Operations[:my_operation] = Operation(op = my_callback,
                                                  every = 1e-3,
                                                  cycle = "post");

run(mysim.model; T = 1.0, dt = 1e-3);
```

As you can imagine, there are many ways to implement this, but the general formula will always involve using an operation to hook into the step-cycle (either pre- or post-step) of the model at runtime.

For more information, please see the examples and documentation.

## @todo:
- Implement model I/O for saving and loading
- Write examples and documentation
- Add unit tests
