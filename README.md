# <img src="./docs/src/assets/logo.png" alt="" width=90 /> Spike.jl

[![Build Status](https://github.com/FabulousFabs/Spike.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/FabulousFabs/Spike.jl/actions/workflows/CI.yml?query=branch%3Amain)&nbsp;&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;&nbsp;**[DOCUMENTATION](https://fabulousfabs.github.io/Spike.jl.Docs/)**

This package allows simulations of spiking neural networks in Julia. Experimenting with spiking networks typically means working with many changing parts in equations (e.g., most projects will differ slightly in the models used), which can create a lot of unnecessary overhead in (re-)implementing models.

Spike.jl addresses this issue by providing extreme flexibility in model specification by utilising a heavy meta-programming approach that allows models to be specified through expressions, whereupon Spike.jl will handle everything else behind the scenes.

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

For more information, please see the examples and documentation.
