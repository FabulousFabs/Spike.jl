#  Guide
Spike allows quick and easy implementations of spiking network models through a lot of flexibility offered by the way it utilises meta-programming for designing models - fortunately, it is also really quite easy to use. The following is a brief introduction.

## Installation
Spike is not currently available as a released package, but may be installed using the Julia package manager. From your Julia REPL, type `]` to enter the Pkg mode and run
```Julia
  pkg> add https://github.com/FabulousFabs/Spike.jl.git
```
That's it!

## Magic models
Generally speaking, all Spike models boil down to three fundamental components: [`Model`](@ref)s that contain and track all specifications, [`NeuronGroup`](@ref)s and [`Synapses`](@ref) where the former contains a collection of neurons governed by the same equations, whereas the latter is much the same but for synapses between [`NeuronGroup`](@ref)s.

To make things a bit easier, we generally don't have to specify [`Model`](@ref) structures ourselves. Instead, we can [`cast_magic()`](@ref) and let Spike figure out the rest behind the scenes for us. Once we have cast magic, we can specify structures and connections in our main file.

## Adding neurons
For starters, let us create some neurons, then:
```Julia
using Spike;

# Casting magic ensures that the back-end will put together a model structure for us.
cast_magic();

# Next, let's actually specify some neurons
population = NeuronGroup(N = 10, method = euler, events = Dict(),
                         eq = :(dv_dt = (.-(v - E_L) .+ I) ./ 𝜏_v;
                                I_t = 3e-3 .* rand(N) ./ dt;),
                         parameters = Dict(:v => zeros(10),
                                           :v_th => ones(10),
                                           :v_reset => -0.1 * ones(10),
                                           :E_L => zeros(10),
                                           :I => zeros(10),
                                           :𝜏_v => 100e-3 * ones(10)));
```

Now, the most pertinent part here is the specification of a set of equations that will govern the behaviour of these neurons. We specified these equations as an expression and, as we see through the use of ```rand(N)```, we can also use regular Julia functions in them if need be. Similarly, however, we also need to obey to Julia's way of writing maths. For example, all parameters (e.g., `v`, `E_L`, ...) are vectors and, as such, we need to use dot notation `.-`, `.+`, `.*`, `./` in our equations, too. More on expressions and how to phrase your equations later.

## Running a magic model
For now, we have 10 neurons in our group who receive random input current. Let's simulate this population and see what some of the currents and voltages look like. To do this, we will use a [`StateMonitor`](@ref) and [`run`](@ref) where, given that we do not supply a [`Model`](@ref) structure, [`run`](@ref) will build a magic network from our scope for us. Like so:

```Julia
# Let's add a monitor that measures v and I
state_monitor = StateMonitor(obj = population, vars = [:v, :I]);

# Let's do a magic run
run(; T = 1.0, dt = 1e-3);

# ...and plot the results for our neurons
using Plots;
plot(repeat(state_monitor.t, outer = [1, 10]), transpose(state_monitor.states[:I]))
plot(repeat(state_monitor.t, outer = [1, 10]), transpose(state_monitor.states[:v]))
```

Nice, that was pretty easy! Okay, so the input is all over the place, as we expected. But wait, what's going on with the membrane potential? Why isn't it being reset but keeps on rising?

## Adding events
Well, it turns out that we specified `events = Dict()` for our neurons. We can create arbitrary events that will govern some conditional behaviour of our neurons. Note that, importantly, events are also always broadcast to synapses. More on this later. For the time being, let us go back to our neuron group and add some events to our call, like so:

```Julia
# let's change the definition of our neuron group
population = NeuronGroup(N = 10, method = euler,
                         eq = :(dv_dt = (.-(v - E_L) .+ I) ./ 𝜏_v;
                                I_t = 3e-3 .* rand(N) ./ dt;),
                         events = Dict(:spike => (:(v .> v_th), :(v = v_reset;))),
                         parameters = Dict(:v => zeros(10),
                                           :v_th => ones(10),
                                           :v_reset => -0.1 * ones(10),
                                           :E_L => zeros(10),
                                           :I => zeros(10),
                                           :𝜏_v => 100e-3 * ones(10)));

...

# let's also use an event monitor to track spikes directly
spike_monitor = EventMonitor(obj = population, event = :spike);

...

# and let's also plot them
scatter(spike_monitor.t, spike_monitor.i, markershape = :vline)
```

That's more like it, we have quite some spiking behaviour now. Very nice! 
