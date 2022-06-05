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
For now, we have 10 neurons in our group who receive random input current. Let's simulate this population and see what some of the currents and voltages look like. To do this, we will use a [`StateMonitor`](@ref) and run where, given that we do not supply a [`Model`](@ref) structure, run will build a magic network from our scope for us. Like so:

```Julia
...

# Let's add a monitor that measures v and I
state_monitor = Spike.StateMonitor(obj = population, vars = [:v, :I]);

# Let's do a magic run
Spike.run(; T = 1.0, dt = 1e-3);

# ...and plot the results for our neurons
using Plots;
plot(repeat(state_monitor.t, outer = [1, 10]), transpose(state_monitor.states[:I]))
plot(repeat(state_monitor.t, outer = [1, 10]), transpose(state_monitor.states[:v]))
end
```

Nice, that was pretty easy! Okay, so the input is all over the place, as we expected. But wait, what's going on with the membrane potential? Why isn't it being reset but keeps on rising?

## Adding events
Well, it turns out that we specified `events = Dict()` for our neurons. We can create arbitrary events that will govern some conditional behaviour of our neurons. Note that, importantly, events are also always broadcast to synapses. More on this later. For the time being, let us go back to our neuron group and add some events to our call, like so:

```Julia
...

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

## Adding synapses
Now, let's move on and create a slightly more complex model where groups of neurons actually communicate with each other. To do this, let us change our population such that its input currents are no longer random but simply decay in time. This means adjusting

```Julia
population = NeuronGroup(N = 100, method = euler,
                         eq = :(dv_dt = (.-(v - E_L) .+ I) ./ 𝜏_v;
                                dI_dt = .-I ./ 𝜏_I;),
                         events = Dict(:spike => (:(v .> v_th), :(v = v_reset;))),
                         parameters = Dict(:v => zeros(100),
                                           :v_th => ones(100),
                                           :v_reset => -0.1 * ones(100),
                                           :E_L => zeros(100),
                                           :I => zeros(100),
                                           :𝜏_v => 50e-3 * ones(100),
                                           :𝜏_I => 10e-3 * ones(100)));
pop_states = StateMonitor(obj = population, vars = [:v, :I]);
pop_spikes = EventMonitor(obj = population, event = :spike);
```

Excellent. Notice that we also increased the population size and adjusted some parameters to make it more interesting. Next, let's create a second population that we will use as a random spiking population, like so:

```Julia
rd_spikers = NeuronGroup(N = 100, method = euler,
                         eq = :(v_t = rand(N);),
                         events = Dict(:spike => (:(v .> p), :(v = v_reset;))),
                         parameters = Dict(:v => zeros(100),
                                           :v_reset => zeros(100),
                                           :p => 0.975 * ones(100)));
rds_spikes = EventMonitor(obj = rd_spikers, event = :spike);
```

And now, let's create some [`Synapses`](@ref) like so:

```Julia
forward = Synapses(pre = rd_spikers, post = population,
                   cond = :((i .- 25) .<= j .&& (i .+ 25) .>= j),
                   prob = 0.25,
                   on_pre = Dict(:spike => :(post_I = post_I .+ w;)),
                   parameters = Dict(:w => rand));
```

Now, there are a few important things to take note of here:
1. We are using an expression, again, as a condition for creating synapses that, effectively, amounts to a kernel of size 50.
2. On top of that, we are asking that synapses be created with a probability of .25 only, effectively yielding 12-13 synapses per neuron.
3. We are now using the event `:spike` that we created earlier and are using `on_pre` to create a conditional equation. In this equation, we have full access to all parameters of the synapse, the pre- and the postsynaptic neuron by using no prefix, `pre_*` or `post_*`, respectively. Note that `on_post` is also available.
4. We are specifying parameters for the synapse, but we are using a new syntax that is not vectorised that applies only to synapses. Parameters can either be of type `Number` or `Function`.
5. If we wanted to, we could also specify `eq` and `method` for these synapses, just like we did for neurons.

For completeness sake, here is the full example so far:

```Julia
using Spike;
using Plots;

# Casting magic ensures that the back-end will put together a model structure for us.
cast_magic();

# Let's create some random spikers
rd_spikers = NeuronGroup(N = 100, method = euler,
                         eq = :(v_t = rand(N);),
                         events = Dict(:spike => (:(v .> p), :(v = v_reset;))),
                         parameters = Dict(:v => zeros(100),
                                           :v_reset => zeros(100),
                                           :p => 0.975 * ones(100)));

# Let's create our target population
population = NeuronGroup(N = 100, method = euler,
                         eq = :(dv_dt = (.-(v - E_L) .+ I) ./ 𝜏_v;
                                dI_dt = .-I ./ 𝜏_I;),
                         events = Dict(:spike => (:(v .> v_th), :(v = v_reset;))),
                         parameters = Dict(:v => zeros(100),
                                           :v_th => ones(100),
                                           :v_reset => -0.1 * ones(100),
                                           :E_L => zeros(100),
                                           :I => zeros(100),
                                           :𝜏_v => 50e-3 * ones(100),
                                           :𝜏_I => 10e-3 * ones(100)));

# And let's add some synapses
forward = Synapses(pre = rd_spikers, post = population,
                   cond = :((i .- 25) .<= j .&& (i .+ 25) .>= j),
                   prob = 0.25,
                   on_pre = Dict(:spike => :(post_I = post_I .+ w;)),
                   parameters = Dict(:w => rand));

# Let's add a monitor that measures v and I
pop_states = StateMonitor(obj = population, vars = [:v, :I]);
pop_spikes = EventMonitor(obj = population, event = :spike);
rds_spikes = EventMonitor(obj = rd_spikers, event = :spike);

# Let's do a magic run
Spike.run(; T = 1.0, dt = 1e-3);

# ...and plot the results for our neurons
plot(repeat(pop_states.t, outer = [1, 10]), transpose(pop_states.states[:I]))
plot(repeat(pop_states.t, outer = [1, 10]), transpose(pop_states.states[:v]))
scatter(rds_spikes.t, rds_spikes.i, markershape = :vline)
scatter(pop_spikes.t, pop_spikes.i, markershape = :vline)
```

Before we move on, feel free to toy with this example a little bit to familiarise yourself with these basic functions.

## On expressions
There are a couple more rules to follow when writing expressions that have not been mentioned yet and should now be mentioned (perhaps you have even discovered some of them while playing around with the model above).

- Expressions can have one of four forms:
  1. differential equations: `dx_dt = .- x ./ 𝜏;`
  2. linear equations: `x_t = rand(N) .* y;`
  3. assignment: `x = x_r;`
  4. conditional: `x .== y`
- Subexpressions should always be terminated with a semicolon.
- Currently, only `eq::Expr` parameters (i.e., in [`NeuronGroup`](@ref) and [`Synapses`](@ref)) support differential equations. Hooks, like `on_pre`, for example, do not.

With that in mind, let's move to the final section of this introduction.

## Adding operations
Finally, let's bring it all together. Say we want to have some kind of time-varying input to our model that we simulated otherwise (for example, we are simulating an agent's decision in an environment). To do this, we will need to be able to run our model in parallel to our simulation. We can do this by using an [`Operation`](@ref). Take the following model, for example:

```Julia
using Spike;
using Plots;

# Casting magic ensures that the back-end will put together a model structure for us. # hide
cast_magic();

# Let's create our target population
N_pop = 100
population = NeuronGroup(N = N_pop, method = euler,
                         eq = :(dv_dt = (.-(v - E_L) .+ I) ./ 𝜏_v;
                                dI_dt = .-I ./ 𝜏_I;),
                         events = Dict(:spike => (:(v .> v_th), :(v = v_reset;))),
                         parameters = Dict(:v => zeros(N_pop),
                                           :v_th => ones(N_pop),
                                           :v_reset => -0.1 * ones(N_pop),
                                           :E_L => zeros(N_pop),
                                           :I => zeros(N_pop),
                                           :𝜏_v => 50e-3 * ones(N_pop),
                                           :𝜏_I => 10e-3 * ones(N_pop)));

# Let's capture spikes
pop_spikes = EventMonitor(obj = population, event = :spike);

# An arbitrary function
function my_function()
       # here's where I can do something else, at the same clock speed as the model
       # and then feed something back into the model.
       # for example:
       global population;
       population.parameters[:I] .= 5 .* rand(population.N);
end

# register the operation
callback = Operation(op = my_function);

# Run everything
Spike.run(; T = 1.0, dt = 1e-3);

# ...and plot the results for our neurons
scatter(pop_spikes.t, pop_spikes.i, markershape = :vline)
```

Now, we can easily simulate the model and, for example, an agent in an environment whose actions the model determines at the same time.

## Advanced usage
For more information, refer to the [Examples](@ref) or [Library](@ref).
