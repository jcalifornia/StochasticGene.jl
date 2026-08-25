# simulator Function

Simulate stochastic gene expression models using the Gillespie algorithm.

## Syntax

```julia
simulator(rin, transitions, G, R, S, insertstep; kwargs...)
```

## Arguments

### Required Arguments

- `rin::Vector{Float64}`: Initial transition rates
- `transitions::Tuple`: Tuple of vectors specifying state transitions
- `G::Int`: Number of gene states
- `R::Int`: Number of pre-RNA steps
- `S::Int`: Number of splice sites (must be ≤ R - insertstep + 1)
- `insertstep::Int`: Step where reporter becomes visible

### Optional Keyword Arguments

#### Simulation Parameters
- `warmupsteps::Int = 0`: Number of warmup steps before recording
- `nalleles::Int = 1`: Number of alleles
- `nhist::Int = 20`: Number of histogram bins
- `bins::Vector{Float64} = Float64[]`: Custom histogram bins
- `traceinterval::Float64 = 0.0`: Time interval for trace recording (0 = no traces)

#### Model Configuration
- `coupling::Tuple = tuple()`: Coupling parameters for multi-unit models
- `onstates::Vector{Int} = Int[]`: States where transcription is active
- `splicetype::String = ""`: Splicing configuration ("", "offeject", "offdecay")

#### Observation Model
- `probfn::Function = prob_Gaussian`: Probability function for observations
- `noise::Vector{Float64} = Float64[]`: Noise parameters
- `reportersteps::Vector{Int} = Int[]`: Steps where reporter is visible

#### Output Control
- `tspan::Tuple{Float64, Float64} = (0., 1000.)`: Time span for simulation
- `ntrials::Int = 1`: Number of simulation trials
- `resultfolder::String = ""`: Output folder for results

## Returns

- `histogram::Vector{Float64}`: Steady-state RNA count distribution
- `traces::Vector{Vector{Float64}}`: Intensity traces (if `traceinterval > 0`)
- `dtimes::Vector{Float64}`: Dwell times (if applicable)

## Examples

### Basic Two-State Model

```julia
using StochasticGene

# Simple two-state telegraph model
rates = [0.1, 0.2]  # G1->G2, G2->G1
transitions = ([1,2], [2,1])
G, R, S = 2, 0, 0
insertstep = 1

# Simulate RNA histogram
histogram = simulator(
    rates, transitions, G, R, S, insertstep,
    nhist = 50,
    ntrials = 1000
)
```

### GRS Model with Traces

```julia
# Gene-Reporter-Splice model
rates = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
transitions = ([1,2], [2,1])
G, R, S = 2, 3, 2
insertstep = 1

# Simulate with intensity traces
histogram, traces = simulator(
    rates, transitions, G, R, S, insertstep,
    traceinterval = 1.0,    # 1 minute intervals
    tspan = (0., 500.),     # 500 minute simulation
    noise = [10.0, 5.0],    # Gaussian noise parameters
    probfn = prob_Gaussian
)
```

### Coupled Model

```julia
# Two model definitions and two simulation units.
transitions = (([1, 2], [2, 1]), ([1, 2], [2, 1]))
G, R, S, insertstep = (2, 2), (1, 1), (0, 0), (1, 1)

# Unit 1 in state 1 modifies transition 1 of unit 2.
coupling = ((1, 2), [(1, 1, 2, 1)])

# Complete rates for model 1, then model 2, then one gamma per connection.
rates = [
    0.30, 0.20, 0.40, 0.50, 1.0,
    0.25, 0.35, 0.45, 0.55, 1.0,
    0.20,
]

residence = simulator_ss(
    rates, transitions, G, R, S, insertstep;
    coupling=coupling,
    noiseparams=[0, 0],
    warmuptime=100.0,
    totalsteps=200_000,
)
```

The canonical coupling form is `(unit_model, connections[, sign_modes])`.
Each connection `(beta, s, alpha, t)` means that source unit `beta` in state
`s` affects transition `t` of target unit `alpha`. The first tuple maps units
to model definitions; repeated mappings such as `(1, 1)` are supported.

For reporter-state coupling:

```julia
# Rsum: every R position is a separate connection.
rsum = make_coupling("R5", (3, 3), (3, 3))

# Rany: one sentinel source state, active if any R position is occupied.
rany = ((1, 2), [(1, G[1] + R[1] + 1, 2, 5)])
```

With `Rsum`, `n` occupied positions and a tied strength produce the factor
`1 + n*gamma`. With `Rany`, any nonzero occupancy produces `1 + gamma` once.
Direct simulator vectors contain one gamma per connection, so repeat a tied
`Rsum` gamma across all expanded R-position connections. A tied `Rsum` over
`m` simultaneously occupiable positions requires `gamma > -1/m`; `Rany`
requires `gamma > -1`.

### Hierarchical Model

```julia
# Simulate data for hierarchical fitting
rates = [0.1, 0.2, 0.3]
transitions = ([1,2], [2,1])
G, R, S = 2, 1, 0
insertstep = 1

# Generate multiple traces for hierarchical analysis
ntraces = 10
all_traces = Vector{Vector{Float64}}()

for i in 1:ntraces
    # Add noise to rates for each trace
    noisy_rates = rates .* (1 .+ 0.1 * randn(length(rates)))
    
    histogram, traces = simulator(
        noisy_rates, transitions, G, R, S, insertstep,
        traceinterval = 0.5,
        tspan = (0., 200.),
        noise = [20.0, 10.0]
    )
    
    append!(all_traces, traces)
end
```

### Custom Observation Model

```julia
# Define custom observation function
function prob_Poisson(y, μ, σ)
    return pdf(Poisson(μ), round(Int, y))
end

# Simulate with Poisson observation noise
histogram = simulator(
    rates, transitions, G, R, S, insertstep,
    probfn = prob_Poisson,
    noise = [15.0],  # Poisson rate parameter
    nhist = 30
)
```

## Rate Ordering

The rate vector `rin` must follow this specific ordering:

1. **G transitions**: Rates between gene states
2. **R transitions**: Rates between pre-RNA steps
3. **S transitions**: Splicing rates
4. **Decay rates**: mRNA decay rates
5. **Noise parameters**: Observation noise parameters

For coupled models, concatenate each **model definition's complete block** in
model order. A block contains that model's kinetic rates followed by its noise
parameters. Append one coupling strength per connection. This is model order,
not unit order: if `unit_model == (1, 1)`, model 1's block appears once and is
reused by both units.

The simulator validates the exact number of model, noise, and coupling values.
It also validates the unit-model map and optional coupling sign modes, so a
malformed coupled vector fails before simulation rather than shifting rate
indices silently.

### Example Rate Ordering

For a model with G=2, R=3, S=2:
```julia
rates = [
    # G transitions
    0.1,    # G1 -> G2
    0.2,    # G2 -> G1
    
    # R transitions
    0.3,    # R1 -> R2
    0.4,    # R2 -> R3
    0.5,    # R3 -> eject
    
    # S transitions
    0.6,    # Splice site 1
    0.7,    # Splice site 2
    
    # Decay
    0.05    # mRNA decay
]
```

## Performance Notes

1. **Memory Usage**: Large `nhist` values require more memory
2. **Simulation Time**: Longer `tspan` and smaller rates increase runtime
3. **Trace Recording**: `traceinterval > 0` significantly increases memory usage
4. **Parallel Processing**: Use multiple calls for embarrassingly parallel simulations

## Common Use Cases

### Parameter Estimation Validation
```julia
# Generate synthetic data for testing parameter estimation
true_rates = [0.1, 0.2, 0.3]
synthetic_data = simulator(
    true_rates, transitions, G, R, S, insertstep,
    nhist = 100,
    ntrials = 5000
)

# Use synthetic_data to validate fitting algorithms
```

### Model Comparison
```julia
# Compare different model structures
models = [
    (2, 0, 0),  # Simple telegraph
    (2, 1, 0),  # With pre-RNA
    (3, 1, 0),  # Three gene states
]

for (G, R, S) in models
    histogram = simulator(rates, transitions, G, R, S, insertstep)
    # Analyze and compare histograms
end
```

### Burst Analysis
```julia
# Simulate for burst size analysis
histogram, traces = simulator(
    rates, transitions, G, R, S, insertstep,
    traceinterval = 0.1,    # High temporal resolution
    tspan = (0., 1000.),
    onstates = [2]          # State 2 is transcriptionally active
)

# Analyze burst properties from traces
```

## Error Handling

The function includes validation for:
- Rate vector length consistency
- Valid state transitions
- Proper model dimensions (G, R, S relationships)
- Positive rate values
- Valid time spans

## See Also

- `simulate_trace`: Generate traces only
- [`simulate_trials`](analysis.md#Validate-theory-against-simulation): Compare
  theoretical and simulated correlations, including threaded repeated
  experiments and ON-ON validation
- `fit`: Fit models to data
- `prob_Gaussian`: Gaussian observation model
