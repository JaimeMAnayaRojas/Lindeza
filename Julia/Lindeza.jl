#Import Turing, Distributions and DataFrames
using Turing, Distributions, DataFrames, Distributed, StatisticalRethinking

# Import MCMCChain, Plots, and StatsPlots for visualizations and diagnostics.
using MCMCChains, Plots, StatsPlots
using CSV, CategoricalArrays
# Set a seed for reproducibility.
using Random 
using StatsModels


Random.seed!(12);

# Turn off progress monitor.
Turing.setprogress!(false)
# Turn off progress monitor.
pwd()

cd("Dropbox/Jaime M/Projects_JM/Muenster/Lindeza/Lindeza") # set the working directory

readdir() # check the files in the data/ folder



data = CSV.read("data/Larvae_adults.csv", DataFrame) # read the data

Gmat  = string.("G", data.Generation)


# data.Rep =replace(data.Rep, Pair.(unique(data.Rep), axes(unique(data.Rep), 1))...)
df = select(data, 
    :larvae,
    :Regime => (x-> recode(x, "Control" => 0, "Matched" => 1, "Unmatched" => 0 )) => :Matched,
    :Regime => (x-> recode(x, "Control" => 0, "Matched" => 0, "Unmatched" => 1 )) => :Unmatched,
    :Replicate
)



df.G4 = Int64.(data.Generation .== 4)
df.G5 = Int64.(data.Generation .== 5)
df.G6 = Int64.(data.Generation .== 6)
df.G7 = Int64.(data.Generation .== 7)
df.G8 = Int64.(data.Generation .== 8)
df.G9 = Int64.(data.Generation .== 9)
log(mean(df.larvae))


@model Larvae_mod(larvae, G4, G5, G6, G7, G8, G9,  Match, Unmatch, Replicate) = begin
    
    Y = larvae
    M = Match
    U = Unmatch
    R = Replicate




   # Set variance prior.
   # σ ~ truncated(Normal(0, 1), 0, Inf)
    σᵣ ~ truncated(Normal(0, 2), 0, Inf)
    αᵣ ~ filldist(Normal(0, σᵣ), length(unique(R)))
    α ~ Normal(4,8) 
    βm ~ Normal(0,1)
    βu ~ Normal(0,1)
    β4 ~ Normal(0,1)
    β5 ~ Normal(0,1)
    β6 ~ Normal(0,1)
    β7 ~ Normal(0,1)
    β8 ~ Normal(0,1)
    β9 ~ Normal(0,1)


    β4m ~ Normal(0,1)
    β5m ~ Normal(0,1)
    β6m ~ Normal(0,1)
    β7m ~ Normal(0,1)
    β8m ~ Normal(0,1)
    β9m ~ Normal(0,1)


    β4u ~ Normal(0,1)
    β5u ~ Normal(0,1)
    β6u ~ Normal(0,1)
    β7u ~ Normal(0,1)
    β8u ~ Normal(0,1)
    β9u ~ Normal(0,1)

    λ = α .+  β4 .* G4  .+ β5 .* G5  .+ β6 .* G6  .+ β7 .* G7  .+ β8 .* G8  .+ β9 .* G9  .+ # Effects on the slope
             (β4m .* G4  .+ β5m .* G5  .+ β6m .* G6  .+ β7m .* G7  .+ β8m .* G8  .+ β9m .* G9 .+ βm) .* M  .+
             (β4u .* G4  .+ β5u .* G5  .+ β6u .* G6  .+ β7u .* G7  .+ β8u .* G8  .+ β9u .* G9 .+ βu) .* U  .+
             αᵣ[R] # Random effects on the intercept (i.e., replicate)

    Y .~ Poisson.(exp.(λ)) 

end    

# Sample using NUTS(n_iters::Int, n_adapts::Int, δ::Float64), where:
# n_iters::Int : The number of samples to pull.
# n_adapts::Int : The number of samples to use with adapatation.
# δ::Float64 : Target acceptance rate.



model = Larvae_mod(df.larvae, df.G4, df.G5, df.G6, df.G7, df.G8, df.G9, df.Matched, df.Unmatched, df.Replicate)
num_chains = 2

chain = sample(model, NUTS(2000, 0.95), MCMCThreads(), 4_000, num_chains; discard_adapt=true);



# Removing the first 200 values of the chains.
chains_new = chain[1001:4000,:,:]
(chains_new)

plot(chains_new)

post = DataFrame(chains_new)

println(names(post)[2:27])

post2 = post[:,3:27]


ci = (mapcols(x -> quantile(x, [0.025, 0.975]),post2))
sum_par =DataFrame(parameter = names(post2), mean = round.(mean.(eachcol(post2)), digits = 3), 
                    lc = round.(Vector(ci[1, :]), digits = 3), up = round.(Vector(ci[2, :]), digits = 3))


println(sum_par)