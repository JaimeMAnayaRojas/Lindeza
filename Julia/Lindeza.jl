#Import Turing, Distributions and DataFrames
using Turing, Distributions, DataFrames, Distributed

# Import MCMCChain, Plots, and StatsPlots for visualizations and diagnostics.
using MCMCChains, Plots, StatsPlots
using CSV, CategoricalArrays
# Set a seed for reproducibility.
using Random
Random.seed!(12);

# Turn off progress monitor.
Turing.setprogress!(false)
# Turn off progress monitor.
pwd()

cd("/home/jaime/Dropbox/Projects_JM/Muenster/Lindeza/") # set the working directory
cd("/home/mao/Dropbox/Projects_JM/Muenster/Lindeza/") # set the working directory
readdir() # check the files in the data/ folder



data = CSV.read("data/Larvae_adults.csv", DataFrame) # read the data

data.Rep = string.(data.Regime, "-",  data.Replicate)
data.Rep =replace(data.Rep, Pair.(unique(data.Rep), axes(unique(data.Rep), 1))...)


df = select(data, 
    :larvae => :Y,
    :Regime => (x-> recode(x, "Control" => 0, "Matched" => 1, "Unmatched" => 0 )) => :M,
    :Regime => (x-> recode(x, "Control" => 0, "Matched" => 0, "Unmatched" => 1 )) => :U,
    :Rep => :R,
    :Generation 

)


@model AR(Y, G,  M, U, R) = begin
    
    G1 = G .- 4
    G2 = G.^2 .- 4^2
    G3 = G.^3 .- 4^3
    # Set variance prior.
   # σ ~ truncated(Normal(0, 1), 0, Inf)
    σᵣ ~ truncated(Normal(0, 1), 0, Inf)
    αᵣ ~ filldist(Normal(0, σᵣ), length(unique(R)))
    α ~ Normal(4,3) 
    βₘ ~ Normal(0,1)
    βᵤ ~ Normal(0,1)
    β_G1 ~ Normal(0,1)
    β_G2 ~ Normal(0,1)
    β_G3 ~ Normal(0,1)
    β_GM ~ Normal(0,1)
    β_GU ~ Normal(0,1)

    λ = α .+ β_G1 .* G1 .+ β_G2 .* G2 .+ β_G3 .* G3 .+ βₘ .* M .+ βᵤ .* U .+ # Effects on the intercept
        G1 .* (β_GM .* M  + β_GU .* U) .+ # Effects on the slope 
         αᵣ[R] # Random effects on the intercept (i.e., replicate)

    Y .~ Poisson.(exp.(λ)) 

end    

# Sample using NUTS(n_iters::Int, n_adapts::Int, δ::Float64), where:
# n_iters::Int : The number of samples to pull.
# n_adapts::Int : The number of samples to use with adapatation.
# δ::Float64 : Target acceptance rate.



model = AR(df.Y,df.Generation, df.M, df.U, df.R)
num_chains = 4

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