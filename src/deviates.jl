# module RandomDeviates

# export bernoulli_deviate, binomial_deviate, exponential_deviate
# export gamma_deviate, inverse_gaussian_deviate, normal_deviate
# export poisson_deviate, t_deviate, weibull_deviate

using SpecialFunctions, Missings

"""Generates a Bernoulli random deviate with success probability p."""

function bernoulli_deviate(p::T) where T <: Real
  if rand(T) <= p
    return 1
  else
    return 0
  end
end

"""Generates a binomial random deviate with success probability p
and n trials."""

function binomial_deviate(p::T, n::Int) where T <: Real
  sucesses = 0
  for i = 1:n
    if rand(T) <= p
      sucesses = sucesses + 1
    end
  end
  return sucesses
end

"""Generates an exponential random deviate with mean mu."""

function exponential_deviate(mu::T) where T <: Real
  return -mu * log(rand(T))
end

"""Generates a gamma deviate with shape parameter alpha
and intensity lambda."""

function gamma_deviate(alpha::T, lambda::T) where T <: Real
  n = floor(Int, alpha)
  z = - log(prod(rand(T, n)))
  beta = alpha - n
  if beta <= one(T) / 10^6
    y = zero(T)
  elseif beta < one(T) / 10^3
    y = (beta / rand(T))^(one(T) / (one(T) - beta))
  else
    (r, s) = (one(T) / beta, beta - one(T))
    for i = 1:1000
      u = rand(T, 2)
      y = - log(one(T) - u[1]^r)
      if u[2] <= (y / (one(T) - exp(- y)))^s
        exit
      end
    end
  end
  return (z + y) / lambda
end

"""Generates an inverse Gaussian deviate with mean mu and
scale lambda."""

function inverse_gaussian_deviate(lambda::T, mu::T) where T <: Real
  w = mu * randn()^2
  c = mu / (2 * lambda)
  x = mu + c * (w - sqrt(w * (4 * lambda + w)))
  p = mu / (mu + x)
  if rand() < p
    return x
  else
    return mu^2 / x
  end
end 

"""Generates a normal random deviate with mean mu and 
standard deviation sigma."""

function normal_deviate(mu::T, sigma::T) where T <: Real
  return sigma * randn(T) + mu
end

"""Generates a Poisson random deviate with mean mu."""

function poisson_deviate(mu::T) where T <: Real
  (x, p, k) = (exp(-mu), one(T), 0)
  while p > x
    k = k + 1
    p = p * rand(T)
  end
  return k - 1
end

"""Generates a t random deviate with, location mu, scale sigma,
and degrees of freedom df."""

function t_deviate(mu::T, sigma::T, df::T) where T <: Real
  x = randn(T)
  y = 2 * gamma_deviate(df / 2, one(T))
  return sigma * (x / sqrt(y / df)) + mu
end

"""Generates a Weibull random deviate with scale lambda
and shape alpha."""

function weibull_deviate(lambda::T, alpha::T) where T <: Real
  return lambda * (- log(rand()))^(one(T) / alpha)
end

# end # module RandomDeviates


# module InverseLinks

# export cauchit_inverse_link, cloglog_inverse_link, identity_inverse_link
# export inverse_inverse_link, logit_inverse_link, log_inverse_link
# export probit_inverse_link, sqrt_inverse_link

"""inverse cauchit link."""

function cauchit_inverse_link(x::t) where t <: Real
  return atan(x) / pi + one(t) / 2
end

"""inverse cloglog link."""

function cloglog_inverse_link(x::t) where t <: Real
  return one(t) - exp(-exp(x))
end 

"""inverse identity link."""

function identity_inverse_link(x::t) where t <: Real
  return x
end

"""inverse inverse link."""

function inverse_inverse_link(x::t) where t <: Real
  return one(t) / x
end

"""inverse logit link."""

function logit_inverse_link(x::t) where t <: Real
  return one(t) / (one(t) + exp(-x))
end

"""inverse log link."""

function log_inverse_link(x::t) where t <: Real
  return exp(x)
end

"""inverse probit link."""

function probit_inverse_link(x::t) where t <: Real
  return (one(t) + erf(x / sqrt(2 * one(t)))) / 2
end

"""inverse sqrt link."""

function sqrt_inverse_link(x::t) where t <: Real
  return x * x
end

# end # module InverseLinks

# module SimulateData

# export ResponseType
# export missing_entries, simulate_random_effects, simulate_glm_trait

mutable struct ResponseType
  family :: String
  inverse_link :: String
  location :: Float64
  scale :: Float64
  shape :: Float64
  df :: Float64 # degrees of freedom
  trials :: Int
end

"""Apply inverse link."""

function apply_inverse_link(μ::Vector{T}, dist::ResponseType) where T <: Real
  if dist.inverse_link == "Cauchit"
    μ = cauchit_inverse_link.(μ)
  elseif dist.inverse_link == "Cloglog"
    μ = cloglog_inverse_link.(μ)
  elseif dist.inverse_link == "Identity"
    μ = identity_inverse_link.(μ)
  elseif dist.inverse_link == "Inverse"
    μ = inverse_inverse_link.(μ)
  elseif dist.inverse_link == "Logit"
    μ = logit_inverse_link.(μ)
  elseif dist.inverse_link == "Log"
    μ = log_inverse_link.(μ)
  elseif dist.inverse_link == "Probit"
    μ = probit_inverse_link.(μ)
  elseif dist.inverse_link == "Sqrt"
    μ = sqrt_inverse_link.(μ)
  else
    error("Link function not supported!")
  end
end

"""Simulate GLM model traits."""

function simulate_glm_trait(μ::Vector{T}, dist::ResponseType) where T <: Real
  if dist.family == "bernoulli"
    x = bernoulli_deviate.(μ)
  elseif dist.family == "binomial"
    x = binomial_deviate.(μ, dist.trials)
    if dist.trials < 0
      error("Trials cannot be negative for a binomial distribution!")
    end
  elseif dist.family == "exponential"
    x = exponential_deviate.(μ)
  elseif dist.family == "gamma"
    if dist.shape < zero(T)
      error("Shape cannot be negative for a gamma distribution!")
    end
    x = gamma_deviate.(dist.shape, μ)
  elseif dist.family == "inverse gaussian"
    if dist.scale <= zero(T)
      error("Scale must be positive for an inverse Gaussian distribution!")
    end
    x = inverse_gaussian_deviate.(dist.scale, μ)
  elseif dist.family == "normal"
    if dist.scale < zero(T)
      error("Scale cannot be negative for a normal distribution!")
    end
    x = normal_deviate.(μ, dist.scale)
  elseif dist.family == "poisson"
    x = poisson_deviate.(μ)
  elseif dist.family == "t"
    if dist.scale < zero(T)
      error("Scale cannot be negative for a t distribution!")
    end
    x = t_deviate.(μ, dist.scale, dist.df)
  elseif dist.family == "weibull"
    if dist.shape <= zero(T)
      error("Shape must be positive for a Weibull distribution!")
    end
    x = weibull_deviate.(μ, dist.shape)
  else
    error("Response distribution not supported!")
  end
  return x
end

#"""Simulate random effects."""

# function simulate_random_effects(vc::Vector{VarianceComponent}, 
#   people::Int, traits::Int)
# 
#   X = zeros(people, traits)
#   for i = 1:length(vc) # components
#     if isa(vc[i].var_comp, Number) || isa(vc[i].var_comp, Vector)
#       A = diagm(vc[i].var_comp)
#     else # isa(vc[i].var_comp, Matrix)
#       A = vc[i].var_comp
#     end
#     B = vc[i].cov_mat
#     Y = randn(people, traits) # standard normal deviates
#     X = X + chol(B)' * (Y * chol(A)) # Roth's Kronecker lemma
#   end
#   return X
# end

using DataFrames

"""Recreate missing pattern of original columns in target columns 
and create additional random missing entries."""

function missing_entries!(df::DataFrame, original::Vector{Int}, 
  target::Vector{Int}, missing_rate::T) where T <: Real
  (rows, m, n) = (size(df, 1), length(original), length(target))
  @assert m == n || m == 0 "In trait simulation the number of original 
    traits should match the number of simulated traits."
  if m == n
    for j = 1:n
      for i = 1:rows
        if ismissing(df[i, original[j]]) || rand(T) <= missing_rate
          df[i, target[j]] = missing
        end
      end
    end
  else
    for j = 1:n
      for i = 1:rows
        if rand(T) <= missing_rate
          df[i, target[j]] = missing
        end
      end
    end
  end
end

# end # module SimulateData

using StatsBase
using SpecialFunctions
using DataFrames
using Missings
# using RandomDeviates
# using InverseLinks
# using SimulateData
(trials, mu, lambda, p, n, df) = (1000000, 2.0, 2.0, 0.5, 10, 5.0)
 x = zeros(trials);
 for i = 1:trials
   x[i] = bernoulli_deviate(p)
   x[i] = binomial_deviate(p, n)
   x[i] = exponential_deviate(mu)
   x[i] = poisson_deviate(mu)
   x[i] = bernoulli_deviate(p)
   #x[i] = t_deviate(df)
  x[i] = inverse_gaussian_deviate(lambda, mu)
   x[i] = weibull_deviate(lambda, 2.0)
 end
# describe(x)
# var(x)
# println(lambda*sqrt(pi)/2)
# x = 2 * randn(10);
# sort!(x)
# x = CauchitLink.(x)
# x = CloglogLink.(x)
# x = IdentityLink.(x)
# x = InverseLink.(x)
# x = LogitLink.(x)
# x = LogLink.(x)
# x = ProbitLink.(x)
# x = SqrtLink.(x)
# df = DataFrame(A=1:4, B=rand(4), C=randstring.([3,3,3,3]), D=rand(4))
# for j = 1:size(df,2)
#   df[:, j] = allowmissing(df[:, j])
# end
# original = [2] # zeros(Int, 0)
# target = [4]
# missing_rate = 0.5
# df[2,2] = missing
# missing_entries!(df, original, target, missing_rate)
# 
μ = rand(5);


dist = ResponseType(" ", "", 0.0, 0.0, 0.0, 0.0, 0)
# dist.family = "bernoulli"
# dist.family = "binomial"
# dist.family = "exponential"
# dist.family = "gamma"
# dist.family = "inverse gaussian"
# dist.family = "normal"
# dist.family = "poisson"
# dist.family = "t"
dist.family = "weibull"
dist.scale = 1.0
dist.shape = 1.0
dist.df = 3.0
dist.trials = 2
x = simulate_glm_trait(μ, dist)
println(μ)
println(x)