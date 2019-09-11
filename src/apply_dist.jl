# try to emulate the apply_link_new.jl supertype so I can bring in mean and scale

# This super type of all response distribution types

abstract type DistributionFunction end

##distributions##

"""Binomial Distribution"""

struct BinomialDist <: DistributionFunction
end

function binomial_dist(n,s,x)
  return Binomial(n,x)
end

"""Negative Binomial Distribution"""

struct NegativeBinomialDist <: DistributionFunction 
end

function negative_binomial_dist(n,r,x)
  return NegativeBinomial(r,x)
end 

"""Poisson Distribution"""

struct PoissonDist <: DistributionFunction
end

function poisson_dist(n,s,x)
  return Poisson(x)
end

"""Exponential Distribution"""
struct ExponentialDist <: DistributionFunction
end

function exponential_dist(n,s,x)
  return Exponential(x)
end

"""Gamma Distribution"""

struct GammaDist <: DistributionFunction
end

function gamma_dist(n,s,x)
   y = x/s
  return Gamma(s,y)
end

"""Inverse Gaussian Distribution."""

struct InverseGaussianDist <: DistributionFunction
end

function inverse_gaussian_dist(n,s,x)
     lamda = (1.0/s)^2
  return InverseGaussian(x,lamda)
end

"""Logistic Distribution"""
 
struct LogisticDist <: DistributionFunction
end


function logistic_dist(n,s,x)
  return Logistic(x,s)
end

"""Lognormal Distribution"""

struct LognormalDist <: DistributionFunction
end

function lognormal_dist(n,s,x)
  return LogNormal(x,s)
end

"""Normal Distribution"""

struct NormalDist <: DistributionFunction
end

function normal_dist(n,s,x)
  return Normal(x,s)
end


""" Weibull Distribution"""
struct Weibull <: DistributionFunction
end

function weibull(n,s,x)
   return Weibull(x,s)
end

##APPLY DISTRIBUTIONS

apply_dist(n,s,x,dist::BinomialDist) = binomial_dist.(n,s,x)

apply_dist(n,s,x,dist::NegativeBinomialDist) = negative_binomial_dist.(n,s,x)

apply_dist(n,s,x,dist::PoissonDist) = poisson_dist.(n,s,x)

apply_dist(n,s,x,dist::ExponentialDist) = exponential_dist.(n,s,x)

apply_dist(n,s,x,dist::GammaDist) = gamma_dist.(n,s,x)

apply_dist(n,s,x,dist::InverseGaussianDist) = inverse_gaussian_dist.(n,s,x)

apply_dist(n,s,x,dist::LogisticDist) = logistic_dist.(n,s,x)

apply_dist(n,s,x,dist::LognormalDist) = lognormal_dist.(n,s,x)

apply_dist(n,s,x,dist::NormalDist) = normal_dist.(n,s,x)

apply_dist(n,s,x,dist::Weibull) = weibull_dist.(n,s,x)
