module BoundaryCrossingProbabilities

greet() = print("Thank you for installing BoundaryCrossingProbabilities")

using PyPlot # Plotting package
using HyperDualNumbers # Used for automatic differentiation
using CPUTime # Timing of functions
using Distributions # used for path simulations
using FastGaussQuadrature # used to get Gaussian integrals
using SpecialFunctions # used to get sine integral
using TaylorSeries # higher order derivatives
using Dierckx # 2D Interpolation for unstructured grids

struct MeshParams
    x0::Float64  # Initial point
    T::Float64 # Time horizon
    b::Function # Drift coeffcient
    σ::Function # Diffusion coefficient
    V::Function # Potential
    target_set_bool::Bool # target set
    target_set::Array{Float64,1} # target set vector of form [a,b], a<b
    correction::String # Bridge correction type
    absorb::Bool # One sided boundary absorption
	n::Int64 # Number of nodes in time partition
    δ::Float64 # 1/2 + δ power of first to second last point
    p2::Float64 # power of last point
    γ::Float64 # extra scaling
    integration_method::String # Method of integration {"trapezoidal", "simpsons", "doubleExponential"}
end

include("./TransitionDensities.jl")
include("./IntegrationWeights.jl")
include("./SpaceTimeGrid.jl")
include("./MarkovChainApproximation.jl")

export 
    MeshParams
    transitionDensity
    FKE
    BKE
end