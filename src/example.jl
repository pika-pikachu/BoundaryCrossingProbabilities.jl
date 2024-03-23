include("./BoundaryCrossingProbabilities.jl")

import .BoundaryCrossingProbabilities
using PyPlot

b(t,x) = 0 # Drift coefficient
σ(t,x) = 1 # Diffusion coefficient
V(t,x) = (1im)*x^2  # Potential
T = 1.0 # Terminal Time
x0 = 0 # Initial condition

p = BoundaryCrossingProbabilities.MeshParams(
	25, # n
	0, # δ, 1/2 + δ is the space step power before the final time
	1, # pn power of the space step at the final time
	T, # T # Terminal time
	x0, # x0 Initial condition
	2, #γ,
	"Brownian", #bridge correction,
	false, #one sided boundary
	b, # Drift coefficient
	σ, # Diffusion coefficient
    V, # Potential
	false, # no target set
	[1.2,3], # Target set
	"trapezoidal" # integration scheme
	);

function gU(t, a=2, theta = 2)
    if t <= 0.00634616586979
        return theta/2
    else
        return t/theta*acosh(a*exp(theta^2/(2t)))
    end
end;

function gL(t, a=2, theta = 2)
    if t <= 0.00634616586979
        return -theta/2
    else
        return -t/theta*acosh(a*exp(theta^2/(2t)))
    end
end;

soltn_BKE, v = BoundaryCrossingProbabilities.BKE(p, t-> 2, t-> -2, true, true);
soltn_FKE = BoundaryCrossingProbabilities.FKE(p, t-> 2, t-> -2, true);

v(0,0)
