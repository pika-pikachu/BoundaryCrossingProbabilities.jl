include("./BoundaryCrossingProbabilities.jl")

import .BoundaryCrossingProbabilities
using PyPlot

b(t,x) = 0 # Drift coefficient
σ(t,x) = 1 # Diffusion coefficient
V(t,x) = (1im)*x^2  # Potential
T = 1.0 # Terminal Time
x0 = 0 # Initial condition

p = BoundaryCrossingProbabilities.MeshParams(
    x0, # x0 Initial condition
    T, # Terminal time	
    b, # Drift coefficient
	σ, # Diffusion coefficient
    V, # Potential
    false, # no target set
	[1.2,3], # Target set X_T \in [a,b]
    "Brownian", #bridge correction,
	false, # one sided boundary
    25, # number of time steps 
	0, # δ, 1/2 + δ is the space step power before the final time
	1, # pn power of the space step at the final time
	2, # γ, constant space scaling
	"trapezoidal" # integration scheme
	);

function gU(t, a = 2, theta = 2)
    if t <= 0.00634616586979
        return theta/2
    else
        return t/theta*acosh(a*exp(theta^2/(2t)))
    end
end;

function gL(t, a = 2, theta = 2)
    if t <= 0.00634616586979
        return -theta/2
    else
        return -t/theta*acosh(a*exp(theta^2/(2t)))
    end
end;

plotFlag = true

soltn_BKE, v = BoundaryCrossingProbabilities.BKE(p, t -> gU(t), t -> gL(t), true, plotFlag);
soltn_FKE, u = BoundaryCrossingProbabilities.FKE(p, t -> gU(t), t -> gL(t), true, plotFlag);


soltn_BKE
soltn_FKE

v(0,0)

u(1,2)
