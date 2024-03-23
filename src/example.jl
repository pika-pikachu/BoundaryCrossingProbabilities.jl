include("./BoundaryCrossingProbabilities.jl")

import .BoundaryCrossingProbabilities
using PyPlot

b(t,x) = 0 # Drift coefficient
σ(t,x) = 1 # Diffusion coefficient
V(t,x) = (1im)*x^2  # Potential
T = 1.0 # Terminal Time
x0 = 0 # Initial condition

p = BoundaryCrossingProbabilities.MeshParams(
	70, # n
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

y_mesh = range(-4,stop = 4,length = 1000)
x0 = 1
s = 0
t = 1
figure()
plt.plot(y_mesh,[BoundaryCrossingProbabilities.transitionDensity(s, x0, t, y, p)  for y in y_mesh]);
plt.grid();
plt.xlabel("y")
plt.ylabel(L"p(s,x_0,t,y)");
plt.title("Transition density");
display(gcf())

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

soltn_BKE = BoundaryCrossingProbabilities.BKE(p, t-> 2, t-> -2, true);

display(gcf())

using Interpolations

real(soltn_BKE)[:,1]
real(soltn_BKE)[:,2]
real(soltn_BKE)[:,3]

real(soltn_BKE)[:,1:2]

tuple(real(soltn_BKE)[:,1:2])

tuple(real(soltn_BKE)[:,1],real(soltn_BKE)[:,2])
(x,y)
tuple(real(soltn_BKE)[:,1],real(soltn_BKE)[:,2])
real(soltn_BKE)[:,3]

itp = LinearInterpolation(tuple(real(soltn_BKE)[:,1],real(soltn_BKE)[:,2]),real(soltn_BKE)[:,3])

itp = LinearInterpolation(real(soltn_BKE)[:,1],real(soltn_BKE)[:,2],
real(soltn_BKE)[:,3])

x = range(-2, 3, length = 10)
y = range(-2, 3, length = 10)
z = @. cos(x) + sin(y')
z
itp = LinearInterpolation((x,y),z)
(x,y)
itp(0.2,0.5)