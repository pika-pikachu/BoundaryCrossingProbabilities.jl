# BoundaryCrossingProbability.jl
Computes the boundary crossing probability for a general diffusion process and time-dependent boundary.

To install the BoundaryCrossingProbabilities Julia package, open Julia,  press "]" and type 

```
add BoundaryCrossingProbabilities
```

To import the package, type 

```
using BoundaryCrossingProbabilities
```

To set up the algorithm, we need to specify parameters of the diffusion process (initial position, drift, diffusion and potential), and the time interval. Let's take the Brownian motion example with a complex potential term.

```
x0 = 0 # Initial condition
b(t,x) = 0 # Drift coefficient
σ(t,x) = 1 # Diffusion coefficient
V(t,x) = (1im)*x^2  # Potential
T = 1.0 # Terminal Time
```

Then we set up a Julia Type called MeshParams, which takes in all the diffusion process parameters.

```
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
```

Then we define the upper and lower boundaries

```
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
```

Now we can obtain the solution to the problem. 

```
plotFlag = true

soltn_BKE, v = BoundaryCrossingProbabilities.BKE(p, t -> gU(t), t -> gL(t), true, plotFlag);
soltn_FKE, u = BoundaryCrossingProbabilities.FKE(p, t -> gU(t), t -> gL(t), true, plotFlag);
```


