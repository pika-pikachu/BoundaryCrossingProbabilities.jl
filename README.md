# BoundaryCrossingProbabilities.jl

This is a julia package for computing accurate approximations of the time-dependent boundary crossing probability for a general diffusion process. The implemented algorithm is based on a Brownian bridge corrected Markov chain algorithm proposed in Liang & Borovkov (2023).

## Background

For $T>0,$ $x_0 \in \mathbb{R},$ let $X$ be a solution to the following stochastic differential equation:

$$ \begin{cases}
dX_t = \mu(t,X_t)\,dt +\sigma(t,X_t)\,dW_t, \quad t \in (0,T),\\ 
X_0 = x_0.
\end{cases}$$

where we assume $\mu$ and $\sigma$ satisfy the usual sufficient for the existence and uniqueness of a strong solution to the above SDE. 

For two continuous functions $g_-$ and $g_+,$ this package computes an approximation of the non-crossing probability

$$ 
	\mathbf{P}(g_-(t) < X_t < g_+(t) , t\in [0,T]).
$$

Moreover, the tool can actually compute expressions of the form

$$ v(t,x) = \mathbf{E}[e^{-\int_t^TV(s,X_s)\,ds}\psi(X_T);g_{-}(s) < X_s < g_+(s), s \in [t,T] | X_t = x], $$

and the (discounted) taboo transition density

$$ u(t,x) = \frac{\partial}{\partial x}\mathbf{E}[e^{-\int_0^tV(s,X_s)\,ds}\mathbf{1}( X_t \leq x, g_{-}(s) < X_s < g_+(s), s \in [0,t])|X_0 = x_0], $$

which are known to be probabilistic solutions to the following parabolic PDEs:

$$ Lv = 0, \quad v(t,g_{\pm}(t)) =0, \quad v(T,x)= \psi(x). $$

$$ L^*u = 0, \quad u(t,g_{\pm}(t)) =0, \quad u(0,x)= \delta_{x_0}(x) $$

where 

$$L f(s,x) := \dot{f}(s,x) + \mu(s,x) f'(s,x) + \frac{1}{2}\sigma^2(s,x)f''(s,x).$$


## Code

To install and import the BoundaryCrossingProbabilities package, run the following at the Julia REPL

```
Pkg.add("BoundaryCrossingProbabilities")
using BoundaryCrossingProbabilities
```

To set up the boundary crossing probability algorithm, we need to specify parameters of the diffusion process (initial position, drift, diffusion and potential), and the time interval. Let's take the Brownian motion example with a complex potential term.

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
gU(t) = 4 - t^2
gL(t) = -4 + t^2
```

Now we can obtain the solution to the problem. 

```
plotFlag = true

soltn_BKE, v = BoundaryCrossingProbabilities.BKE(p, t -> gU(t), t -> gL(t), true, plotFlag);
soltn_FKE, u = BoundaryCrossingProbabilities.FKE(p, t -> gU(t), t -> gL(t), true, plotFlag);
```

![Screenshot](complex_potential.png)

