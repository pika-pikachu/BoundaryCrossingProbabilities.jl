@doc """
	transitionDensity(s, x, t, y, p)

(s,x): starting position
(t,y): ending position

Returns Edgeworth expansion of the transition density of the diffusion process.
Reference: Sahalia, A (2001) 
""" -> 
function transitionDensity(s, x, t, y, p::MeshParams) 
	Δ = t - s
	δt, δx = set_variables("δt δx", order = 3)
	b_taylor_x = p.b(s, x + Taylor1(Float64, 4)) # x partial derivative
	b_taylor_t = p.b(s + Taylor1(Float64, 2), x) # t partial derivative
	b_taylor_tx = p.b(s + δt, x + δx) # x,t mixed partial derivative
	σ_taylor_x = p.σ(s, x + Taylor1(Float64, 4)) # x partial derivative
	σ_taylor_t = p.σ(s + Taylor1(Float64, 2), x) # t partial derivative
	σ_taylor_tx = p.σ(s + δt, x + δx) # x, t mixed partial derivative 
	if typeof(b_taylor_t) != Taylor1{Float64} # Drift not time-dependent
		b_t = b_tt = b_tx = b_txx = 0
	else
		b_t = b_taylor_t[1]
		b_tt = b_taylor_t[2]*2
		b_tx = b_taylor_tx[2][2]
		b_txx = b_taylor_tx[3][3]*2
	end
	if typeof(σ_taylor_t) != Taylor1{Float64} # Diffusion not time-dependent
		σ_t = σ_tt = σ_tx = σ_txx = 0
	else
		σ_t = σ_taylor_t[1]
		σ_tt = σ_taylor_t[2]*2
		σ_tx = σ_taylor_tx[2][2]
		σ_txx = σ_taylor_tx[3][3]*2
	end
	if typeof(b_taylor_x) != Taylor1{Float64} # Drift Not space-dependent
		b_x = b_xx = b_xxx = b_xxxx = 0
	else
		b0 = b_taylor_x[0]
		b_x = b_taylor_x[1]
		b_xx = b_taylor_x[2]*2
		b_xxx = b_taylor_x[3]*3*2
		b_xxxx = b_taylor_x[4]*4*3*2
	end
	if typeof(σ_taylor_x) != Taylor1{Float64} # Non-space dependent diffusion Gaussian case
		σ_x = σ_xx = σ_xxx = σ_xxxx = 0
		if typeof(b_taylor_x) != Taylor1{Float64} # Brownian motion case
			return exp(-(y-x)^2/2/Δ)/sqrt(2*pi*Δ)
		else
		b_taylor = p.b(s, x + Taylor1(Float64, 4))
		σ0 = p.σ(s, x)
		b0 = b_taylor[0]
		b_x = b_taylor[1]
		b_xx = b_taylor[2]*2
		b_xxx = b_taylor[3]*3*2
		b_xxxx = b_taylor[4]*4*3*2
		μ1 = Δ*b0 + Δ^2/2*(b_x * b0 + b_xx/2*σ0^2) + Δ^3/24*(4*b0^2*b_xx + 4*b0*(b_x^2 + σ0*(σ0*b_xxx)) + σ0^2*(6*b_x*b_xx + σ0^2*b_xxxx))
		μ2 = Δ*σ0^2 + Δ^2*(b0^2 + 0.5*σ0^2*(2*b_x)) +
					  Δ^3/12*(4*b0^2*(3*b_x) + 
					  			14*b0*σ0^2*b_xx + 
								σ0^2*( 8*b_x^2 + σ0^2*4*b_xxx))
		σ2 = μ2 - μ1^2
		z = (y - x - μ1)/sqrt(σ2)
		end
		return exp(-z^2/2)/sqrt(2 * pi * σ2)
	else # Edgeworth expansion
		σ0 = σ_taylor_x[0]
		σ_x = σ_taylor_x[1]
		σ_xx = σ_taylor_x[2]*2
		σ_xxx = σ_taylor_x[3]*3*2
		σ_xxxx = σ_taylor_x[4]*4*3*2
		μ1 = Δ*b0 + Δ^2/2*(b_x * b0 + b_xx/2*σ0^2 + b_t) + 
				Δ^3/24*(4*b0^2*b_xx + 
							4*b0*(b_x^2 + σ0*(σ_x * b_xx + σ0*b_xxx) + 2*b_tx) + 
							σ0^2*(6*b_x*b_xx + 2*σ_x^2*b_xx + 2*σ0*b_xx*σ_xx + 
									4*σ0*σ_x *b_xxx + σ0^2*b_xxxx +4*b_txx) +
							4*σ0*b_xx*σ_t 
						)
		μ2 = Δ*σ0^2 + Δ^2*(b0^2 + σ0*(b0*σ_x + σ_t) + 1/2*σ0^2*(2*b_x + σ_x^2 + σ0*σ_xx)) +
				  Δ^3/12*(4*b0^2*(3*b_x + σ_x^2 + σ0*σ_xx) + 
				  			2*b0*σ0*(6*b_x*σ_x + 2*σ_x^3 + 8*σ0*σ_x*σ_xx + σ0*(7*b_xx + 2*σ0*σ_xxx)) + 
							σ0^2*( 8*b_x^2 + 2*σ_x^4 + 16*σ0*σ_x^2*σ_xx + 
									8*b_x*(σ_x^2 + σ0*σ_xx) + 
									2*σ0*σ_x*(7*b_xx + 4*σ0*σ_xxx) + 
									σ0^2*(5*σ_xx^2 + 4*b_xxx + σ0*σ_xxxx) +
									8*σ_xx*σ_t + 8*b_tx + 8*σ_x*σ_tx
									) +
							2*b0*(6*b_t + 4*σ_t*σ_x + 4*σ0*σ_tx) +
							4*σ0^3*σ_txx +
							4*σ0*(σ_x*b_t + 2*b_x*σ_t + σ_x^2*σ_t + σ_tt) +
							4*σ_t^2
							)
		σ2 = μ2 - μ1^2
		if σ2 <= 0
			return 0
		end
		μ3 = (3*Δ^2*σ0^3*σ_x + 
				Δ^3*σ0^2*(4*b0*σ_x^2 + 
					σ0*(5*b_x*σ_x + 2*(2*σ_x^3 + b0*σ_xx + σ_tx)) + 
					σ0^2*(b_xx + 7*σ_x*σ_xx) + 
					4*σ_x*σ_t - 3*b_t +
					σ0^3*σ_xxx))/sqrt(σ2)^3
		μ4 = (3*Δ^2*σ0^4 + 
				Δ^3*σ0^3*(6*(b0*σ_x + σ_t) + σ0*(6*b_x + 19*σ_x^2) + 7*σ0^2*σ_xx))/σ2^2
		z = (y - x - μ1)/sqrt(σ2)
		c1 = 1/6*(z^3 - 3*z)*μ3
		c2 = 1/24*(z^4 - 6*z^2 + 3)*(μ4 - 3) + 1/72*(μ3)^2*(-15 + 45*z^2 - 15*z^4 + z^6)
		# c2 = 1/24*(z^4 - 6*z^2 + 3)*(μ4 - 3) 
		return exp(-z^2/2)/sqrt(2 * pi * σ2)*(1 + c1 + c2)
	end
end;


function bridgeCrossing(s, x, t, y, g::Function, p::MeshParams)
	if p.correction == "LargeDeviation"
		Δ = t - s
		inv_σ(x) = 2/(p.σ(s,x) + p.σ(t,x))
		Fx = (g(s) - x)/2*sum(gaussWeights .* (inv_σ.((g(s)-x)/2*gaussNodes .+ (g(s) + x)/2)))
		Fy = (g(s) - y)/2*sum(gaussWeights .* (inv_σ.((g(s)-y)/2*gaussNodes .+ (g(s) + y)/2)))
		# Fx = numIntegrate(inv_σ, x, g(s), gaussWeights, gaussNodes)
		# Fy = numIntegrate(inv_σ, y, g(s), gaussWeights, gaussNodes)
		# dg = ε₁part(g(Hyper(s, 1.0, 0.0, 0.0)))
		dg = (g(t) - g(s))/(t-s)
		Cf = exp(-2 * inv_σ(g(s)) * dg * Fx)
		return Cf*exp(-2/Δ*Fx*Fy)
	elseif p.correction == "Brownian"
		return exp(-2/(t-s)/((p.σ(s,x) + p.σ(t,x))/2)^2*(g(s) - x)*(g(t) - y))
	else
		return 0
	end
end

@doc """
	bridgeCorrection(x, y, s, t, gU, gL, p)

(s,x): starting position
(t,y): ending position

Returns the probability of not hitting either the upper boundary gU or gL 
""" -> 
function bridgeCorrection(s, x, t, y, gU::Function, gL::Function, p::MeshParams)
	if p.correction != "NoCorrection" # Bridge correction
		if p.absorb == true # one sided boundary
			death_prob_vec = bridgeCrossing(s, x, t, y, gU, p)
		else
			death_prob_vec = bridgeCrossing(s, x, t, y, gU, p) +
									bridgeCrossing(s, x, t, y, gL, p)
		end	
		# return max(0,min(1,1 - death_prob_vec)) # makes almost no difference for remote boundaries.
		return 1 - death_prob_vec
	else
		if (gL(s) < x < gU(s)) && (gL(t) < y < gU(t))
			return 1
		else 
			return 0
		end
	end
end;

function potential(s,t,x,y,p::MeshParams) 
    # Trapezoidal rule; equivalent to averaging
    return -(p.V(s,x) + p.V(t,y))/2
end;

function potentialCorrection(s,x,t,y,p::MeshParams) 
    # Trapezoidal rule; equivalent to averaging
    trapezoidal_V = (p.V(s,x) + p.V(t,y))/2
    return exp(-trapezoidal_V*(t-s))
end;