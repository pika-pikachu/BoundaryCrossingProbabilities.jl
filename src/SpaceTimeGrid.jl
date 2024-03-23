@doc """
	t_n(k,p, gU, gL)

returns the time slice at time k
""" -> 
function t_n(k::Int64, p::MeshParams, gU::Function, gL::Function)
	return k/p.n*p.T # uniform slicing
	# Non uniform slicing
	# if k == 0
	# 	return 0
	# else
	# 	a = Int(floor(p.n*2/10))
	# 	b = Int(p.n - a)
	# 	scaling = vcat(2*ones(a),ones(b))
	# 	dtVector = scaling/sum(scaling)
	# 	tVector = cumsum(dtVector)
	# 	return tVector[k]
	# end
end

@doc """
	h(k, p, gU, gL)

Returns the grid span at time k. The last grid span is finer due to the Brownian bridge correction
only occuring on the left hand side (in time). This reduces the decay rate the boundary at the last point.
""" -> 
function h(k::Int64, p::MeshParams, gU::Function, gL::Function)
t = t_n(k, p, gU, gL)
dt = t - t_n(k-1, p, gU, gL)
p_vec = (p.δ + 1/2) * ones(p.n) # power for 1 to (n-1)th 
p_vec[end] = p.p2 # power for last transition
if k == p.n && (p.target_set == true)
	a = p.target_set[1]
	b = p.target_set[2]
	beta_k = b - a # distance between boundaries
	x = beta_k/dt^p_vec[k]
	γ_k = x/floor(p.γ * 2 * x)
	return γ_k*dt^p_vec[k]
elseif k == 0  
	σ0 = minimum([p.σ(t,x) for x in range(gU(t),gL(t),length = 20)])
	return gU(0)/floor(p.γ * gU(0)/σ0/dt^p_vec[1])
else
	# return (gU(t)[1] - gL(t)[1])/floor(Int, p.γ * (gU(t)[1] - gL(t)[1])/dt^p_vec[k])
	# return (gU(t) - gL(t))/floor(p.γ * (gU(t) - gL(t))/p.σ(0,0)/dt^p_vec[k])
	σ0 = minimum([p.σ(t,x) for x in range(gU(t),gL(t),length = 20)])
	# γ should be set around 1.5-2, otherwise the exponentially decaying error becomes noticeable
	return (gU(t) - gL(t))/floor(p.γ * (gU(t) - gL(t))/σ0/dt^p_vec[k])
end
end


@doc """
	gridLattice(k, p, gU, gL)

Returns equally spaced points between the functions gU and gL at time point point t_k with spacing parameters p.

""" -> 
function gridLattice(k::Int64, p::MeshParams, gU::Function, gL::Function)
if k == 0 
	#return p.x0 # Initial point
	grid = range(gU(0), stop = p.x0, length = round(Int64,(gU(0)-p.x0)/h(0, p, gU, gL)) + 1 ) # trapezoidal rule
	return grid
end	
t_k = t_n(k, p, gU, gL) 
h_k = h(k, p, gU, gL) # grid spacing at time t_k
	if p.absorb == true # one sided boundary
		# grid = (gU(t_k)[i] - h_k/2):(-h_k):(gL(t_k)[i]) # floating midpoint rule
		#grid = gU(t_k):(-h_k):gL(t_k) # floating trapezdoial grid
		grid = range(gU(t_k), stop = gL(t_k), length = round(Int64, (gU(t_k)- gL(t_k))/h_k)) # trapezoidal grid
		# grid = range(gU(t_k)[i] - h_k/2, stop = gL(t_k)[i] + h_k/2, length = floor(Int64,(gU(t_k)[i] - gL(t_k)[i])/h_k)) # midpoint rule
	else # two sided boundary
		if k == p.n && p.target_set_bool == true # TARGET SET
			a = p.target_set[1]
			b = p.target_set[2]
			grid = range(b, stop = a, length = round(Int64,2*(b-a)/h_k)+1) 
		else
			grid = range(gU(t_k), stop = gL(t_k), length = round(Int64, (gU(t_k) - gL(t_k))/h_k) + 1 ) # trapezoidal rule
			# grid = range(gU(t_k) - h_k/2, stop = gL(t_k) + h_k/2, length = round(Int64,(gU(t_k)- gL(t_k))/h_k)) # midpoint rule
		end
	end
return grid
end;
