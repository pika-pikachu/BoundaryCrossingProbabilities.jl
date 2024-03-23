@doc """
	transitionMatrix(k, p, gU, gL)

Returns the taboo transition matrix of the diffusion from time t_{k-1} to t_k.

""" -> 
function transitionMatrix(k::Int64, p::MeshParams, gU::Function, gL::Function)
s = t_n(k-1, p, gU, gL)
t = t_n(k, p, gU, gL)
grid_s = gridLattice(k-1, p, gU, gL) # grid at time s
grid_t = gridLattice(k, p, gU, gL) # grid at time t
dh = grid_t[1] - grid_t[2] # spacing of grid. Using implied spacing to avoid numerical instability
M = zeros(Complex{Float64},length(grid_s), length(grid_t)) # Transition matrix of taboo transition probabilities
	for i in 1:length(grid_s)
		x = grid_s[i]
		for j in 1:length(grid_t)
			y = grid_t[j]
			M[i, j] = dh * transitionDensity(s, x, t, y, p) *
						 bridgeCorrection(s, x, t, y, gU, gL, p) *
                        potentialCorrection(s, x, t, y, p)
		end
		if p.absorb == true # only for one dimensional
			#M[i, length(grid_t)] = bridgeCorrection(s, x, t, grid_t[end], gU, gL, p) * absorb_prob(s, x, t, grid_t[end], dh, p) 
			M[i, length(grid_t)] = absorb_prob(s, x, t, grid_t[end], dh, p) 
		end
	end
	if p.absorb == true && k > 1
		M[length(grid_s), :] .= 0.0
		M[length(grid_s), length(grid_t)] = 1.0 # Once at the boundary the Markov chain will remain there
	end
return M 
end;



@doc """
	FKE(p, gU, gL, plotOnFlag, plot3DFlag, matrixOutputFlag)

Returns approximation of u(t,x) is the solution of the forward Kolmogorov equation PDE:

∂_tu = (L^*)u 
u(t,g+(t)) = 0
u(t,g-(t)) = 0
u(0,x) = δ(x-x0) for x ∈ (g-(t),g+(t))

""" -> 
function FKE(p::MeshParams, gU::Function, gL::Function, plot3DFlag = false)
#try 
    P = zeros(length(gridLattice(0,p,gU,gL)))' # initial distribution
    P[end] = 1 # Dirac mass at zero
    pointCount = 1 # count number of grid points in the current time-slice
    u_xyz = Array{Complex{Float64}}(undef, 0, 3) # initialising the PDE solution surface
        for k in 1:p.n
            grid_s = gridLattice(k-1, p, gU, gL) 
            grid_t = gridLattice(k, p, gU, gL)
            λ = [grid_s, # used in doubleExponential integration
                gL(t_n(k,p,gU,gL)), 
                gU(t_n(k,p,gU,gL))] # collection of points grid in original space
            w = weightIntegration(
                length(grid_t), 
                p.integration_method, 
                λ, 
                h(k, p, gU, gL)) # weight vector
            P = (P * transitionMatrix(k,p,gU,gL)) .* w' 
            if plot3DFlag == true
                dh = grid_t[1] - grid_t[2]
                if p.absorb == true
                    pointCount = vcat(pointCount, length(grid_t) - 1) # drop the lower boundary point
                    u_xyz = solutionGenerator(u_xyz, grid_t[1:(end-1)], k, P[1:(end-1)]/dh, gU, gL, p)
                else
                    pointCount = vcat(pointCount, length(grid_t)) # drop the lower boundary point
                    u_xyz = solutionGenerator(u_xyz, grid_t, k, P/dh, gU, gL, p)
                end
            end
        end
    nonCrossingProbability = sum(P)
        if plot3DFlag == true
            plotEngine(u_xyz, pointCount, "forward")
        end
        if p.target_set_bool == true
            return sum(P[2:(end-1)]) + (P[1] + P[end])/2
        else
            return P
        end
#catch err
#    return 0 
#end
end;


@doc """
	BKE(p, gU, gL, plotOnFlag, plot3DFlag, matrixOutputFlag)

Returns the solution of the backward Kolmogorov equation PDE at v(0,0)

(∂_t + L)v = 0 
v(t,U(t)) = 0
v(t,L(t)) = 0
v(T,x) = 1 for x ∈ (L(t),U(t))

""" -> 
function BKE(p::MeshParams, gU::Function, gL::Function, plot3DFlag = false, returnInterpolation = false, method = "backward")
P = ones(length(gridLattice(p.n, p, gU, gL)))'
pointCount = 1 # initialise number of points for a given t_k
u_xyz = Array{Float64}(undef, 0, 3)
	for k in range(p.n, 1, step = -1)
		# P = P*transitionMatrix(k,p,gU,gL)'
		grid_s = gridLattice(k-1, p, gU, gL) # grid at time t
		grid_t = gridLattice(k, p, gU, gL) # grid at time t
		P = P .* weightIntegration(length(grid_t), p.integration_method)' * transitionMatrix(k,p,gU,gL)'
        if plot3DFlag == true 
            if (p.absorb == true) & (k > 1)
                pointCount = vcat(pointCount, length(grid_s) - 1) # drop the lower boundary point
                u_xyz = solutionGenerator(u_xyz, grid_s[1:(end-1)], k - 1, P[1:(end-1)], gU, gL)
            else
                pointCount = vcat(pointCount, length(grid_s)) # drop the lower boundary point
                u_xyz = solutionGenerator(u_xyz, grid_s, k - 1, P, gU, gL, p)
            end
		end
	end
    if plot3DFlag == true
        plotEngine(u_xyz, pointCount, method)
    end
	if returnInterpolation == true
		x,y,z = real(u_xyz)[:,1], real(u_xyz)[:,2], real(u_xyz)[:,3] 
		spl = Spline2D(x,y,z; kx=3, ky=3, s=1e-4)
		v(s,x) = evalgrid(spl,[s],[x])[1]	
		return u_xyz, v
	else
		return u_xyz
	end
end;

@doc """
	solutionGenerator(u_xyz, pointCount)

Returns a typle containing:
	1. a vector of (x, y, z) coordinates for the solution of the partial differential equation.
	2. number of points corresponding to the kth time slice

""" -> 
function solutionGenerator(u_xyz, space_mesh, k, P, gU, gL, p::MeshParams, scale = false)
	n_length = length(P')
	if (n_length == 1) & (k == 0) 
		x0 = 0
		u_xyz = vcat(u_xyz, hcat(0, x0, P[:]))
		return u_xyz
	end
	x = t_n(k, p, gU, gL) * ones(n_length)
	y = space_mesh
	z = P'[:]
	u_xyz = vcat(u_xyz, hcat(x, y, z))
	return u_xyz
end

@doc """
	scatter3d(x, y, z, cs, ax, colorsMap)

Initiates a 3D scatter plot for 3 vectors, x, y and z. A colour map is generated based on the z variable.
Default perceptually uniform "plasma" colourmap.

Reference: https://stackoverflow.com/questions/8891994/matplotlib-3d-scatter-plot-with-color-gradient

""" -> 
function scatter3d(x, y, z, cs, ax, colorsMap= "plasma")
    cm = plt.get_cmap(colorsMap)
    cNorm = matplotlib.colors.Normalize(vmin=minimum(cs), vmax=maximum(cs))
    scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=cm)
    ax.scatter(x, y, z, c = scalarMap.to_rgba(cs), s = 3, alpha = 1)
    scalarMap.set_array(cs)
    #colorbar(scalarMap)
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
end

@doc """
	line3d(u_xyz, pointCount, alpha)

Plots 3D plot based on splits of vector of x,y,z coordinates. The splits are determined by the pointCount.
alpha controls the opacity of the 3D lines.
data2 = bcpBKE(p, gU,t->-3, true, true, true)
u_xyz = data2[1]
pointCount = data2[2]
line3D(u_xyz, pointCount)
""" -> 
function line3D(u_xyz, pointCount, polarity = "real", alpha = 0.6)
	cumsum_point_count = cumsum(pointCount[1:(end-1)])
	u_xyz = u_xyz[1:(end-1),:,:]
	for i in 1:(length(cumsum_point_count)-1)
        sub_u_xyz = u_xyz[cumsum_point_count[i]:(cumsum_point_count[i+1]-1),:]
        if polarity == "real"
            plot3D(real(sub_u_xyz[:,1]), real(sub_u_xyz[:,2]), real(sub_u_xyz[:,3]), color = "black", alpha = alpha, marker = " ",linewidth = 1)
        elseif polarity == "imag"
            plot3D(real(sub_u_xyz[:,1]), real(sub_u_xyz[:,2]), imag(sub_u_xyz[:,3]), color = "black", alpha = alpha, marker = " ",linewidth = 1)
        end
    end
end;


@doc """
	plotEngine(u_xyz, pointCount, p, value, title1, title2, zaxislabel, style = "color")

title1: Main title of the plot
title2: Output name
value: Numerical value corresponding to title2
zaxislabel: String label for z-axis

Helper plotting engine. 
""" -> 
function plotEngine(u_xyz, pointCount, direction = "forward")
	#fig = plt.figure(figsize = (4.5/2, 4.5/2), dpi = 100)
	#fig = plt.figure(figsize = (6, 6), dpi = 300)    
    if (direction == "forward") | (direction == "backward")
        fig = plt.figure(figsize=(10, 5), dpi = 300)
        ax_1 = subplot(1,2,1,projection="3d")
        scatter3d(real(u_xyz[:,1]), real(u_xyz[:,2]), real(u_xyz[:,3]), real(u_xyz[:,3]), ax_1)
        line3D(u_xyz, pointCount, "real")
        if direction == "forward"
            ax_1.view_init(40., 20.)
            ax_1.set_zlabel("Re u(t,x)")          
        else # backward equation
            ax_1.view_init(40., 200.)
            ax_1.set_zlabel("Re v(t,x)")                  
        end
        ax_1.set_xlabel("t (Time)")
        ax_1.set_ylabel("x (Space)")    
        ax_1.set_title("Real part")
        ax_2 = subplot(1,2,2,projection="3d")
        scatter3d(real(u_xyz[:,1]), real(u_xyz[:,2]), imag(u_xyz[:,3]), imag(u_xyz[:,3]), ax_2)
        line3D(u_xyz, pointCount, "imag")    
        if direction == "forward"
            ax_2.view_init(40., 20.) 
            ax_2.set_zlabel("Im u(t,x)") 
        else
            ax_2.view_init(40., 200.)    
            ax_2.set_zlabel("Im v(t,x)") 
        end
        ax_2.set_xlabel("t (Time)")
        ax_2.set_ylabel("x (Space)")    
        ax_2.set_title("Imaginary part")  
    else
        fig = plt.figure(figsize=(10, 5), dpi = 300)
        ax_1 = subplot(1,1,1,projection="3d")
        scatter3d(real(u_xyz[:,1]), real(u_xyz[:,2]), real(u_xyz[:,3]), real(u_xyz[:,3]), ax_1)
        line3D(u_xyz, pointCount, "real")
        if direction == "forward_real"
            ax_1.view_init(40., 20.)
            ax_1.set_zlabel("u(t,x)")          
        else # backward equation
            ax_1.view_init(40., 200.)
            ax_1.set_zlabel("v(t,x)")                  
        end
        ax_1.set_xlabel("t (Time)")
        ax_1.set_ylabel("x (Space)")
    end
end;