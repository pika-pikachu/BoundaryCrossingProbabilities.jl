@doc """

weightIntegration(m, method, x)

m: number of nodes
method: integration weight method
x: vector of initial points in original space [gL,gU]

Returns quadrature weights.

""" ->
function weightIntegration(m, method = "trapezoidal", x =[1,1,1], h_k = 1)
	if (method == "trapezoidal") | (method == "trapezoidalAdjusted")
		return ones(m)
	elseif method == "simpsons"
		if m == 1 
			return 1
		elseif m == 2 #Trapezoidal rule# 
			return [1,1]/2 
		elseif m == 3 #Simpsons 1/3 rule# 
			return [1,4,1]/3
		elseif m == 4 #Simpsons 3/8 rule#
			return [3,9,9,3]/8
		elseif isodd(m) # simpson 1/3 rule
			return vcat(1,repeat([4,2], Int((m-1)/2))[1:end-1], 1)/3
		else # composite simpson 1/3, 3/8 rule
			return vcat(1/3,repeat([4,2]/3, Int(m/2-2))[1:end-1],1/3 + 3/8, 9/8,9/8, 3/8)
		end
	elseif method == "boole"
		return repeat([32,12,32,14],3)
	elseif method == "doubleExponential"
		if m == 1
			return 1
		end
		L = 2.8 # any higher and the inversion loses signifinace
		Z = h_k*(ceil(L/h_k):(-1):floor(-(L/h_k))) # Everything is in decreasing order
		return replace!([Dψ(z, x[2], x[3]) for z in Z], NaN => 0)
	end
end;
