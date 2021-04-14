### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ c19dbeac-9ab0-11eb-1aa3-a7c673a3e0b4
begin
	using LinearAlgebra, SparseArrays, LowRankApprox, FillArrays, DataFrames, Plots
	theme(:bright)
end

# ╔═╡ ae7cded4-9ab7-11eb-135a-9d0399a2f464
atol = 1e-3; rtol = 1e-3;

# ╔═╡ 82020c5a-9ad8-11eb-188a-2f698b962e5f
function lr_compress(A::LinearOperator{T}; kest, r, atol, rtol) where T
	m, n = size(A)
	maxiter = 10
	Ω = randn(n, kest)
	S = A*Ω
	i = 0
	failed=true
	
	F = pqrfact!(S, sketch=:none)
	Q = F.Q
	
	while failed && kest ≤ min(m,n)
		i = i+1
		
		Ωtest = randn(n,r)
		Stest = A*Ωtest;
    	nrm = sqrt(1/r)*norm(Stest)
		Stest .= Stest .- Q*Q'*Stest;
    	nrm_est = sqrt(1/r)*norm(Stest)
    	failed = nrm_est > atol && nrm_est > rtol*nrm
		if failed
			F = pqrfact!(Stest, sketch=:none)
			Q = [Q F.Q]
			kest = kest+r
		end
	end
	AP = A'*Q
	return Q, AP
end

# ╔═╡ cc2c56ee-9adf-11eb-0e29-3531d7174dbd
begin
	U = randn(100,20); V = randn(100,20)
	Lop = (y, _, x) -> U*(V'*x)
	Rop = (y, _, x) -> V*(U'*x)
	Aop = LinearOperator{Float64}(100,100, Lop, Rop, nothing)
	Q, AP = lr_compress(Aop; kest = 10, r=5, atol, rtol)
	norm(Q*AP' - U*V')
end

# ╔═╡ 07450b00-9af3-11eb-1fbc-35d8e054ac6d
P = prange(Aop; sketch=:randn, atol, rtol)

# ╔═╡ d9e3251c-9ab0-11eb-2790-aba54b11d1cb
begin
	df = DataFrame()
	df[!, "n"] = [64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072]
	df[!, "average nnz per row"] = Vector{Float64}(undef, length(df[:,1]))
	df[!, "application time"] = Vector{Float64}(undef, length(df[:,1]))
	df[!, "pqr factorization time"] = Vector{Float64}(undef, length(df[:,1]))
	df[!, "prange factorization time"] = Vector{Float64}(undef, length(df[:,1]))
	df[!, "sketch factorization time"] = Vector{Float64}(undef, length(df[:,1]))
end;

# ╔═╡ 22f93552-9ab1-11eb-0b37-e52edd94cfd4
for r = 1:length(df[:,1])
	n = df.n[r]
	kest = 400
	ρ = kest/n^2
	println("Running it for n = ", n)
	#A = sprandn(n,k,ρ)
	A = sprandn(n,n,ρ)
	df[r, "average nnz per row"] = nnz(A)/n
	X = randn(n,10)
	# build operator for sampling
	Lop = (y, _, x) -> A*x
	Rop = (y, _, x) -> (x'*A)'
	Aop = LinearOperator{Float64}(size(A)..., Lop, Rop, nothing)
	_, df[r, "application time"] = @timed Aop*X
	_, df[r, "pqr factorization time"] = @timed pqrfact(Aop, sketch=:randn, atol=atol, rtol=rtol)
	_, df[r, "prange factorization time"] = @timed prange(Aop, sketch=:randn, atol=atol, rtol=rtol)
	_, df[r, "sketch factorization time"] = @timed sketchfact(Aop; atol=atol, rtol=rtol)
end

# ╔═╡ b1b5cec2-9ab1-11eb-38a1-9f048dda12cc
df

# ╔═╡ f10f0106-9ab1-11eb-0c75-dd2e54a04667
begin
	n = df[:,"n"]
	m = df[:, "average nnz per row"]
	kest = n.*m
	
	plot(xaxis=:log, yaxis=:log, legend=:topleft)
	
	plot!(df[:,"n"], df[:,"application time"], marker=true, label="application times")
	plot!(df[:,"n"], df[:,"pqr factorization time"], marker=true, label="pqr factorization times")
	plot!(df[:,"n"], df[:,"prange factorization time"], marker=true, label="prange factorization times")
	plot!(df[:,"n"], df[:,"sketch factorization time"], marker=true, label="sketch factorization time")
	
	plot!(n, 0.0000000005.*n.*kest, label=nothing)
	plot!(n, 0.000000005.*n.*(kest.^2), label=nothing)
	plot!(n, 0.000000005.*kest.*(n.^2), label=nothing)
	#plot!(df.n, 0.00000000001.*df.n.*log.(df.n).^8, label=nothing)
end

# ╔═╡ b757419c-9af4-11eb-38d6-e55c899d9668


# ╔═╡ Cell order:
# ╠═c19dbeac-9ab0-11eb-1aa3-a7c673a3e0b4
# ╠═ae7cded4-9ab7-11eb-135a-9d0399a2f464
# ╠═82020c5a-9ad8-11eb-188a-2f698b962e5f
# ╠═cc2c56ee-9adf-11eb-0e29-3531d7174dbd
# ╠═07450b00-9af3-11eb-1fbc-35d8e054ac6d
# ╠═d9e3251c-9ab0-11eb-2790-aba54b11d1cb
# ╠═22f93552-9ab1-11eb-0b37-e52edd94cfd4
# ╠═b1b5cec2-9ab1-11eb-38a1-9f048dda12cc
# ╠═f10f0106-9ab1-11eb-0c75-dd2e54a04667
# ╠═b757419c-9af4-11eb-38d6-e55c899d9668
