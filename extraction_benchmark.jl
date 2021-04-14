### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 9976914c-9bc5-11eb-3b13-abe73dd8f7c4
using SparseArrays, BenchmarkTools, DataFrames, Random, Plots

# ╔═╡ a8f932f0-9c6a-11eb-3d59-fbaba6e0f95f
md"
# Extraction of submatrices from sparse matrices in CSC format
"

# ╔═╡ 9f8aec8c-9bd8-11eb-26ec-197ad80da9f2
theme(:bright)

# ╔═╡ a69cb3b0-9bc5-11eb-28c8-eb8af4cc451c
begin
	df = DataFrame()
	df[!, "n"] = @. 2^(8:25)
	df[!, "average nnz per row"] = Vector{Float64}(undef, length(df[:,1]))
	df[!, "extraction time"] = Vector{Float64}(undef, length(df[:,1]))
	df[!, "memory used"] = Vector{Float64}(undef, length(df[:,1]))
	df[!, "extraction time (sorted)"] = Vector{Float64}(undef, length(df[:,1]))
	df[!, "memory used (sorted)"] = Vector{Float64}(undef, length(df[:,1]))
end;

# ╔═╡ 233cec96-9bd0-11eb-0945-ab11d417df3a
begin
	nrow = 100
	k = 100
end;

# ╔═╡ 05952304-9bd0-11eb-2680-4be8ef22e8fc
for r = 1:length(df[:,1])
	n = df.n[r]
	ρ = nrow/n
	println("Running it for n = ", n)
	A = sprandn(n,n,ρ)
	df[r, "average nnz per row"] = nnz(A)/n
	i = randperm(n)[1:k]; j = randperm(n)[1:k]
	b = @benchmark $A[$i,$j]
	println(b)
	df[r, "extraction time"] = mean(b.times)
	df[r, "memory used"] = b.memory
	i = collect(100:200)
	b = @benchmark  SparseArrays.getindex_I_sorted_bsearch_I($A, $i, $j)
	df[r, "extraction time (sorted)"] = mean(b.times)
	df[r, "memory used (sorted)"] = b.memory
end

# ╔═╡ 3f715068-9bd1-11eb-169b-cbad9b57fe1f
df

# ╔═╡ 7df7f0ae-9bd8-11eb-2dde-effd6140752d
begin
	n = df[:,"n"]
	m = df[:, "average nnz per row"]
	kest = n.*m
	
	plot(xaxis=:log, yaxis=:log, legend=:topleft)
	
	plot!(n, df[:,"extraction time"], marker=true, label="extraction time")
	plot!(n, df[:,"memory used"], marker=true, label="memory usage")
	plot!(n, df[:,"extraction time (sorted)"], marker=true, label="extraction time (sorted)")
	plot!(n, df[:,"memory used (sorted)"], marker=true, label="memory usage (sorted)")
	
	plot!(n,10*n)
end

# ╔═╡ 5acbe6ee-9c89-11eb-184f-03356faa27ba
begin
	nn = 10000
	ii = randperm(nn)[1:k]
	jj = randperm(nn)[1:k]
	X = sprandn(nn, nn, nrow/nn)
	b = @benchmark SparseArrays.getindex_I_sorted_bsearch_I(X, ii, jj)
end

# ╔═╡ e1b01fec-9caa-11eb-353c-e52aeba86885
@. 2^(8:25)

# ╔═╡ bcafbc24-9d0a-11eb-0ad9-55f67d4d508f


# ╔═╡ Cell order:
# ╟─a8f932f0-9c6a-11eb-3d59-fbaba6e0f95f
# ╠═9976914c-9bc5-11eb-3b13-abe73dd8f7c4
# ╠═9f8aec8c-9bd8-11eb-26ec-197ad80da9f2
# ╠═a69cb3b0-9bc5-11eb-28c8-eb8af4cc451c
# ╠═233cec96-9bd0-11eb-0945-ab11d417df3a
# ╠═05952304-9bd0-11eb-2680-4be8ef22e8fc
# ╠═3f715068-9bd1-11eb-169b-cbad9b57fe1f
# ╠═7df7f0ae-9bd8-11eb-2dde-effd6140752d
# ╠═5acbe6ee-9c89-11eb-184f-03356faa27ba
# ╠═e1b01fec-9caa-11eb-353c-e52aeba86885
# ╠═bcafbc24-9d0a-11eb-0ad9-55f67d4d508f
