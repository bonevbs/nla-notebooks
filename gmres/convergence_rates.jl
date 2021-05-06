### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ a1d3cf28-ad70-11eb-1515-252b85369c50
begin
	using LinearAlgebra, SparseArrays, IterativeSolvers, Plots, Random, Distributions
	theme(:bright)
	Random.seed!(123)
	pgfplotsx()
end;

# ╔═╡ fb965b16-ad70-11eb-2c9a-5bf8d0002e5a
md"
## Setup
"

# ╔═╡ 75cf84e8-ad71-11eb-297b-0351d933895b
n=200;

# ╔═╡ 7151ded8-ad86-11eb-0fde-a781cb6ae0cf
md"Generate some random spectrum"

# ╔═╡ 02e67324-ad71-11eb-198e-116a05e7daf4
Λ = randn(ComplexF64, n) .+ 3.0

# ╔═╡ 540f84b0-ad86-11eb-14cd-2187093cf8d7
begin
	fig_spectrum = scatter(Λ, xlim=(0,5), ylim=(-2.5,2.5))
	savefig(fig_spectrum, "spectrum.tex")
	fig_spectrum
end

# ╔═╡ 6ca595a6-ad85-11eb-1897-8d33b04bf1b0
begin
	V = randn(n,n)
	V, R = qr(V);
	A = V*diagm(Λ)*V'
	b = randn(n,1)
end;

# ╔═╡ 8b759070-ad86-11eb-10af-d3cf79c2a590
md"
## run GMRES
"

# ╔═╡ 49b677bc-ad86-11eb-3f26-478c2d8b816b
D = eigvals(A)

# ╔═╡ 86f51672-ad86-11eb-2514-23e5123d7947
_, hist = gmres(A, b; reltol=1e-9, restart=100, log=true, maxiter=100)

# ╔═╡ 2ccc46be-ad87-11eb-2fc5-2f75a7558d2b
r = hist[:resnorm]./hist[:resnorm][1]

# ╔═╡ d169e590-ad86-11eb-2bf0-0be4e2189136
begin
	fig_gmres = plot(r, yaxis=:log, marker=true, ylim=(1E-10, 1))
	savefig(fig_gmres, "gmres_iters.tex")
	fig_gmres
end

# ╔═╡ Cell order:
# ╟─fb965b16-ad70-11eb-2c9a-5bf8d0002e5a
# ╠═a1d3cf28-ad70-11eb-1515-252b85369c50
# ╠═75cf84e8-ad71-11eb-297b-0351d933895b
# ╟─7151ded8-ad86-11eb-0fde-a781cb6ae0cf
# ╠═02e67324-ad71-11eb-198e-116a05e7daf4
# ╠═540f84b0-ad86-11eb-14cd-2187093cf8d7
# ╠═6ca595a6-ad85-11eb-1897-8d33b04bf1b0
# ╟─8b759070-ad86-11eb-10af-d3cf79c2a590
# ╠═49b677bc-ad86-11eb-3f26-478c2d8b816b
# ╠═86f51672-ad86-11eb-2514-23e5123d7947
# ╠═2ccc46be-ad87-11eb-2fc5-2f75a7558d2b
# ╠═d169e590-ad86-11eb-2bf0-0be4e2189136
