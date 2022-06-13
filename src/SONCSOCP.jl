module SONCSOCP

using JuMP
using Mosek
using MosekTools
using LinearAlgebra
using MultivariatePolynomials
using ECOS

export soncsocp

include("heuristic.jl")
include("SONC.jl")

end
