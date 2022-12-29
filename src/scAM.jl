module scAM

# first list all the dependencies

using ForwardDiff # remove this later
using MultivariateStats
using Plots
using DelimitedFiles
using DataFrames
using LinearAlgebra
using Statistics
using StatsPlots
using UMAP
using Clustering
using NearestNeighbors
using VegaLite
using Random
using XLSX
using HypothesisTests

# Write your package code here
# then add all the files (include all the files that include internal + external functions)
include("extra_file.jl")

# amke sure to only export user functions [do not export iternal functions]
export normalize_matrix

end
