module scAM

    # first list all the dependencies essentially same as below
    # include("Include.jl")

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
    # include("extra_file.jl")

    # make sure to only export user functions [do not export internal functions]
    # export normalize_matrix

    # these files contain all the internal functions [no need to export these functions]
    include("Preprocessing.jl")
    include("Tools.jl")
    include("Utilities.jl")
    include("Analysis.jl")

    # These files contain all the external functions
    include("Structure.jl")
    include("Operations.jl")

    # export custom object type
    export scAMobj

    # export the external functions
    export create_scAMobj
    export percentage_set_feature
    export MetricsPlot
    export modify_scAMobj

    export merge_scAMobjs

    export Normalize
    export Scale

    export Calculate_UMAP
    export Cluster
    export plot_UMAP
    export FeaturePlot

    export FindAllDGE
    export FindDGE

end