using scAM
using Test

# using DelimitedFiles

@testset "scAM.jl" begin
    # Write your tests here.

    # println(pwd())
    data_folder = joinpath(pwd(),"data")
    # Read the files and matrices
    files = readdir(data_folder)

    # tests to make scAMobj from a single dge file
    file_path = joinpath(data_folder, files[1])
    test_obj = create_scAMobj(file_path)
    @test typeof(test_obj) == scAMobj

    # tests to make merge scAMobjs from multiple dge files
    # initialize a vector to store the scAMobjs
    scAMobjs = Vector{scAMobj}()

    # parse the directory
    for file in files

        # get the file_path
        file_path = joinpath(data_folder, file)

        # create the corresponding object
        obj = create_scAMobj(file_path)

        # add it to the vector
        push!(scAMobjs, obj)

    end

    # merge them and test validity
    scAMobj_merged = merge_scAMobjs(scAMobjs)
    @test typeof(scAMobj_merged) == scAMobj

    # Add percent mito and check if it has atleast one non zero element [init is all zeros]
    scAMobj_merged = percentage_set_feature(scAMobj_merged, "MT-")
    @test iszero(scAMobj_merged.pct_feature) == false

    # Subset the object and test validity
    scAMobj_merged = modify_scAMobj(scAMobj_merged,500, 10000, 800, 4000, 40)
    @test typeof(scAMobj_merged) == scAMobj

    # Normalize and check if it has atleast one non zero element [init is all zeros]
    scAMobj_merged = Normalize(scAMobj_merged, 10000)
    @test iszero(scAMobj_merged.normalized_counts) == false

    # Test scaling (not used though...) and check if it has atleast one non zero element [init is all zeros]
    scAMobj_merged = Scale(scAMobj_merged)
    @test iszero(scAMobj_merged.scaled_counts) == false

    # Calculate UMAP embedding and check if it has atleast one non zero element [init is all zeros]
    scAMobj_merged = Calculate_UMAP(scAMobj_merged)
    @test iszero(scAMobj_merged.umap_embedding) == false

    # Clustering and check if it has atleast one non zero element [init is all zeros]
    scAMobj_merged = Cluster(scAMobj_merged, 8)
    @test iszero(scAMobj_merged.clusters) == false

    # How to test plots??
    # Plot UMAP [if no metadata is provided, it will plot by clusters]
    p = plot_UMAP(scAMobj_merged)
    @test string(typeof(p)) == "VegaLite.VLSpec"

    p = plot_UMAP(scAMobj_merged, scAMobj_merged.sample_id)
    @test string(typeof(p)) == "VegaLite.VLSpec"

    # Test feature plot
    p = FeaturePlot(scAMobj_merged, "ALB")
    @test string(typeof(p)) == "VegaLite.VLSpec"

    # test plot_proportions
    p = plot_proportions(scAMobj_merged.clusters, scAMobj_merged.sample_id)
    @test string(typeof(p)) == "VegaLite.VLSpec"

    # test violin plots
    p = ViolinPlot(scAMobj_merged, "ALB", scAMobj_merged.clusters)
    @test string(typeof(p)) == "VegaLite.VLSpec"

    p = ViolinPlot(scAMobj_merged, "CLEC4F", scAMobj_merged.clusters)
    @test string(typeof(p)) == "VegaLite.VLSpec"

    # Run DGE for all clusters and check the results dictionary has been populated
    dge_results = FindAllDGE(scAMobj_merged)
    @test isempty(dge_results) == false

    # Run DGE for cluster 3 and check if the data frame has at least one row
    dge_res = FindDGE(scAMobj_merged, "3")
    @test isempty(dge_res) == false
end