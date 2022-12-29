# create the custom object for this tool

mutable struct scAMobj

    # Think of all the information it needs to hold

    # The matrix
    matrix::Matrix

    # cells
    cells::Vector{String}

    # gene names
    genes::Vector{String}
    
    # The counts
    # raw counts
    counts::Array{Float64,2}
    # normalized counts
    normalized_counts::Array{Float64,2}
     # scaled counts
    scaled_counts::Array{Float64,2}

    # sample_id(s) [maybe be {Any}]
    sample_id::Vector{String}
    
    # UMAP embedding
    umap_embedding::Array{Float64,2}

    # Clustering assignments
    #clustering_result::KmeansResult{Matrix{Float64}, Float64, Int64}
    clusters::Vector{Int64}

    # number of features
    nFeauture::Vector{Float64}

    # number of counts
    nCount::Vector{Float64}

    # percent_feature
    pct_feature::Vector{Float64}

    # Do this as well later

    function scAMobj(matrix, cells, genes, counts, sample_id)

        # initialize the normalized and scaled counts with zeros # we'll have to check the size later
        normalized_counts = zeros(Float64, length(genes), length(cells))
        scaled_counts = zeros(Float64, length(genes), length(cells))

        # intialize the UMAP embeddings and kmeans clusters
        umap_embedding = zeros(Float64, 2, length(cells)) # umap is a 2 dimsension reduction
        clusters = zeros(Int64, length(cells))

        # intialize the nFeauture, nCount and pct_feature
        nFeauture = zeros(Float64, length(cells))
        nCount = zeros(Float64, length(cells))
        pct_feature = zeros(Float64, length(cells))

        new(matrix, cells, genes, counts, normalized_counts, scaled_counts, sample_id, umap_embedding, clusters, nFeauture, nCount, pct_feature)

    end    

end
