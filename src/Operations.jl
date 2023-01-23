# Create scAMobj
function create_scAMobj(file_path)

    # we just need the file path
    # read the file_path and create_matrix
    mtx = readdlm(file_path)

    # do the pre filtering
    mtx = modify_matrix(mtx)

    # gather all the information needed to create the scAMobj
    (cells, genes) = get_cells_and_genes(mtx)
    counts = get_raw_counts(mtx)

    sample_file_name = splitpath(file_path)[end]
    sample_id = fill(split(sample_file_name, "_dge.txt")[1],length(cells))

    # create the object using the custructor
    obj = scAMobj(mtx, cells, genes, counts, sample_id)

    # fill in the info
    obj.matrix = mtx
    obj.cells = cells
    obj.genes = genes
    obj.counts = counts
    obj.sample_id = sample_id

    obj.nFeauture, obj.nCount = add_feature_counts_metrics(obj.matrix)

    return(obj)

end

# function for adding pct feature set [e.g. peercent mito genes]
function percentage_set_feature(obj::scAMobj, prefix)

    obj.pct_feature = add_gene_metrics(obj.matrix, prefix)

    return(obj)

end

# function to plot counts, features and prefix_mito
function MetricsPlot(obj::scAMobj)

    p1 = violin(obj.nCount, linewidth = 0, title = "Counts")
    p2 = violin(obj.nFeauture, linewidth = 0, title = "Features")
    p3 = violin(obj.pct_feature, linewidth = 0, title = "Percent Mito")

    plot(p1, p2, p3, layout = (1,3), legend = false, xaxis = nothing)

end

# Merge scAMobjects
function merge_scAMobjs(objs::Vector{scAMobj})

    # find the number of objects
    # n_objs = length(objs)

    # create df dictionary
    sample_df_dictionary = Dict{String, DataFrame}()
    sample_number_mapping_dictionary = Dict{Int64, String}() 

    # start a counter 
    i_obj = 1
    # go through it object by object
    for obj in objs

        # get the matrix from the object 
        mtx = obj.matrix

        # get the sample name [pick the first element because they're all the same]
        sample_name = obj.sample_id[1]

        # make the sample_df
        sample_df = make_sample_df_from_matrix(mtx, i_obj)

        # append data frame to sample dictionaries
        # add sample to dictionary; key = sample_name; value = sample_df
        sample_df_dictionary[sample_name] = sample_df

        # add index to  mapping dictionary; key = sample_name; value = index
        sample_number_mapping_dictionary[i_obj] = sample_name

        # +1
        i_obj = i_obj + 1

    end

    # create merged df with all samples
    merged_df = merge_sample_df(sample_df_dictionary, sample_number_mapping_dictionary)

    # create final merged matrix
    final_mtx = create_merged_mtx_from_df(merged_df)

    # fixed_matrix
    fixed_mtx = modify_matrix(final_mtx)

    # get the cells, genes and counts
    (cells, genes) = get_cells_and_genes(final_mtx)
    counts = get_raw_counts(final_mtx)

    # get the sample ids
    sample_ids = add_sample_ident(fixed_mtx, sample_number_mapping_dictionary)

    # create the final object
    final_obj = scAMobj(fixed_mtx, cells, genes, counts, sample_ids)

    final_obj.nFeauture, final_obj.nCount = add_feature_counts_metrics(final_obj.matrix)

    return(final_obj)

end



# Modify after QC
function modify_scAMobj(obj::scAMobj, nCounts_low, nCounts_high, nFeatures_low, nFeatures_high, pct_mito)

    # do the subsetting
    (obj.matrix, good_cell_indices) = subset_dge(obj.matrix, nCounts_low, nCounts_high, nFeatures_low, nFeatures_high, pct_mito)

    # recreate the object (keep the indices of the good cells)
    (obj.cells, obj.genes) = get_cells_and_genes(obj.matrix)
    obj.counts = get_raw_counts(obj.matrix)
    obj.sample_id = obj.sample_id[good_cell_indices]

    obj.nFeauture, obj.nCount = add_feature_counts_metrics(obj.matrix)
    
    # This will work too
    # obj.nFeauture = obj.nFeauture[good_cell_indices]
    # obj.nCount = obj.nCount[good_cell_indices]

    obj.pct_feature = obj.pct_feature[good_cell_indices]

    return(obj)

end


# Normalize
function Normalize(obj::scAMobj, scale)

    obj.normalized_counts = normalize_data(obj.counts, scale)

    return(obj)

end

# Scale
function Scale(obj::scAMobj)

    obj.scaled_counts = scale_data(obj.normalized_counts)

    return(obj)

end

# Calulate_UMAP_reduction
function Calculate_UMAP(obj::scAMobj)

    obj.umap_embedding = calulate_UMAP_reduction(obj.normalized_counts)

    return(obj)

end

# Cluster
function Cluster(obj::scAMobj, n_clusters)

    obj.clusters = run_clustering(obj.umap_embedding, n_clusters)

    return(obj)

end

# Plot UMAP
function plot_UMAP(obj::scAMobj)

    make_UMAP(obj.umap_embedding, obj.clusters)

end

function plot_UMAP(obj::scAMobj, metadata)

    make_UMAP(obj.umap_embedding, metadata)

end

# feature plot
function FeaturePlot(obj::scAMobj, gene_of_interest)

    make_feature_plot(obj.umap_embedding, gene_of_interest, obj.genes, obj.normalized_counts)

end

# Run DGE
function FindAllDGE(obj::scAMobj)

    # get all the info
    cells = obj.cells
    genes = obj.genes
    counts = obj.normalized_counts
    clusters = string.(obj.clusters)

    # make the cluster df
    df = make_cluster_df(cells, clusters)

    # create empty dict for storing dge results
    dge_results_dictionary = Dict{String, DataFrame}()

    # loop through the clusters
    for cluster in unique(clusters)

        # do dge for that cluster
        test_df = perform_dge_by_cluster(df, genes, counts, cluster)

        # save the dge results; key = cluster, value = dge data frame
        dge_results_dictionary[cluster] = test_df

    end
    
    return(dge_results_dictionary)

end

function FindDGE(obj::scAMobj, cluster_of_interest::String)

    # test for all clusters
    cells = obj.cells
    genes = obj.genes
    counts = obj.normalized_counts
    clusters = string.(obj.clusters)

    # make the data frame
    df = make_cluster_df(cells, clusters)

    # do dge
    test_df = perform_dge_by_cluster(df, genes, counts, cluster_of_interest)
    
    # return as df
    return(test_df)
end