function make_cluster_df(cells, clusters)

    # create a df with the cells and clusters
    df = DataFrame(cells = cells, cluster = clusters, index = collect(1:length(cells)))

    return(df)

end

function perform_dge_by_cluster(df, genes, counts, cluster_of_interest)

    # subset on the cluster of interest
    df_of_interest = df[df.cluster .== cluster_of_interest, :]
    # get the cells
    # cells_of_interest = df_of_interest.cells
    # get the indices of interest
    indices_of_interest = df_of_interest.index

    # subset on the other clusters
    # df_other = subset(df, :cluster => x -> x .!= cluster_of_interest)
    df_other = df[df.cluster .!= cluster_of_interest, :]
    # get the other cells
    # scells_other = df_other.cells
    # get the other indices
    indices_other = df_other.index

    # get the corresponding matrices
    matrix_of_interest = counts[:, indices_of_interest]
    matrix_other = counts[:, indices_other]

    # initialze an empty vector to store log2FC and p_val
    log2FC_vector = Vector{Float64}()
    pval_vector = Vector{Float64}()

    # initialze an empty vector to pct expression
    pct_of_interest_vector = Vector{Float64}()
    pct_other_vector = Vector{Float64}()

    # loop through the Genes
    for i in 1:length(genes)

        # get the gene expression of interest and other
        gene_exp_of_interest = matrix_of_interest[i,:]
        gene_exp_other = matrix_other[i,:]

        # calculate the coressponding average gene expression 
        avg_gene_exp_of_interest = mean(gene_exp_of_interest)
        avg_gene_exp_other = mean(gene_exp_other)

        # calculate the log2FC
        log2FC = log2(avg_gene_exp_of_interest / avg_gene_exp_other)

        # append it to the log2FC_vector
        push!(log2FC_vector, log2FC)

        # calculate the p_val
        pval = pvalue(MannWhitneyUTest(gene_exp_of_interest, gene_exp_other))

        # append it to the pval_vector
        push!(pval_vector, pval)

        # calculate pct expressions
        pct_of_interest = (length(gene_exp_of_interest) - sum(iszero.(gene_exp_of_interest))) / length(gene_exp_of_interest)
        pct_other = (length(gene_exp_other) - sum(iszero.(gene_exp_other))) / length(gene_exp_other)

        # append them
        push!(pct_of_interest_vector, pct_of_interest)
        push!(pct_other_vector, pct_other)

    end

    # create a data frame of the results
    dge_df = DataFrame(Gene = genes, Log2FC = log2FC_vector, p_value = pval_vector, PctExp1 = pct_of_interest_vector, PctExp2 = pct_other_vector)

    # Filter and remove all the useless values
    dge_df = filter(:Log2FC => x -> !any(f -> f(x), (ismissing, isnothing, isnan, isinf)), dge_df)

    # sort by log2FC by default
    sort!(dge_df, order(:Log2FC, rev = true))

    return(dge_df)

end
