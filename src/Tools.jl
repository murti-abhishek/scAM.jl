# function to normalize data
function normalize_data(raw_counts, scale)

	# get the total UMI count for each cell (essentially sum by column)
	total_cell_counts = sum(raw_counts, dims = 1)

	# get divide each count by the total UMI of the cell
	normalized_matrix = raw_counts ./ total_cell_counts

	# multiply by the scaling factor of your choice
	normalized_matrix = normalized_matrix .* scale

	# add 1
	normalized_matrix = normalized_matrix .+ 1

	# take natural log
	normalized_matrix = log.(normalized_matrix)
	
    # return the normalized matrix
	return(normalized_matrix)
end		


# find highly variable features
function get_gene_mean_and_var(normalized_counts)

    # get the number of genes
	(n_gene, _) = size(normalized_counts)

    # initialize the arrays for holding the mean and variance of each gene
	all_gene_var = zeros(Float64, n_gene)
	all_gene_mean = zeros(Float64,n_gene)

    # loop through the number of genes
	for row in 1:n_gene

        # get the raw counts of one gene
		gene_exp = normalized_counts[row,:]

        # calculate the mean and variance of that gene and add it to the respective array
		all_gene_mean[row] = mean(gene_exp)
		all_gene_var[row] = var(gene_exp)
	end

    # return the mean and variance vectors respectively
	return(all_gene_mean, all_gene_var)
end


# Plot variable genes
function plot_variable_genes(genes, exp_avg, exp_var)

    # create a data frame frame with genes means and variances
    df = DataFrame(Genes = genes, Avg_Exp = exp_avg, Variance = exp_var)

    # sort the data frame by decreasing order of variance
    sort!(df, [:Variance], rev = true)

    # first plot the top 10 most variable genes with names
    gr()
    @df df[1:10,:] scatter(
        :Avg_Exp,
        :Variance,
        group = :Genes,
        color = "red"
    )
    # then plot the rest
    @df df[11:end,:] scatter!(
        :Avg_Exp,
        :Variance,
        color = "black"
    )
    
end


# Scale data
function scale_data(normalized_counts)

    # believe me this works; dims = 2 is to do it by row i.e. genes
    scaled_matrix = (normalized_counts .- mean(normalized_counts,dims=2)) ./ std(normalized_counts,dims=2)

    # add functionality to clip the values to a max of 10

    # return the scaled matrix
    return(scaled_matrix)
end

## Run PCA
function run_PCA(counts, dimensions)

    # fit PCA model
    M = fit(PCA, counts; maxoutdim = dimensions)

    # return the PCA object
    return(M)
end


# function for elbow plot
function elbow_plot(PCA_obj)

    # Calculate the expalined variance of each PC
    pct_contribution = principalvars(PCA_obj) ./ tvar(PCA_obj) * 100

    # Get the cumulative variance of each PC
    pct_contribution_cumulative = cumsum(pct_contribution)

    # do the thing
    scatter(1:length(pct_contribution_cumulative), 
    pct_contribution_cumulative, 
    xlabel = "Principal Component", 
    ylab = "Cumulative Explained Variance")

end

# Function to run clustering from raw counts [make it after normalization....]
function calulate_UMAP_reduction(normalized_matrix)

    # UMAP
    Random.seed!(1)
    UMAP_reduction = umap(normalized_matrix; n_neighbors=10, min_dist=0.005, n_epochs=300)

    return(UMAP_reduction)

end

# Function to run clustering from raw counts [make it after normalization....]
function run_clustering(UMAP_reduction, n_clusters)

    # calculate the neighbours from the UMAP_reduction
    Random.seed!(1)
    clustering_result = kmeans(UMAP_reduction, n_clusters)

    return(clustering_result.assignments)

end

# Function to plot the UMAP
function make_UMAP(UMAP_reduction, clustering_result)

    # make the data frame for plotting
    df = DataFrame(UMAP_1 = UMAP_reduction[1,:], UMAP_2 = UMAP_reduction[2,:], cluster = string.(clustering_result))

    # make the plot
    p = df |> @vlplot(:point, x=:UMAP_1, y=:UMAP_2, color=:cluster)

    return(p)
end


# Function to make feature plots
function make_feature_plot(UMAP_reduction, gene_of_interest, genes, gene_exp_mtx)

    # get the index of the gene you want to plot
    index_of_interest = findall(x->x==gene_of_interest, genes)

    # make the data frame
    df = DataFrame(UMAP_1 = UMAP_reduction[1,:], UMAP_2 = UMAP_reduction[2,:], Expression = vec((gene_exp_mtx[index_of_interest,:])))

    # make the plot
    p = df |> @vlplot(:point, x=:UMAP_1, y=:UMAP_2, color=:Expression, title = gene_of_interest)

    return(p)

end

# function to plot proportions

# function to add scores

# function to create modified counts assay

# function to plot histograms