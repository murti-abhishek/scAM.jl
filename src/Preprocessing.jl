
# function to remove cells with zero counts
function check_dge(mtx)

    # Get the list of genes and cells
    (cells, genes) = (mtx[1,2:end], mtx[2:end,1])

    # Conver them to a string
    cells = convert(Array{String,1}, cells)
    genes = convert(Array{String,1}, genes)

    # get the gene expression matrix and convert to Float64
    gene_exp_mtx = convert(Array{Float64,2}, mtx[2:end, 2:end])

    # create a boolean vector to keep track of the good cells
    good_cells_index = Vector{Bool}()

    # get the number of cells
    (_, n_cells) = size(gene_exp_mtx)

    # loop through the cells aka columns
    for col in 1:n_cells

        # get the total UMI of the cell
        total_umi = sum(gene_exp_mtx[:,col])

        # if the total UMI > 0; good cell :)
        if total_umi != 0
        
            push!(good_cells_index, true)

        # if the total UMI = 0; bad cell :(
        else
        
            push!(good_cells_index, false)
        end

    end

    # get the fixed matrix
    fixed_gem = gene_exp_mtx[:,good_cells_index]

    # and the good cells
    cells_to_keep = cells[good_cells_index]

    # return the fixed gene expression matrix, the genes and the good cells
    return(fixed_gem, genes, cells_to_keep)
    
end

# function to set minimum cells and genes? (Maybe not, because there is subset functionality)

# function to merge counts with genes and cells and create matrix
function create_matrix(gem, cells, genes)

    # get the number of cells and genes
    (n_genes, n_cells) = size(gem)

    # Initialize a matrix of Any type with nothing
    mod_mtx = Array{Any}(nothing,n_genes+1,n_cells+1)

    # this is irrelevant but important
    mod_mtx[1,1] = "Cells_by_Genes"

    # cols as cells and rows as genes
    mod_mtx[1,2:end] = cells
    mod_mtx[2:end,1] = genes

    # the actual matrix
    mod_mtx[2:end,2:end] = gem

    # send it back
    return(mod_mtx)

end


# create new modified matrix
function modify_matrix(mtx)

    # get the fixed matrix
    (gem, genes, cells) = check_dge(mtx)

    fixed_mtx = create_matrix(gem, cells, genes)

    return(fixed_mtx)

end

# get cells and genes from matrix
function get_cells_and_genes(mtx)

    # Get the list of genes and cells
    (cells, genes) = (mtx[1,2:end], mtx[2:end,1])

    # Convert them to a string
    cells = convert(Array{String,1}, cells)
    genes = convert(Array{String,1}, genes)

    return(cells, genes)

end

function get_raw_counts(mtx)

    # get the gene expression matrix and convert to Float64
    gene_exp_mtx = convert(Array{Float64,2}, mtx[2:end, 2:end])

    return(gene_exp_mtx)

end

# funtion to display metrics
function get_info(mtx)

    (cells, genes) = get_cells_and_genes(mtx)

    # print dge info
    print(string(length(cells)), " cells")
    print(" X ")
    print(string(length(genes)), " genes")

end


# function to merge data 


# function to add feature and count metrics
function add_feature_counts_metrics(mtx)

    # get the info
    (cells, _) = get_cells_and_genes(mtx)
    gene_exp_mtx = get_raw_counts(mtx)

    # initialze an empty vector to store Counts and Features
    count_vector = Vector{Float64}()
    feature_vector = Vector{Float64}()

    # Loop through the cells
    for i in 1:length(cells)

        # get the raw counts of the cell
        cts = gene_exp_mtx[:,i]

        # get the total counts for that cell
        n_counts = sum(cts)

        # append it to the vector
        push!(count_vector, n_counts)

        # get features for that cell
        n_features = length(cts) - sum(iszero.(cts))

        # append it to the vector
        push!(feature_vector, n_features)

    end

    return(count_vector, feature_vector)

end

# function to do percentage feature set for mito genes (or any genes with the given prefix)
function add_gene_metrics(mtx, gene_prefix)

    # get the gene list
    #(gene_exp_mtx, genes, cells) = check_dge(mtx)

    # get the info
    (cells, genes) = get_cells_and_genes(mtx)
    gene_exp_mtx = get_raw_counts(mtx)

    # generate a list of all gene indices
    all_indices = collect(1:length(genes))

    # find the indices (and genes) that start with the given gene_prefix
    prefix_indices = findall(x -> startswith(x,gene_prefix), genes)
    prefix_genes = genes[prefix_indices]

    # find the other indices (and genes)
    other_genes_indices = setdiff(all_indices, prefix_indices)
    other_genes = genes[other_genes_indices]

    # initialze an empty vector to store percentage values
    prefix_pct_vector = Vector{Float64}()

    # Loop through the cells
    for i in 1:length(cells)

        # get the prefix gene counts and total counts for that cell
        prefix_counts = sum(gene_exp_mtx[prefix_indices,i])
        total_counts = sum(gene_exp_mtx[:,i])

        # get the percentage
        prefix_pct = (prefix_counts / total_counts) * 100

        # append it to the vector
        push!(prefix_pct_vector, prefix_pct)

    end

    # return the vector
    return(prefix_pct_vector)

end

# function to plot counts, features and prefix_mito
function plot_metrics(mtx)

    # get features and counts
    (n_counts, n_features) = add_feature_counts_metrics(mtx)

    # get_percent_mito
    percent_mito = add_gene_metrics(mtx, "MT-")

    # plot features and counts
    p1 = violin(n_counts, linewidth = 0, title = "Counts")
    p2 = violin(n_features, linewidth = 0, title = "Features")
    p3 = violin(percent_mito, linewidth = 0, title = "Percent Mito")
    plot(p1, p2, p3, layout = (1,3), legend = false, xaxis = nothing)

end

# function to set cutoffs and subset
function subset_dge(mtx, nCounts_low, nCounts_high, nFeatures_low, nFeatures_high, pct_mito)

    # get all info
    #(gene_exp_mtx, genes, cells) = check_dge(mtx)

    # get the info
    (cells, genes) = get_cells_and_genes(mtx)
    gene_exp_mtx = get_raw_counts(mtx)

    # generate a list of all cell indices
    all_cell_indices = collect(1:length(cells))

    # get the metrics for all cells
    (n_counts, n_features) = add_feature_counts_metrics(mtx)
    percent_mito = add_gene_metrics(mtx, "MT-")

    # collect all the bad cell indices
    low_counts_indices = findall(x -> x < nCounts_low, n_counts)
    high_counts_indices = findall(x -> x > nCounts_high, n_counts)

    low_feature_indices = findall(x -> x < nFeatures_low, n_features)
    high_feature_indices = findall(x -> x > nFeatures_high, n_features)

    low_mito_indices = findall(x -> x > pct_mito, percent_mito)

    # merge all the bad cell indices
    bad_cell_indices = union(low_counts_indices, high_counts_indices, 
    low_feature_indices, high_feature_indices, 
    low_mito_indices)

    # find the good cell indices (and genes)
    good_cell_indices = setdiff(all_cell_indices, bad_cell_indices)
    good_cells = cells[good_cell_indices]
    
    # do the subset basd on the chosen criteria
    subsetted_matrix = gene_exp_mtx[:,good_cell_indices]

    # create the new matrix
    subsetted_mtx = create_matrix(subsetted_matrix, good_cells, genes)

    # return the fixed gene expression matrix, the genes and the good cells
    # return(subsetted_matrix, genes, good_cells)
    return(subsetted_mtx, good_cell_indices)

end
