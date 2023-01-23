include("Preprocessing.jl")

# function to make a date frame from one sample *_dge.txt file
function make_sample_df_from_path(file_path, sample_number)

    # read the matrix
    mtx = readdlm(file_path)
    #sample_name = split(split(file_path,"/")[end], "_dge.txt")[1]

    # read the counts [and convert them to Float64]
    counts = convert(Array{Float64,2}, mtx[2:end, 2:end])

    # make data frame
    df = DataFrame(counts, :auto)

    # Change colnames to cell names
    symbols = Array{String,1}(mtx[1,2:end])

    # add the sample number to the end of the column (cells)
    symbols = string.(symbols,"_",sample_number)

    # change the colnames
    rename!(df, symbols)

    # Add a column for gene name
    df[!,:Gene] = mtx[2:end,1]

    return(df)

end

# function to create a dictionary that hold the sample name as well as the corresponding data frame + a dictionary for bookkeeping (holding sample indices)
function create_sample_df_dictionary(data_directory)

    # parse the data directory and get all the files
    files = readdir(data_directory)

    # create sample dict
    sample_df_dictionary = Dict{String, DataFrame}()

    # create sample to number mapping dictionary
    sample_number_mapping_dictionary = Dict{Int64, String}()

    # initialize the sample number
    i = 1 

    # loop through the samples
    for sample in files

        # remove _dge.txt from samples (file) name
        sample_name = split(sample,"_dge.txt")[1]

        # get the sample file_path
        sample_path = joinpath(data_directory, sample)

        # create sample_df
        sample_df = make_sample_df_from_path(sample_path, i)

        # add sample to dictionary; key = sample_name; value = sample_df
        sample_df_dictionary[sample_name] = sample_df

        # add index to  mapping dictionary; key = sample_name; value = index
        sample_number_mapping_dictionary[i] = sample_name

        # +1
        i = i + 1

    end

    return(sample_df_dictionary, sample_number_mapping_dictionary)

end


# function to merge all the data frames by the Gene column
function merge_sample_df(sample_df_dictionary, sample_number_mapping_dictionary)

    # how many samples are there?
    n_samples = length(sample_number_mapping_dictionary)

    if n_samples == 1

        print("What are you trying to merge?")
        merged_df = sample_df_dictionary[sample_number_mapping_dictionary[1]]

    elseif n_samples == 2

        # merge them
        merged_df = outerjoin(sample_df_dictionary[sample_number_mapping_dictionary[1]],
        sample_df_dictionary[sample_number_mapping_dictionary[2]], 
        on = :Gene, 
        makeunique=true)
    
    else

        # if there are more than 2 samples
        # merge the first two
        merged_df = outerjoin(sample_df_dictionary[sample_number_mapping_dictionary[1]],
        sample_df_dictionary[sample_number_mapping_dictionary[2]], 
        on = :Gene, 
        makeunique=true)

        # merge the rest
        for i in 3:n_samples

            merged_df = outerjoin(merged_df, 
            sample_df_dictionary[sample_number_mapping_dictionary[i]], 
            on = :Gene, 
            makeunique=true)

        end

    end

    # replace the missing values with zeros
    for col in eachcol(merged_df)
        col = replace!(col, missing => 0.0)
    end

    return(merged_df)

end

# create a matrix from the merged data frame [in a format that Preprocessing.jl can understand]
function create_merged_mtx_from_df(merged_df)

    # get the gene column from df
    genes = merged_df.Gene

    # remove the gene name column
    only_counts = merged_df[!, Not(:Gene)]

    # get the cells
    cells = names(only_counts)

    # make sure theyre Float64
    only_counts = Matrix(only_counts)
    only_counts = convert(Array{Float64,2}, only_counts)

    # create the matrix
    mtx = create_matrix(only_counts, cells, genes)

    return(mtx)

end

# add sample ident from the
function add_sample_ident(mtx, sample_number_mapping_dictionary)

    # get the cells
    (cells, _) = get_cells_and_genes(mtx)

    # initilaize an empty vector
    sample_ids = Vector{String}()

    # loop through the cells
    for cell in cells

        # get the sample number
        sample_number = split(cell,"_")[2]

        # the sample number is a string, convert to an integer
        sample_number = parse(Int64, sample_number)

        # get the sample name associated with that sample number
        sample_id = sample_number_mapping_dictionary[sample_number]

        # append to sample_ids
        push!(sample_ids, sample_id)

    end

    return(sample_ids)

end


function make_sample_df_from_matrix(mtx, sample_number)

    # read the counts [and convert them to Float64]
    counts = convert(Array{Float64,2}, mtx[2:end, 2:end])

    # make data frame
    df = DataFrame(counts, :auto)

    # Change colnames to cell names
    symbols = Array{String,1}(mtx[1,2:end])

    # add the sample number to the end of the column (cells)
    symbols = string.(symbols,"_",sample_number)

    # change the colnames
    rename!(df, symbols)

    # Add a column for gene name
    df[!,:Gene] = mtx[2:end,1]

    return(df)

end