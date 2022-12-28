function my_func(x,y) 
    return(2 * x + y)
end

function x_derivative(x,y)

    derivative = ForwardDiff.derivative(x -> my_func(x,y),x)

    return(derivative)
end

function normalize_matrix(raw_counts::Array{Float64,2}, scale::Float64)

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

