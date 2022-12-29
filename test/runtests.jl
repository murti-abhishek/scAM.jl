using scAM
using Test

using DelimitedFiles

@testset "scAM.jl" begin
    # Write your tests here.
    data_folder = string(pwd(),"/data/")
    # Read the files and matrices
    files = readdir(data_folder)

    file_path = string(data_folder, files[1])
    test_matrix = readdlm(file_path)
    test_matrix = test_matrix[2:end, 2:end]
    result = normalize_matrix(test_matrix, 10)
    @test typeof(result) == Matrix{Float64}

    file_path = string(data_folder, files[2])
    test_matrix = readdlm(file_path)
    test_matrix = test_matrix[2:end, 2:end]
    result = normalize_matrix(test_matrix, 10)
    @test typeof(result) == Matrix{Float64}

    file_path = string(data_folder, files[3])
    test_matrix = readdlm(file_path)
    test_matrix = test_matrix[2:end, 2:end]
    result = normalize_matrix(test_matrix, 10)
    @test typeof(result) == Matrix{Float64}

end