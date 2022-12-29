using scAM
using Test

@testset "scAM.jl" begin
    # Write your tests here.

    # Read the matrix
    test_matrix = [[1,2,3] [4,5,6] [7,8,9]]
    result = normalize_matrix(test_matrix, 10)
    @test typeof(result) == Matrix{Float64}

end