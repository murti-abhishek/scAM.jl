using scAM
using Test
#using Revise

# scAM.my_func(2,1)

@testset "scAM.jl" begin
    # Write your tests here.
    @test my_func(2,1) == 5
    @test my_func(2,2) == 6

    @test x_derivative(1,2) == 2
    @test x_derivative(2,2) == 2

end