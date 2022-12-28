using scAM
using Test
#using Revise

# scAM.my_func(2,1)

@testset "scAM.jl" begin
    # Write your tests here.
    @test my_func(2,1) == 5
    @test my_func(2,2) == 6
end