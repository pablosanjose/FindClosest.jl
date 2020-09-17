using FindClosest
using StaticArrays, Random
using Test

@testset "FindClosest.jl" begin
    for D in 1:10, l in 2:100, i in 1:100
        Random.seed!(i)
        pts = rand(SVector{D,Float64},l)
        @test findclosest(pts) == FindClosest.bruteforce(pts)
    end
end