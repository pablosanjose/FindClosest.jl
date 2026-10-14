using FindClosest
using StaticArrays, Random
using Test

@testset "FindClosest.jl" begin
    for D in 1:10, l in 2:100, i in 1:100
        Random.seed!(i)
        pts = rand(SVector{D,Float64},l)
        # The distances come from separately compiled `@simd` reductions, which may
        # round differently, so only the indices are compared exactly.
        d, ij = findclosest(pts)
        d_ref, ij_ref = FindClosest.bruteforce(pts)
        @test ij == ij_ref
        @test d ≈ d_ref
    end
end