# FindClosest.jl

[![Build Status](https://github.com/pablosanjose/FindClosest.jl/workflows/CI/badge.svg)](https://github.com/pablosanjose/FindClosest.jl/actions)

The Findclosest.jl package exports a single function `findclosest(pts[, dist])` that finds the closest pair of points in `pts`. It employs a standard multi-dimensional divide-and-conquer recursive algorithm with optimal scaling `O(n log(n))` in the number of points `n`.

## Example
```julia
julia> using FindClosest, StaticArrays

julia> pts = rand(SVector{3,Float64}, 10^5);

julia> @time findclosest(pts)
  0.083620 seconds (29 allocations: 5.049 MiB)
(0.00018927587118480567, (4517, 34872))
```
The returned tuple has the form `(d, (i, j))`, where `(i, j)` are the indices of the two closest points in `pts` and `d` is their distance `dist(pᵢ, pⱼ)`.

Possible future developments include multithreading, non-recursive algorithms and `O(1)` memory usage