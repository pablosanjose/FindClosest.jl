module FindClosest

using Distances

export findclosest

const BRUTEFORCE_THRESHOLD = 33

"""
    findclosest(pts[, dist = Euclidean()]) -> (d, (i, j))

Find the distance `d` and indices `(i, j)` of the two closest points `pᵢ` and `pⱼ` in `pts`,
according to distance `dist(pᵢ, pⱼ)` (see the Distance.jl package). If `length(pts) < 2`,
return `nothing`.

The search is performed using a single-threaded, multidimensional, divide-and-conquer
algorithm, with runtime complexity `O(n log(n))`.

# Examples
```
julia> using StaticArrays; pts = rand(SVector{3,Float64}, 10^3);

julia> findclosest(pts)
(0.007292429825826281, (75, 535))
```
"""
function findclosest(pts, dist = Euclidean())
    scratch = scratch_arrays(pts)
    pts´ = last(scratch)
    copy_itr!(pts´, enumerate(pts))

    len = length(pts´)
    len < 2 && return nothing
    len <= BRUTEFORCE_THRESHOLD && return bruteforce(pts´, dist)

    dim = length(scratch)
    sort_dim!(pts´, dim)
    T = typeof(dist_inds(pts´, (firstindex(pts´), firstindex(pts´)), dist))
    return findclosest_sorted(pts´, dist, T, scratch)
end

function findclosest_sorted(pts, dist, ::Type{T}, scratch::NTuple{D})::T where {D,T}
    len = length(pts)
    if len <= BRUTEFORCE_THRESHOLD
        d = bruteforce(pts, dist)
    elseif D == 1
        d0 = dist_inds(pts, (firstindex(pts), lastindex(pts)), dist)
        d = findclosest_sorted_1D(pts, dist, d0)
    else
        half = len ÷ 2
        ptsL = view(pts, 1:half)
        ptsR = view(pts, half + 1:len)
        dL = findclosest_sorted(ptsL, dist, T, scratch)
        dR = findclosest_sorted(ptsR, dist, T, scratch)
        dLR = min_dinds(dL, dR)
        δ = first(dLR)
        inter_range = inter_range_dim(ptsL, ptsR, δ, D, half, len)
        if length(inter_range) < 2
            d = dLR
        else
            scratch´ = Base.front(scratch)
            ptsLR = copy_view!(last(scratch´), pts, inter_range)
            sort_dim!(ptsLR, D - 1)
            if D - 1 == 1
                dinter = findclosest_sorted_1D(ptsLR, dist, dLR)
            else
                dinter = findclosest_sorted(ptsLR, dist, T, scratch´)
            end
            d = min_dinds(dLR, dinter)
        end
    end
    return d
end

function findclosest_sorted_1D(pts, dist, dLR)
    d = dLR
    ifirst, ilast = firstindex(pts), lastindex(pts)
    for i in ifirst:ilast-1
        j = i + 1
        while j <= ilast && first(last(pts[j])) - first(last(pts[i])) < first(d)
            d = min_dinds(d, dist_inds(pts, (i, j), dist))
            j += 1
        end
    end
    return d
end

function bruteforce(ps, dist = Euclidean())
    ifirst, ilast = firstindex(ps), lastindex(ps)
    d = dist_inds(ps, (ifirst, ilast), dist)
    for i in ifirst:ilast-1, j in i+1:ilast
        d = min_dinds(d, dist_inds(ps, (i, j), dist))
    end
    return d
end

## Auxiliary functions

function scratch_arrays(pts)
    p = first(pts)
    T = typeof(p)
    D = length(p)
    scratch = ntuple(i -> Vector{Tuple{Int,T}}(undef, 0), D)
    return scratch
end

function dist_inds(ps::AbstractVector{<:Tuple}, (i, j), dist)
    (i1, p1), (i2, p2) = ps[i], ps[j]
    return dist(p1, p2), sortedtuple(i1, i2)
end

dist_inds(ps, (i, j), dist) = dist(ps[i], ps[j]), sortedtuple(i, j)

min_dinds((d, is), (d´, is´)) = ifelse(d < d´, (d, is), (d´, is´))

sortedtuple(i, j) = ifelse(i < j, (i, j), (j, i))

sort_dim!(pts, dim) = sort!(pts; alg = QuickSort, by = p -> last(p)[dim])

center_dim((i, p), (i´, p´), dim) = (p[dim] + p´[dim])/2

function inter_range_dim(ptsL, ptsR, δ, D, half, len)
    center = center_dim(last(ptsL), first(ptsR), D)
    rmax = findfirst(p -> last(p)[D] > center + δ, ptsR)
    rmin = findlast( p -> last(p)[D] < center - δ, ptsL)
    inter_range =
        (rmin === nothing ? 1 : rmin + 1):(rmax === nothing ? len : half + rmax - 1)
    return inter_range
end

function copy_itr!(dst, src)
    resize!(dst, 0)
    @inbounds for s in src
        push!(dst, s)
    end
    return dst
end

function copy_view!(dst, src, inter_range)
    len = length(inter_range)
    resize!(dst, max(length(dst), len))
    dstview = view(dst, 1:len)
    copy!(dstview, view(src, inter_range))
    return dstview
end

end # module