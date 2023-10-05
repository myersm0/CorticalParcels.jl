
import CorticalSurfaces: vertices
export Parcel, vertices, size, length, density
export intersect, union, setdiff, getindex, setindex, setindex!
export overlap, complement

struct Parcel
	membership::BitVector
end

"""
    Parcel(n::Int)

Make an empty `Parcel` within a representational space of `n` vertices
"""
function Parcel(n::Int)
	return Parcel(falses(n))
end

"""
    Parcel(surface::SurfaceSpace, args...)

Make an empty `Parcel` where `surface` dictates the length of the representational space
"""
function Parcel(surface::SurfaceSpace, args...)
	return Parcel(size(surface, args...))
end

"""
    Parcel(verts::Vector{Int}; n::Int)

Make a `Parcel`, given its vertex indices within a representational space of length `n`
"""
function Parcel(verts::Vector{Int}; n::Int)
	temp = falses(n)
	@inbounds temp[verts] .= true
	return Parcel(temp)
end

"""
    Parcel(coords::Matrix, tree::KDTree)

Given a `Matrix` of arbitrary x, y, z coordinates and a `KDTree` representing the 
positions of defined cortical vertex indices, make a `Parcel` by mapping those 
coordinates to the set of defined indices via nearest neighbor search
"""
function Parcel(coords::AbstractMatrix, tree::KDTree)
	inds, dists = knn(tree, coords, 1)
	inds = [x[1] for x in inds] # flatten result to just a vector
	nverts = size(tree.data, 1)
	all(inds .> 0) || return Parcel(nverts)
	return Parcel(inds; n = nverts)
end

"""
    vertices(p::Parcel)

Get the vertex indices belonging to a `Parcel`
"""
vertices(p::Parcel) = findall(p.membership)

"""
    size(p::Parcel)

Get the size (number of non-zero vertices) of a `Parcel`"
"""
Base.size(p::Parcel) = sum(p.membership)

"""
    length(p::Parcel)

Get the length of the representational space in which a `Parcel` is located
"""
Base.length(p::Parcel) = length(p.membership)

"""
    density(p::Parcel)

Get the proportion of member vertices in a `Parcel` relative to the length of its space
"""
density(p::Parcel) = size(p) / length(p)

Base.intersect(p1::Parcel, p2::Parcel) = p1.membership .& p2.membership
Base.union(p1::Parcel, p2::Parcel) = p1.membership .| p2.membership
Base.setdiff(p1::Parcel, p2::Parcel) = p1.membership .& .!p2.membership

Base.intersect!(p1::Parcel, p2::Parcel) = p1.membership .&= p2.membership
Base.union!(p1::Parcel, p2::Parcel) = p1.membership .|= p2.membership
Base.setdiff!(p1::Parcel, p2::Parcel) = p1.membership .&= .!p2.membership

Base.intersect(p::Parcel, x::BitVector) = p.membership &= x
Base.union(p::Parcel, x::BitVector) = p.membership |= x
Base.setdiff(p::Parcel, x::BitVector) = p.membership &= .!x

Base.intersect!(p::Parcel, x::BitVector) = p.membership .&= x
Base.union!(p::Parcel, x::BitVector) = p.membership .|= x
Base.setdiff!(p::Parcel, x::BitVector) = p.membership .&= .!x

Base.intersect!(p::Parcel, x::Vector{T}) where T <: Integer = p.membership .*= x
Base.union!(p::Parcel, x::Vector{T}) where T <: Integer = p.membership[x] .= true
Base.setdiff!(p::Parcel, x::Vector{T}) where T <: Integer = p.membership[x] .= false

Base.getindex(p::Parcel, args...) = getindex(p.membership, args...)
Base.setindex!(p::Parcel, args...) = setindex!(p.membership, args...)
Base.dotview(p::Parcel, args...) = view(p.membership, args...)

"""
    overlap(p1::Parcel, p2::Parcel)

Compute the number of member vertices shared between two `Parcel`s `p1`, `p2`
"""
overlap(p1::Parcel, p2::Parcel) = p1.membership' * p2.membership

"""
    complement(p1::Parcel, p2::Parcel)

Compute the number of member vertices in `Parcel` `p1` not shared by those of `p2`
"""
complement(p1::Parcel, p2::Parcel) = p1.membership' * .!p2.membership

"""
    centroid(p, dmat)

Find the centroid of a parcel (the vertex that has the least summed distance
to all other vertices in the parcel). `dmat` is expected to be
a square distance matrix of dimensions (length(p), length(p)).
"""
function centroid(p::Parcel, dmat::AbstractMatrix)
	all(size(dmat) .== length(p)) || error(DimensionMismatch)
	verts = vertices(p)
	dists = @inbounds sum(dmat[verts, verts]; dims = 1)
	return verts[argmin(dists)]
end


