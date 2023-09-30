
import CorticalSurfaces: vertices
export Parcel, vertices, size, length, intersect, union, setdiff, getindex, setindex

struct Parcel
	membership::SparseVector{Bool, Int}
end

"""
    Parcel(n::Int)

Make an empty `Parcel` within a representational space of `n` vertices
"""
function Parcel(n::Int)
	return Parcel(spzeros(Bool, n))
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
	temp = spzeros(Bool, n)
	@inbounds temp[verts] .= true
	return Parcel(temp)
end

"""
    Parcel(coords::Matrix, tree::KDTree)

Given a `Matrix` of arbitrary x, y, z coordinates and a `KDTree` representing the 
positions of defined cortical vertex indices, make a `Parcel` by mapping those 
coordinates to the set of defined indices via nearest neighbor search
"""
function Parcel(coords::Matrix, tree::KDTree)
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
Base.size(p::Parcel) = sum(p.membership.nzval)

"""
    nnz(p::Parcel)

Get the size (number of non-zero vertices) of a `Parcel`"
"""
SparseArrays.nnz(p::Parcel) = size(p)


"""
    length(p::Parcel)

Get the length of the representational space in which a `Parcel` is located
"""
Base.length(p::Parcel) = p.membership.n

Base.intersect(p1::Parcel, p2::Parcel) = (p1.membership .& p2.membership).nzind
Base.union(p1::Parcel, p2::Parcel) = (p1.membership .| p2.membership).nzind
Base.setdiff(p1::Parcel, p2::Parcel) = (p1.membership .& .!p2.membership).nzind

Base.getindex(p::Parcel, args...) = getindex(p.membership, args...)
Base.setindex!(p::Parcel, args...) = setindex!(p.membership, args...)
Base.dotview(p::Parcel, args...) = view(p.membership, args...)

function Base.show(io::IO, ::MIME"text/plain", p::Parcel)
	print(io, "Parcel with $(size(p)) vertices")
end

