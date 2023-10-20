
import CorticalSurfaces: vertices
export Parcel, vertices, size, length, density
export intersect, union, setdiff, getindex, setindex, setindex!
export overlap, complement

struct Parcel{T <: Hemisphere}
	surface::T
	membership::BitVector
end

"""
    Parcel(surface::Hemisphere, args...)

Make an empty `Parcel` where `surface` dictates the length of the representational space
"""
function Parcel(surface::Hemisphere)
	return Parcel{typeof(surface)}(surface, falses(size(surface)))
end

"""
	 Parcel(surface, verts)

Make a `Parcel`, given its vertex indices
"""
function Parcel(surface::Hemisphere, verts::Vector{Int})
	membership = falses(size(surface))
	membership[verts] .= true
	return Parcel{typeof(surface)}(surface, membership)
end

"""
    Parcel(surface, coords, tree)

Given a `Matrix` of arbitrary x, y, z coordinates and a `KDTree` representing the 
positions of defined cortical vertex indices, make a `Parcel` by mapping those 
coordinates to the set of defined indices via nearest neighbor search
"""
function Parcel(surface::Hemisphere, coords::AbstractMatrix, tree::KDTree)
	inds, dists = knn(tree, coords, 1)
	inds = [x[1] for x in inds] # flatten result to just a vector
	nverts = size(surface)
	all(inds .> 0) || return Parcel(surface)
	return Parcel(surface, inds)
end

"""
    vertices(p::Parcel)

Get the vertex indices belonging to a `Parcel`
"""
vertices(p::Parcel) = findall(p.membership)

"""
	 vertices(p::Parcel, Exclusive())

Get the vertex indices belonging to a `Parcel` and adjust for exclusion 
of medial wall vertices
"""
vertices(p::Parcel, ::Exclusive) = collapse(findall(p.membership), p.surface)

"""
	 vertices(p::Parcel, Bilateral(), Exclusive())

Get the vertex indices belonging to a `Parcel` and adjust for exclusion 
of medial wall vertices, and make the indexing relative to the whole brain
"""
function vertices(p::Parcel, ::Bilateral, ::Exclusive)
	temp = vertices(p.surface, Bilateral(), Exclusive())
	return temp[collapse(findall(p.membership), p.surface)]
end

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

overlap(p::Parcel, x::Union{Vector{Bool}, BitVector}) = p.membership' * x

overlap(p::Parcel, px::Parcellation) =
	return sum([overlap(p, px[k]) for k in keys(px)])

overlap(px::Parcellation, p::Parcel) =
	return overlap(p, px)

function overlap(px::Parcellation)
	counts = zeros(Int, length(px))
	for k in keys(px)
		counts .+= px[k].membership
	end
	overlap_inds = findall(counts .> 1)
	return sum(counts[overlap_inds] .- 1)
end

"""
    complement(p1::Parcel, p2::Parcel)

Compute the number of member vertices in `Parcel` `p1` not shared by those of `p2`
"""
complement(p1::Parcel, p2::Parcel) = p1.membership' * .!p2.membership

complement(p::Parcel, x::Union{Vector{Bool}, BitVector}) = p.membership' * .!x




