
import Graphs

Base.setindex!(p::Parcel, args...) = setindex!(p.membership, args...)

Base.dotview(p::Parcel, args...) = view(p.membership, args...)

function Base.setindex!(px::HemisphericParcellation{T}, p::Parcel, k::T) where T
	p.surface == px.surface || error("Surface incompatibility")
	px.parcels[k] = p
end
 
function Graphs.Graph(p::Parcel, A::AdjacencyMatrix)
	verts = vertices(p)
	nvertices = length(p)
	temp = spzeros(Bool, nvertices, nvertices)
	temp[verts, verts] .= true
	temp .*= A
	return Graphs.Graph(temp)
end

"""
    cut(p, A)

Cut articulation point(s), if any, from a graph representation of a `Parcel`, and return
a new set of `Parcel`s: one for each connected component remaining after the vertex cut.
"""
function cut(p::Parcel, A::AdjacencyMatrix)
	g = Graphs.Graph(p, A)
	a = Graphs.articulation(g)
	Graphs.rem_vertices!(g, a)
	clusters = filter(x -> length(x) > 1, Graphs.connected_components(g))
	new_parcels = [Parcel(p.surface, c) for c in clusters]
	return new_parcels
end

function cut(p::Parcel)
	haskey(p.surface, :A) || error("Operation requires adjacency matrix `:A`")
	return cut(p, p.surface[:A])
end

# TODO: this function is almost exactly the same as cut(p); not sure how to
# refactor though without type piracy
"""
	 split(p, v)

Remove vertices `v` from a graph representation of a `Parcel`, and return
a new set of `Parcel`s: one for each connected component remaining.
"""
function Base.split(p::Parcel, v::AbstractVector{Integer})
	haskey(p.surface, :A) || error("Operation requires adjacency matrix `:A`")
	g = Graphs.Graph(p, p.surface[:A])
	Graphs.rem_vertices!(g, v)
	clusters = filter(x -> length(x) > 1, Graphs.connected_components(g))
	new_parcels = [Parcel(p.surface, c) for c in clusters]
	return new_parcels
end

"""
    clear!(p)

Zero-out all membership vertices of a `Parcel`.
"""
clear!(p::Parcel) = p.membership .*= false

"""
    delete!(px, k)

Delete `Parcel` with ID `k` from a `Parcellation`.
"""
Base.delete!(px::HemisphericParcellation{T}, k::T) where T = delete!(px.parcels, k)

"""
    append!(p, v)

Add vertex `v::Int` to the `p`'s membership vector.
"""
Base.append!(p::Parcel, v::Integer) = p.membership[v] = true
Base.append!(p::Parcel, v::AbstractVector{Integer}) = p.membership[v] .= true

"""
    merge!(p1, p2, A)

Merge two `Parcel`s by moving the member vertices of `p2` to `p1`.
"""
function Base.merge!(p1::Parcel, p2::Parcel, A::AdjacencyMatrix)
	i = interstices(p1, p2, A)
	sum(i) > 0 || return 0
	union!(p1, i)
	union!(p1, p2)
	clear!(p2)
	return size(p1)
end

function Base.merge!(p1::Parcel, p2::Parcel)
	haskey(p1.surface, :A) || error("Operation requires adjacency matrix `A`")
	p1.surface == p2.surface || error("Surfaces must be the same for both parcels")
	return merge!(p1, p2, p1.surface[:A])
end

"""
    merge!(px, k1, k2, A)

Given a `Parcellation{T}` and two keys of type `T`, merge the two `Parcel`s denoted
by those keys and delete the latter one from the dictionary.
"""
function Base.merge!(px::HemisphericParcellation{T}, k1::T, k2::T, A::AdjacencyMatrix) where T
	p1 = px[k1]
	p2 = px[k2]
	merge!(p1, p2, A)
	delete!(px, k2)
	return size(p1)
end

function Base.merge!(px::HemisphericParcellation{T}, k1::T, k2::T) where T
	haskey(px.surface, :A) || error("Operation requires adjacency matrix `A`")
	return merge!(px, k1, k2, px.surface[:A])
end

"""
    deepcopy(p)

Make a new `Parcel` containing a `deepcopy` of original parcel `p`'s `membership` 
vector. Note however that the surface remains just a reference and is not itself 
copied, since it may be a large object.
"""
function Base.deepcopy(p::Parcel)
	return Parcel(p)
end

"""
    deepcopy(px)

Make a new parcellation containing a `deepcopy` of all parcels from `px`. Note however
that, as with `deepcopy(p::Parcel)`, the surface remains just a reference and is not 
itself copied, since it may be a large object.
"""
function Base.deepcopy(px::HemisphericParcellation{T}) where T
	px′ = HemisphericParcellation{T}(px.surface)
	for k in keys(px)
		px′[k] = deepcopy(px[k])
	end
	return px′
end

