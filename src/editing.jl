
import Graphs
export cut, clear!, append!, merge!
 
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
a new set of `Parcel`s: one for each connected component remaining after the vertex cut
"""
function cut(p::Parcel)
	haskey(p.surface, :A) || error("Operation requires adjacency matrix `A`")
	g = Graphs.Graph(p, p.surface[:A])
	a = Graphs.articulation(g)
	Graphs.rem_vertices!(g, a)
	clusters = filter(x -> length(x) > 1, Graphs.connected_components(g))
	new_parcels = [Parcel(p.surface, c) for c in clusters]
	return new_parcels
end

"""
    clear!(p)

Zero-out all membership vertices of a `Parcel`
"""
clear!(p::Parcel) = p.membership .*= false

"""
    append!(p, v)

Add vertex `v::Int` to the `p`'s membership vector
"""
Base.append!(p::Parcel, v::Int) = p.membership[v] = true
Base.append!(p::Parcel, v::Vector{Int}) = p.membership[v] .= true

"""
    merge!(p1, p2, A)

Merge two `Parcel`s by moving the member vertices of `p2` to `p1`
"""
function Base.merge!(p1::Parcel, p2::Parcel, A::AdjacencyMatrix)
	i = interstices(p1, p2, A)
	sum(i) > 0 || return 0
	union!(p1, i)
	union!(p1, p2)
	clear!(p2)
	return size(p1)
end

"""
    merge!(px, k1, k2, A)

Given a `Parcellation{T}` and two keys of type `T`, merge the two `Parcel`s denoted
by those keys and delete the latter one from the dictionary
"""
function Base.merge!(px::Parcellation{T}, k1::T, k2::T, A::AdjacencyMatrix) where T
	p1 = px[k1]
	p2 = px[k2]
	merge!(p1, p2, A)
	delete!(px.parcels, k2)
	return size(p1)
end


