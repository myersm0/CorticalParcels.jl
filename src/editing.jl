
import Graphs

function Graphs.Graph(p::Parcel, A::SparseMatrixCSC)
	verts = vertices(p)
	nvertices = length(p)
	temp = spzeros(Bool, nvertices, nvertices)
	temp[verts, verts] .= true
	temp .*= A
	return Graphs.Graph(temp)
end

function cut(p::Parcel, A::SparseMatrixCSC)
	g = Graphs.Graph(p, A)
	a = Graphs.articulation(g)
	Graphs.rem_vertices!(g, a)
	clusters = filter(x -> length(x) > 1, Graphs.connected_components(g))
	new_parcels = [Parcel(c; n = length(p)) for c in clusters]
	return new_parcels
end

clear!(p::Parcel) = p.membership .*= false

function Base.merge!(p1::Parcel, p2::Parcel, A::SparseMatrixCSC)
	i = interstices(p1, p2, A)
	sum(i) > 0 || return 0
	union!(p1, i)
	union!(p1, p2)
	clear!(p2)
	return size(p1)
end

function Base.merge!(px::Parcellation{T}, k1::T, k2::T, A::SparseMatrixCSC) where T
	p1 = px[k1]
	p2 = px[k2]
	merge!(p1, p2)
	delete!(px.parcels, k2)
	return size(p1)
end


