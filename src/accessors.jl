
import CorticalSurfaces: vertices
export vertices, size, length, keys, haskey, values, getindex
export vec, union, unassigned, nnz, density

# ===== Parcel functions =====

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

Base.getindex(p::Parcel, args...) = getindex(p.membership, args...)

Base.isequal(p1::Parcel, p2::Parcel) = 
	p1.surface == p2.surface && p1.membership == p2.membership

Base.isequal(p1::Parcel, x::BitVector) = 
	p1.membership == p2.membership

membership(p::Parcel) = p.membership

surface(p::Parcel) = p.surface

# ===== Parcellation functions =====

"""
    size(px::Parcellation)

Get the number of Parcels comprising a Parcellation
"""
Base.size(px::Parcellation) = length(px.parcels)

"""
    length(px::Parcellation)

Get the number of vertices comprising the representation space of a `Parcellation`
"""
Base.length(px::Parcellation) = size(px.surface)

"""
    keys(px::Parcellation)

Get the IDs of all `Parcel`s within a `Parcellation`
"""
Base.keys(px::Parcellation) = keys(px.parcels)

"""
	 haskey(px::Parcellation{T}, k::T)

Check whether `Parcellation{T} px` contains a parcel with key value `k`
"""
Base.haskey(px::Parcellation) = haskey(px.parcels, k)

"""
    values(px::Parcellation)

Access the `Parcel`s in a `Parcellation`
"""
Base.values(px::Parcellation) = values(px.parcels)

"""
    getindex(px::Parcellation{T}, k::T)

Access a single Parcel within a Parcellation via its key of type `T`"
"""
Base.getindex(px::Parcellation{T}, k::T) where T = px.parcels[k]

"""
    vec(px::Parcellation)

Convert a `Parcellation` from its internal `Dict`-based representation into a `Vector{T}`.
`T` must have a `zeros(T, ...)` method. Warning: this is not a sensible representation
in the event that any `Parcel`s overlap.
"""
function Base.vec(px::Parcellation{T}) where T <: Real
	out = zeros(T, length(px))
	for k in keys(px)
		@inbounds out[vertices(px[k])] .= k
	end
	return out
end

function Base.union(px::Parcellation)
	out = falses(length(px))
	for k in keys(px)
		out .|= px[k].membership
	end
	return out
end

"""
    unassigned(px::Parcellation)

Get a `BitVector` identifying unassigned vertices (`1`) in a `Parcellation`
"""
unassigned(px::Parcellation) = .!union(px)

"""
    nnz(px::Parcellation)

Get the number of vertices within a `Parcellation` that are assigned
to at least one `Parcel`
"""
nnz(px::Parcellation) = sum(union(px))

"""
    density(px::Parcellation)

Get the proportion of assigned parcel vertices of a `Parcellation`
relative to the total number of vertices in its surface representation
"""
density(px::Parcellation) = nnz(px) / length(px)

surface(px::Parcellation) = px.surface

