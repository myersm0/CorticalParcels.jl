
import CorticalSurfaces: vertices

# ===== Parcel functions =====

"""
    membership(p)

Get a `BitVector` denoting vertexwise parcel membership
"""
membership(p::Parcel) = p.membership

"""
    vertices(p)

Get the vertex indices belonging to a `Parcel`. Indices will be numbered
inclusive of medial wall by default.
"""
vertices(p::Parcel) = findall(p.membership)

"""
    size(p)

Get the size (number of non-zero vertices) of a `Parcel`".
"""
Base.size(p::Parcel) = sum(p.membership)

"""
    length(p)

Get the length of the representational space in which a `Parcel` is located.
"""
Base.length(p::Parcel) = length(p.membership)

"""
    density(p)

Get the proportion of member vertices in a `Parcel` relative to the length of its space.
"""
density(p::Parcel) = size(p) / length(p)

Base.getindex(p::Parcel, args...) = getindex(p.membership, args...)

Base.isequal(p1::Parcel, p2::Parcel) = 
	p1.surface == p2.surface && p1.membership == p2.membership

Base.isequal(p1::Parcel, x::BitVector) = 
	p1.membership == p2.membership

Base.:(==)(p1::Parcel, p2::Parcel) = isequal(p1, p2)
Base.:(==)(p1::Parcel, x::BitVector) = isequal(p1, p2)

CorticalSurfaces.brainstructure(p::Parcel) = brainstructure(p.surface)

# ===== Parcellation functions =====

"""
    size(px)

Get the number of Parcels comprising a Parcellation.
"""
Base.size(px::HemisphericParcellation) = length(px.parcels)

Base.size(px::BilateralParcellation) = size(px[L]) + size(px[R])

"""
    length(px)

Get the number of vertices comprising the representation space of a `Parcellation`.
"""
Base.length(px::AbstractParcellation) = size(px.surface)

"""
    keys(px)

Get the IDs of all `Parcel`s within a `Parcellation`.
"""
Base.keys(px::HemisphericParcellation) = keys(px.parcels)

Base.keys(px::BilateralParcellation) = union(keys(px[L]), keys(px[R]))

"""
	 haskey(px, k)

Check whether `Parcellation{T} px` contains a parcel with key value `k`.
"""
Base.haskey(px::HemisphericParcellation{T}, k::T) where T = haskey(px.parcels, k)

Base.haskey(px::BilateralParcellation{T}, k::T) where T = 
	haskey(px[L], k) || haskey(px[R], k)

"""
    values(px)

Access the `Parcel`s in a `Parcellation`.
"""
Base.values(px::HemisphericParcellation) = values(px.parcels)

"""
    getindex(px, k)

Access a single Parcel within a Parcellation via its key of type `T`".
"""
Base.getindex(px::HemisphericParcellation{T}, k::T) where T = px.parcels[k]

Base.getindex(px::BilateralParcellation{T}, k::T) where T = 
	haskey(px[L], k) ? px[L][k] : px[R][k]

Base.getindex(px::BilateralParcellation, b::BrainStructure) = px.parcels[b]

"""
    vec(px)

Convert a `HemisphericParcellation` from its internal `Dict`-based representation into 
a `Vector{T}`. `T` must have a `zeros(T, ...)` method. Warning: this is not a sensible
representation in the event that any `Parcel`s overlap.
"""
function Base.vec(px::HemisphericParcellation{T}) where T <: Real
	out = zeros(T, length(px))
	for k in keys(px)
		@inbounds out[vertices(px[k])] .= k
	end
	return out
end

"""
    vec(px)

Convert a `BilateralParcellation` from its internal `Dict`-based representation into 
a `Vector{T}`. `T` must have a `zeros(T, ...)` method. Warning: this is not a sensible
representation in the event that any `Parcel`s overlap.
"""
function Base.vec(px::BilateralParcellation{T}) where T <: Real
	return vcat(vec(px[L]), vec(px[R]))
end

function Base.union(px::HemisphericParcellation)
	out = falses(length(px))
	for k in keys(px)
		out .|= px[k].membership
	end
	return out
end

"""
    unassigned(px)

Get a `BitVector` identifying unassigned vertices (`1`) in a parcellation.
"""
unassigned(px::HemisphericParcellation) = .!union(px)

unassigned(px::BilateralParcellation) = vcat(unassigned(px[L]), unassigned(px[R]))

"""
    nnz(px)

Get the number of vertices within a parcellation that are assigned
to at least one `Parcel`.
"""
nnz(px::HemisphericParcellation) = sum(union(px))

nnz(px::BilateralParcellation) = nnz(px[L]) + nnz(px[R])

"""
    density(px)

Get the proportion of assigned parcel vertices of a parcellation
relative to the total number of vertices in its surface representation.
"""
density(px::AbstractParcellation) = nnz(px) / length(px)

function Base.:(==)(px1::HemisphericParcellation, px2::HemisphericParcellation)
	px1.surface == px2.surface || return false
	all([haskey(px2, k) && px1[k] == px2[k] for k in keys(px1)]) || return false
	return true
end

function Base.:(==)(px1::BilateralParcellation, px2::BilateralParcellation)
	return px1[L] == px2[L] && px1[R] == px2[R]
end



