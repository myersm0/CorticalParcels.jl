
export Parcellation, size, length, keys, values, getindex, vec

struct Parcellation{T}
	surface::SurfaceSpace
	parcels::Dict{T, Parcel}
end

"""
    Parcellation{T}(surface::SurfaceSpace, x::AbstractVector)

Create a `Parcellation` from a `Vector` `x`, the length of which should match the size
of the surface space being supplied. The distinct elements of that `Vector` will 
become the `Parcels` of the resulting `Parcellation` struct. Parcels will be keyed 
by IDs of type `T`; therefore the eltype of the `Vector` you supply must be
coercable to `T`.
"""
function Parcellation{T}(surface::SurfaceSpace, x::AbstractVector) where T
	nverts = size(surface)
	input_size = length(x)
	if input_size != nverts
		input_size == size(surface, Exclusive()) ||
			error("Input Vector doesn't match size of the supplied surface")
		x = pad(x, surface)
	end
	return Parcellation{T}(
		surface,
		Dict(
			[p => Parcel(findall(x .== p); n = nverts) for p in setdiff(x, 0)]
		)
	)
end

"""
    Parcellation{T}(surface::SurfaceSpace, x::AbstractMatrix)

Create a `Parcellation` from a single-column `Matrix` `x`
"""
function Parcellation{T}(surface::SurfaceSpace, x::AbstractMatrix) where T
	size(x, 2) == 1 || error("For Matrix input, column dimension must be singleton")
	Parcellation{T}(surface, x[:])
end

"""
    Parcellation{T}(surface::SurfaceSpace)

Create an empty `Parcellation`
"""
function Parcellation{T}(surface::SurfaceSpace) where T
	Parcellation{T}(surface, Dict{T, Parcel}())
end

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

unassigned(px::Parcellation) = findall(vec(px) .== 0)
nnz(px::Parcellation) = sum(vec(px) .!= 0)

function Base.show(io::IO, ::MIME"text/plain", px::Parcellation)
	print(io, "Parcellation{$(eltype(keys(px)))} with $(size(px)) parcels,")
	print(io, " in a space of $(size(px.surface)) vertices")
	print(io, " ($(size(px.surface, Exclusive())) without medial wall)")
end


