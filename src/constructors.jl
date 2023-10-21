
# ===== Parcel constructors =====

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
    Parcel(p)

Create a new `Parcel` that's a copy of another one `p`
"""
function Parcel(p::Parcel)
	Parcel(p.surface, vertices(p))
end


# ===== Parcellation constructors =====

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
			[p => Parcel(surface, findall(x .== p)) for p in setdiff(x, 0)]
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

