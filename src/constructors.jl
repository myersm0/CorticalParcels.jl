
# ===== Parcel constructors =====

"""
    Parcel(surface::Hemisphere, args...)

Make an empty `Parcel` where `surface` dictates the length of the representational space
"""
function Parcel(surface::Hemisphere)
	return Parcel(surface, falses(size(surface)))
end

"""
	 Parcel(hem, verts)

Make a `Parcel`, given its vertex indices
"""
function Parcel(surface::Hemisphere, verts::Vector{Int})
	membership = falses(size(surface))
	membership[verts] .= true
	return Parcel(surface, membership)
end

"""
    Parcel(hem, coords, tree)

Given a `Matrix` of arbitrary x, y, z coordinates and a `KDTree` representing the 
positions of defined cortical vertex indices, make a `Parcel` by mapping those 
coordinates to the set of defined indices via nearest neighbor search
"""
function Parcel(surface::Hemisphere, coords::AbstractMatrix, tree::KDTree)
	inds, dists = knn(tree, coords, 1)
	inds = [x[1] for x in inds] # flatten result to just a vector
	nverts = size(surface)
	return all(inds .> 0) ? Parcel(surface, inds) : Parcel(surface)
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
    BilateralParcellation{T}(surface::Hemisphere, x::AbstractVector)

Create a `BilateralParcellation` from a `Vector` `x`, the length of which should match 
the size of the surface space being supplied. The distinct elements of that `Vector` 
will become the `Parcels` of the resulting `Parcellation` struct. Parcels will be keyed 
by IDs of type `T`; therefore the eltype of the `Vector` you supply must be
coercable to `T`.
"""
function BilateralParcellation{T}(surface::CorticalSurface, x::AbstractVector) where T
	nverts = size(surface)
	input_size = length(x)
	if input_size != nverts
		input_size == size(surface, Exclusive()) || error(DimensionMismatch)
		x = pad(x, surface)
	end
	parcels = Dict{BrainStructure, HemisphericParcellation{T}}()
	for hem in LR
		verts = vertices(surface[hem], Bilateral(), Inclusive())
		parcels[hem] = HemisphericParcellation{T}(surface[hem], x[verts])
	end
	length(intersect(keys(parcels[L]), keys(parcels[R]))) == 0 ||
		error("Found parcels with membership spanning hemispheres; this is not supported")
	return BilateralParcellation{T}(surface, parcels)
end

function HemisphericParcellation{T}(surface::Hemisphere, x::AbstractVector) where T
	nverts = size(surface)
	input_size = length(x)
	if input_size != nverts
		input_size == size(surface, Exclusive()) || error(DimensionMismatch)
		x = pad(x, surface)
	end
	return HemisphericParcellation{T}(
		surface,
		Dict(p => Parcel(surface, findall(x .== p)) for p in setdiff(x, 0))
	)
end

"""
    BilateralParcellation{T}(surface::CorticalSurface, x::AbstractMatrix)

Create a `BilateralParcellation` from a single-column `Matrix` `x`
"""
function BilateralParcellation{T}(surface::CorticalSurface, x::AbstractMatrix) where T
	size(x, 2) == 1 || error("For matrix input, column dimension must be singleton")
	BilateralParcellation{T}(surface, x[:])
end

"""
    HemisphericParcellation{T}(surface::Hemisphere, x::AbstractMatrix)

Create a `HemisphericParcellation` from a single-column `Matrix` `x`
"""
function HemisphericParcellation{T}(surface::Hemisphere, x::AbstractMatrix) where T
	size(x, 2) == 1 || error("For matrix input, column dimension must be singleton")
	HemisphericParcellation{T}(surface, x[:])
end

"""
    BilateralParcellation{T}(surface::CorticalSurface)

Create an empty `BilateralParcellation`
"""
function BilateralParcellation{T}(surface::CorticalSurface) where T
	return BilateralParcellation{T}(
		surface, 
		Dict(hem => HemisphericParcellation{T}(surface[hem]) for hem in LR)
	)
end

"""
    HemisphericParcellation{T}(surface::Hemisphere)

Create an empty `HemisphericParcellation`
"""
function HemisphericParcellation{T}(surface::Hemisphere) where T
	return HemisphericParcellation{T}(surface, Dict{T, Parcel}())
end

