
export Parcel, Parcellation

struct Parcel{T <: Hemisphere}
	surface::T
	membership::BitVector
end

struct Parcellation{T}
	surface::SurfaceSpace
	parcels::Dict{T, Parcel}
end

