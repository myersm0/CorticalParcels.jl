
struct Parcel
	surface::Hemisphere
	membership::BitVector
end

abstract type AbstractParcellation end

struct HemisphericParcellation{T} <: AbstractParcellation
	surface::Hemisphere
	parcels::Dict{T, Parcel}
end

struct BilateralParcellation{T} <: AbstractParcellation
	surface::CorticalSurface
	parcels::Dict{BrainStructure, HemisphericParcellation{T}}
end

