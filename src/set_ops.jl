
export intersect, union, setdiff, intersect!, union!, setdiff!
export overlap, complement

Base.intersect(p1::Parcel, p2::Parcel) = p1.membership .& p2.membership
Base.union(p1::Parcel, p2::Parcel) = p1.membership .| p2.membership
Base.setdiff(p1::Parcel, p2::Parcel) = p1.membership .& .!p2.membership

Base.intersect!(p1::Parcel, p2::Parcel) = p1.membership .&= p2.membership
Base.union!(p1::Parcel, p2::Parcel) = p1.membership .|= p2.membership
Base.setdiff!(p1::Parcel, p2::Parcel) = p1.membership .&= .!p2.membership

Base.intersect(p::Parcel, x::BitVector) = p.membership &= x
Base.union(p::Parcel, x::BitVector) = p.membership |= x
Base.setdiff(p::Parcel, x::BitVector) = p.membership &= .!x

Base.intersect!(p::Parcel, x::BitVector) = p.membership .&= x
Base.union!(p::Parcel, x::BitVector) = p.membership .|= x
Base.setdiff!(p::Parcel, x::BitVector) = p.membership .&= .!x

Base.intersect!(p::Parcel, x::Vector{T}) where T <: Integer = p.membership .*= x
Base.union!(p::Parcel, x::Vector{T}) where T <: Integer = p.membership[x] .= true
Base.setdiff!(p::Parcel, x::Vector{T}) where T <: Integer = p.membership[x] .= false

"""
    overlap(p1::Parcel, p2::Parcel)

Compute the number of member vertices shared between two `Parcel`s `p1`, `p2`
"""
overlap(p1::Parcel, p2::Parcel) = p1.membership' * p2.membership

overlap(p::Parcel, x::Union{Vector{Bool}, BitVector}) = p.membership' * x

overlap(p::Parcel, px::Parcellation) = return sum([overlap(p, px[k]) for k in keys(px)])

overlap(px::Parcellation, p::Parcel) = return overlap(p, px)

overlap(x::Union{Vector{Bool}, BitVector}, p::Parcel) = overlap(p, x)

function overlap(px::Parcellation)
	counts = zeros(Int, length(px))
	for k in keys(px)
		counts .+= px[k].membership
	end
	overlap_inds = findall(counts .> 1)
	return sum(counts[overlap_inds] .- 1)
end

"""
    complement(p1::Parcel, p2::Parcel)

Compute the number of member vertices in `Parcel` `p1` not shared by those of `p2`
"""
complement(p1::Parcel, p2::Parcel) = p1.membership' * .!p2.membership

complement(p::Parcel, x::Union{Vector{Bool}, BitVector}) = p.membership' * .!x

complement(x::Union{Vector{Bool}, BitVector}, p::Parcel) = complement(p, x)


