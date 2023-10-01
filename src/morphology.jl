
export dilate!, erode!, close!, resize!

"""
    dilate!(p, A; limit = nothing)

Perform a single pass of dilation on `Parcel` `p`, guided by adjacency matrix `A`;
optionally specify a `limit::Int` on the number of new vertices that can be added
"""
function dilate!(
		p::Parcel, A::SparseMatrixCSC; limit::Union{Nothing, Int} = nothing
	)
	parcel_verts = vertices(p)
	border_verts = setdiff(unique(A[:, parcel_verts].rowval), parcel_verts)
	length(border_verts) > 0 || return
	if !isnothing(limit) && length(border_verts) > limit
		border_verts = border_verts[1:limit]
	end
	border = Parcel(border_verts; n = length(p))
	union!(p, border)
	return length(border_verts)
end

"""
    erode!(p, neighbors; limit = nothing)

Perform a single pass of erosion on `Parcel` `p`, guided by adjacency list `neighbors`;
optionally specify a `limit::Int` on the number of vertices that you want to remove
"""
function erode!(
		p::Parcel, neighbors::Vector{Vector{Int}}; limit::Union{Nothing, Int} = nothing
	)
	verts = vertices(p)
	border_verts = verts[
		[any([!(x in verts) for x in neigh[x]]) for x in verts]
	]
	if !isnothing(limit) && length(border_verts) > limit
		border_verts = border_verts[1:limit]
	end
	setdiff!(p, border_verts)
	return length(border_verts)
end

"""
    close!(p, neighbors)

Given a `Parcel` `p` and an adjacency list `neighbors`, perform a morphological
closing to fill in gaps, if any, by finding vertices in `p` where all of its 
neighbors but one belong to `p`. Note: for performance reasons, this may not be
technically quite the same as a true closing operation, `erode!(dilate!(p))`
"""
function close!(p::Parcel, neighbors::Vector{Vector{Int}})
	candidates = union([neigh[v] for v in vertices(p)]...)
	while true
		add_inds = filter(x -> sum(.!p[neigh[x]]) .<= 2, candidates)
		p2 = Parcel(add_inds; n = length(p))
		complement(p2, p) > 0 || break
		union!(p, p2)
	end
end

"""
	 resize!(p, desired_size; A, neighbors)

Resize a `Parcel` `p`, guided by an adjacency matrix and an adjacency list, 
by repeated dilation or erosion until `p` reaches `desired_size`
"""
function Base.resize!(
		p::Parcel, desired_size::Int; A::AbstractMatrix, neigh::Vector{Vector{Int}}
	)
	curr_size = size(p)
	Δ = curr_size - desired_size
	while Δ != 0
		if Δ < 0
			nchanged = dilate!(p, A; limit = abs(Δ))
			Δ += nchanged
		else
			nchanged = erode!(p, neigh; limit = Δ)
			Δ -= nchanged
		end
		if nchanged == 0
			siz = size(p)
			println("Could not achieve size $desired_size; stopped at $siz")
			return siz
		end
	end
	return desired_size
end

