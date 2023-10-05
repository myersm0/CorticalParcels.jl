
export dilate!, erode!, close!, resize!
export dilate, erode, interstices

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

function dilate(p::Parcel, args...)
	p′ = deepcopy(p)
	dilate!(p′, args...)
	return p′
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
		[any([!(x in verts) for x in neighbors[x]]) for x in verts]
	]
	if !isnothing(limit) && length(border_verts) > limit
		border_verts = border_verts[1:limit]
	end
	setdiff!(p, border_verts)
	return length(border_verts)
end

function erode(p::Parcel, args...)
	p′ = deepcopy(p)
	erode!(p′, args...)
	return p′
end

"""
    close!(p, neighbors)

Given a `Parcel` `p` and an adjacency list `neighbors`, perform a morphological
closing to fill in gaps, if any, by finding vertices in `p` where all of its 
neighbors but one belong to `p`. Note: for performance reasons, this may not be
technically quite the same as a true closing operation, `erode!(dilate!(p))`
"""
function close!(p::Parcel, neighbors::Vector{Vector{Int}})
	candidates = union([neighbors[v] for v in vertices(p)]...)
	while true
		add_inds = filter(x -> sum(.!p[neighbors[x]]) .<= 2, candidates)
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
		p::Parcel, desired_size::Int; A::AbstractMatrix, neighbors::Vector{Vector{Int}}
	)
	curr_size = size(p)
	Δ = curr_size - desired_size
	while Δ != 0
		if Δ < 0
			nchanged = dilate!(p, A; limit = abs(Δ))
			Δ += nchanged
		else
			nchanged = erode!(p, neighbors; limit = Δ)
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

"""
    interstices(p1, p2, A)

Find the boundary vertices lying between two parcels
"""
function interstices(p1::Parcel, p2::Parcel, A::SparseMatrixCSC)
	p1′ = dilate(p1, A)
	p2′ = dilate(p2, A)
	return intersect(p1′, p2′)
end

"""
    interstices(px, A)

Iterate through a parcellation and find, for each pair of neighboring parcels
separated by a 1-vertex-wide gap, vertices in the interstitial region that separates them
"""
function interstices(px::Parcellation{T}, A::SparseMatrixCSC) where T
	v = vec(px)
	u = unassigned(px)
	temp = @view A[u, :]

	# find all unassigned vertices that have 2 or more parcels as neighbors
	status = filter(
		x -> length(x) >= 2,
		ThreadsX.map(x -> sort(setdiff(v[x], 0)), eachrow(temp))
	)

	# from all unique parcel-parcel pairs discovered from the above,
	# make a dict in which to store their interstitial vertices
	result = Dict{Tuple{T, T}, Vector{Int}}()
	for parcel_list in status
		for x in parcel_list
			for y in setdiff(parcel_list, x)
				a = min(x, y)
				b = max(x, y)
				haskey(result, (a, b)) && continue
				result[(a, b)] = interstices(px[a], px[b], A)
			end
		end
	end
	return result
end



