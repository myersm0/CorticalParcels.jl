
function fill_in_gaps!(p::Parcel, neigh::Vector{Vector{Int}})
	candidates = union([neigh[v] for v in vertices(p)]...)
	while true
		add_inds = filter(x -> sum(.!p[neigh[x]]) .<= 1, candidates)
		p2 = Parcel(add_inds; n = length(p))
		sum(setdiff(p2, p)) > 0 || break
		intersect!(p, p2)
	end
end

function dilate!(
		p::Parcel, A::SparseMatrixCSC; limit::Union{Nothing, Int} = nothing
	)
	border_verts = setdiff(
		unique(A[:, setdiff(p, baddata)].rowval),
		vertices(p)
	)
	length(border_verts) > 0 || return
	if !isnothing(limit) && length(border_verts) > limit
		border_verts = border_verts[1:limit]
	end
	border = Parcel(border_verts; n = length(p))
	union!(p, border)
	return length(border_verts)
end

function erode!(
		p::Parcel, neigh::Vector{Vector{Int}}; limit::Union{Nothing, Int} = nothing
	)
	good_verts = findall(setdiff(p, baddata))
	border_verts = good_verts[
		[any([!(x in good_verts) for x in neigh[x]]) for x in good_verts]
	]
	if !isnothing(limit) && length(border_verts) > limit
		border_verts = border_verts[1:limit]
	end
	border = Parcel(border_verts; n = length(p))
	setdiff!(p, border)
	return length(border_verts)
end

function resize!(
		p::Parcel, desired_size::Int, A::AbstractMatrix; maxiter::Int = 20
	)
	curr_size = size(p)
	Δ = curr_size - desired_size
	i = 1
	while Δ != 0 && i < maxiter 
		if Δ < 0
			n_added = dilate!(p, A; limit = abs(Δ))
			Δ += n_added
			i += 1
		else
			n_removed = erode!(p, neigh; limit = Δ)
			Δ -= n_removed
			i += 1
		end
		i == maxiter && println("Warning: maxiter condition reached")
	end
	temp = setdiff(rotated_parcel, baddata)
	border_verts = setdiff(A[:, temp].rowval, temp)
	if length(border_verts) / length(temp) >= 3
		p .= false
	end
end

