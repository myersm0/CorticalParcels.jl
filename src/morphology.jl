
function fill_in_gaps!(p::Parcel, neigh::Vector{Vector})
	candidates = # restricted list of neigh
	while true
		add_inds = filter(x -> x .<= 1, [sum(.!p[x]) for x in candidates])
		length(add_inds) > 0 || break
		p[add_inds] .= true
	end
end

function dilate!(p::Parcel, A::SparseMatrixCSC; limit::Int = Inf)
	temp = setdiff(p, baddata)
	border_verts = setdiff(A[:, temp].rowval, temp)
	border_verts = border_verts[1:min(limit, length(border_verts))]
	union!(p, border_verts)
	return length(border_verts)
end

function erode!(p::Parcel, A::SparseMatrixCSC; limit::Int = Inf)
	temp =  setdiff(p, baddata)
	temp2 = findall(temp)
	border_verts = temp2[A[temp, .!temp].rowval]
	border_verts = border_verts[1:min(limit, length(border_verts))]
	setdiff!(p, border_verts)
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
			n_removed = erode!(p, A; limit = Δ)
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

