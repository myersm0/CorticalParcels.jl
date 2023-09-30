
function fill_in_gaps!(p::Parcel, neigh::Vector{Vector})
	# first come up with a prelim list of neighbors for this parcel
	candidates = # restricted list of neigh
	while true
		add_inds = filter(x -> x .<= 1, [sum(.!p[x]) for x in candidates])
		length(add_inds) > 0 || break
		p[add_inds] .= true
	end
end

function resize!(
		p::Parcel, desired_size::Int, A::AbstractMatrix; maxiter::Int = 20
	)
	curr_size = size(p)
	Δ = curr_size - desired_size
	i = 1
	# if rotated parcel is smaller than the real one, grow it:
	while Δ < 0 && i < maxiter 
		# find wherever there's no rot parcel assignment but at least one of the neighbors
		# belongs to the parcel; fill up the rot parcel with as many of these border verts
		# as you can, without making the rot parcel bigger than the real one
		temp = findall(setdiff(p, baddata))
		border_verts = setdiff(A[:, temp].rowval, temp)
		border_verts = border_verts[1:min(abs(Δ), length(border_verts))]
		p[border_verts] .= true
		Δ += length(border_verts)
		i += 1
		i == maxiter && println("Warning: maxiter condition reached")
	end
	# if rotated parcel is bigger than the real one, shrink it:
	while Δ > 0 && i < maxiter 
		# find wherever there's a rotated parcel vert but at least one of its neighbors 
		# doesn't belong to the parcel; get rid of as many of these border verts as you 
		# can, without making the rot parcel smaller than the real one
		temp =  setdiff(p, baddata)
		temp2 = findall(temp)
		border_verts = temp2[A[temp, .!temp].rowval]
		border_verts = border_verts[1:min(Δ, length(border_verts))]
		rotated_parcel[border_verts] .= false
		Δ -= length(border_verts)
		i += 1
		i == maxiter && println("Warning: maxiter condition reached")
	end
	temp = setdiff(rotated_parcel, baddata)
	border_verts = setdiff(A[:, temp].rowval, temp)
	if length(border_verts) / length(temp) >= 3
		p .= false
	end
end

