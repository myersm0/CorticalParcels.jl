

function Base.getindex(
		c::CiftiStruct{E, CIFTI.BRAIN_MODELS(), C}, p::Parcel
	) where {E, C}
	hem = brainstructure(p.surface)
	n_out = length(dts.brainstructure[hem])
	nverts_excl = size(p.surface, Exclusive())
	nverts_excl == n_out || error("Provided parcel doesn't match CIFTI's vertex space")
	verts = c.brainstructure[hem][collapse(vertices(p), p.surface)]
	return c[verts, :]
end

function Base.getindex(
		c::CiftiStruct{E, CIFTI.BRAIN_MODELS(), CIFTI.BRAIN_MODELS()}, p::Parcel
	) where E
	hem = brainstructure(p.surface)
	n_out = length(dts.brainstructure[hem])
	nverts_excl = size(p.surface, Exclusive())
	nverts_excl == n_out || error("Provided parcel doesn't match CIFTI's vertex space")
	verts = c.brainstructure[hem][collapse(vertices(p), p.surface)]
	return c[verts, verts]
end





