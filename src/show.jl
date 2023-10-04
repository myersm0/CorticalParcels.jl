
function Base.show(
		io::IO, ::MIME"text/plain", p::Parcel; label::Any = nothing
	)
	siz = size(p)
	len = length(p)
	dens = round(density(p) * 100; digits = 2)
	id_str = isnothing(label) ? "" : " [$label]"
	print(io, "Parcel")
	printstyled(io, id_str; bold = true)
	print(io, " with $siz non-zero vertices out of $len")
	print(io, " ($dens% dense)")
end

function Base.show(io::IO, mime::MIME"text/plain", px::Parcellation)
	ks = @chain keys(px) collect sample(_, size(px); replace = false)
	dens = Int(round(density(px) * 100; digits = 0))
	print(io, "Parcellation{$(eltype(keys(px)))} with $(size(px)) parcels,")
	print(io, " $dens% dense,")
	print(io, " in a space of $(size(px.surface)) vertices")
	print(io, "\n")
	for i in 1:min(length(ks), 3)
		print(io, "    ⊢ ")
		show(io, mime, px[ks[i]]; label = ks[i])
		print(io, "\n")
	end
	if length(ks) > 4
		println(io, "    ⊢ ⋮")
	end
	if length(ks) > 3
		print(io, "    ⊢ ")
		show(io, mime, px[ks[end]]; label = ks[end])
	end
end

function Base.show(io::IO, mime::MIME"text/plain", pxs::Vector{Parcellation{T}}) where T
	print(io, "Vector of $(length(pxs)) Parcellations with keys of type $T")
end


