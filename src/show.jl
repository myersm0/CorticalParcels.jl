
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

function Base.show(io::IO, mime::MIME"text/plain", px::HemisphericParcellation)
	ks = @chain keys(px) collect sample(_, size(px); replace = false)
	dens = Int(round(density(px) * 100; digits = 0))
	print(io, "HemisphericParcellation{$(eltype(keys(px)))} with $(size(px)) parcels,")
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

function Base.show(io::IO, mime::MIME"text/plain", px::BilateralParcellation)
	println("CORTEX_LEFT:")
	show(io, mime, px[L])
	println("CORTEX_RIGHT:")
	show(io, mime, px[R])
end

function Base.show(
		io::IO, mime::MIME"text/plain", pxs::Vector{HemisphericParcellation{T}}
	) where T
	print(io, "Vector of $(length(pxs)) HemisphericParcellation with keys of type $T")
end

function Base.show(
		io::IO, mime::MIME"text/plain", pxs::Vector{BilateralParcellation{T}}
	) where T
	print(io, "Vector of $(length(pxs)) BilateralParcellation with keys of type $T")
end


