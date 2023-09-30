using CorticalParcels
using CorticalSurfaces
using JLD
using Test
using Chain
using CIFTI

data_dir = joinpath(dirname(@__FILE__), "..", "data")

mw_file = joinpath(data_dir, "mw_verts.jld")
medial_wall = falses(32492)
temp = load(mw_file, "mw_verts")
medial_wall[filter(x -> x <= 32492, temp)] .= true
hem = Hemisphere(medial_wall)

parcel_file = joinpath(data_dir, "test_parcels.dtseries.nii")
cifti_data = CIFTI.load(parcel_file)

types_to_test = [UInt16, Int64]

@testset "CorticalParcels.jl" begin
	for dtype in types_to_test
		px = Parcellation{dtype}(hem, cifti_data[L])
		@test size(px) == length(setdiff(cifti_data[L], 0))
		@test length(px) == 32492
		@test all(trim(vec(px), hem) .== cifti_data[L])
		parcel_sizes = [size(px[p]) for p in keys(px)]
		@test sum(parcel_sizes) == sum(cifti_data[L] .!= 0)
		parcel_vertices = [vertices(px[p]) for p in keys(px)]
		@test all(length.(parcel_vertices) == parcel_sizes)
	end

	# things should work with arbitrary non-numeric types of T as well; test this
	dtype = String
	inds = [9, 99, 999, 9999]
	temp = fill("unassigned", 32492)
	temp[inds] .= "test"
	px = Parcellation{dtype}(hem, temp)
	@test size(px["test"]) == 4
	@test size(px["unassigned"]) == 32492 - 4
	# `vec(px)` is not possible however because we can only do this where T <: Real:
	@test_throws MethodError vec(px)

	# test set ops
	dtype = Int
	px = Parcellation{dtype}(hem, cifti_data[L])
	a = deepcopy(px[6779])
	b = px[18124]
	@test sort(union(a, b)) == sort(union(vertices(a), vertices(b)))
	@test setdiff(a, b) == a.membership.nzind
	a[vertices(b)] .= true
	@test sort(intersect(a, b)) == vertices(b)
	a2 = deepcopy(px[6779])
	a[:] .= false
	@test size(a) == 0
end




