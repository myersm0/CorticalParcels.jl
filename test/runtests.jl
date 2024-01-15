using CorticalParcels
using CorticalSurfaces
using JLD
using Test
using Chain
using CIFTI
using Pkg.Artifacts

data_dir = artifact"CIFTI_test_files"
surf_file = joinpath(data_dir, "MSC01.jld")
temp = load(surf_file)
hems = Dict()
for hem in LR
	coords = temp["pointsets"]["midthickness"][hem]
	mw = temp["medial wall"][hem]
	triangles = temp["triangle"][hem] # required for adjacency calculations below
	hems[hem] = Hemisphere(hem, coords, mw; triangles = triangles)
	make_adjacency_list!(hems[hem])
	make_adjacency_matrix!(hems[hem])
end

c = CorticalSurface(hems[L], hems[R])
hem = c[L] # for most of these tests below we'll just deal with left hem for now

parcel_file = joinpath(data_dir, "test_parcels.dtseries.nii")
cifti_data = CIFTI.load(parcel_file)

types_to_test = [UInt16, Int64]

@testset "CorticalParcels.jl" begin
	for dtype in types_to_test
		px = HemisphericParcellation{dtype}(hem, cifti_data[L])
		@test size(px) == length(setdiff(cifti_data[L], 0))
		@test length(px) == 32492
		@test all(trim(vec(px), hem) .== cifti_data[L])
		parcel_sizes = [size(px[p]) for p in keys(px)]
		@test sum(parcel_sizes) == sum(cifti_data[L] .!= 0)
		parcel_vertices = [vertices(px[p]) for p in keys(px)]
		@test all(length.(parcel_vertices) == parcel_sizes)
		px′ = deepcopy(px)
		@test px′.surface === px.surface
		@test px′ == px
	end

	# things should work with arbitrary non-numeric types of T as well; test this
	dtype = String
	inds = [9, 99, 999, 9999]
	temp = fill("unassigned", 32492)
	temp[inds] .= "test"
	px = HemisphericParcellation{dtype}(hem, temp)
	@test size(px["test"]) == 4
	@test size(px["unassigned"]) == 32492 - 4
	# `vec(px)` is not possible however because we can only do this where T <: Real:
	@test_throws MethodError vec(px)
end

@testset "demo.jl tests" begin
	p1 = Parcel(hem, [17344])
	p2 = Parcel(hem, [8423])
	@test dilate!(p1) == 6
	@test dilate!(p1) == 12

	p1′ = deepcopy(p1)
	dilate!(p1′)
	erode!(p1′)
	@test isequal(p1, p1′)
	@test size(Parcel(hem)) == 0
	
	test_sizes = [1, 100, 500, 1000, 10000, size(hem)]
	for siz in test_sizes
		resize!(p1′, siz)
		@test size(p1′) == siz
	end

	clear!(p1′)
	@test size(p1′) == 0

	p1′ = deepcopy(p1)
	limits = [0, 1, 3, 5, 10, 20, 30]
	limits = [limits; reverse(limits)]
	for limit in limits
		dilate!(p1′; limit = limit)
	end
	@test size(p1′) == size(p1) + sum(limits)

	border_verts = borders(p1)
	@test findall(border_verts) == [17259, 17260, 17261, 17301, 17304, 17342, 17346, 17383, 17386, 17423, 17424, 17425]

	while sum(interstices(p1, p2)) == 0
		dilate!(p2)
	end
	@test overlap(p1, p2) == 0
	@test complement(p1, p2) == size(p1)
	@test complement(p2, p1) == size(p2)

	margin_vertices = findall(interstices(p1, p2))
	@test length(margin_vertices) == 3

	px = HemisphericParcellation{Int}(hem)
	px[1] = deepcopy(p1)
	px[2] = deepcopy(p2)
	@test size(px) == 2

	@test_throws ErrorException px[2] = Parcel(Hemisphere(L, 9999))

	merge!(px, 1, 2)
	@test size(px) == 1 # just one parcel remains now
	@test size(px[1]) == size(p1) + size(p2) + length(margin_vertices) 

	setdiff!(px[1], p2)
	setdiff!(px[1], margin_vertices)
	@test isequal(px[1], p1)

	px[2] = deepcopy(p2)
	append!(px[1], margin_vertices[1])

	p3 = deepcopy(px[1])
	union!(p3, p2)
	@test size(p3) == 1 + size(p1) + size(p2)

	orig_parcels = cut(p3)
	@test isequal(orig_parcels[1], p2)
	@test overlap(orig_parcels[2], p1) == size(p1) - 1

	delete!(px, 1)
	delete!(px, 2)
	@test size(px) == 0

	# load in a real parcellation form a CIFTI file:
	parcel_file = joinpath(data_dir, "test_parcels.dtseries.nii")
	temp = CIFTI.load(parcel_file)
	px = BilateralParcellation{Int}(c, temp)
	pxL = HemisphericParcellation{Int}(c[L], temp[L])
	pxR = HemisphericParcellation{Int}(c[R], temp[R])
	@test vec(px) == vcat(vec(pxL), vec(pxR)) == pad(temp[LR][:], c)

	@test length(keys(px[L])) == 185
	@test density(px[L]) ≈ 0.740613073987443
	@test sum(unassigned(px[L])) == 8428
	@test sum(nnz(px[L])) == length(px[L]) - sum(unassigned(px[L])) == sum(union(px[L]))
end

