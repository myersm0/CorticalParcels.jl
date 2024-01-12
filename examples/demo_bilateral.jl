
using CorticalSurfaces
using CorticalParcels
using CIFTI
using JLD

# first we need to set up a surface to use (refer to CorticalSurfaces.jl for details);
# here we'll just do one hemisphere, the left
data_dir = joinpath(dirname(@__FILE__), "..", "data")
surf_file = joinpath(data_dir, "MSC01.jld")
temp = load(surf_file)
hems = Dict()
for hem in LR
	coords = temp["pointsets"]["midthickness"][hem]
	mw = temp["medial wall"][hem]
	triangle = temp["triangle"][hem] # required for adjacency calculations below
	hems[hem] = Hemisphere(coords, mw; triangles = triangle) 
	hems[hem][:neighbors] = make_adjacency_list(hems[hem])
	hems[hem][:A] = make_adjacency_matrix(hems[hem])
end

c = CorticalSurface(hems[L], hems[R])
	

# note that the two adjacency items created above are required for many parcelwise ops
# (particularly ones that require graph traversal such as erode!() and dilate!()),
# though you could still do some things without them

# given the vertex space defined in the left Hemisphere struct above,
# now make a "parcel" within that space consisting of just a single vertex, 17344
p1 = Parcel(c[L], [17344])

# make another parcel at vertex 8423 (which happens to be 30mm from the previous one)
p2 = Parcel(c[L], [8423])

# grow the first parcel a couple of times, and check the size afterwards ...
dilate!(p1) # 6 vertices are added, so size is now 7
@assert size(p1) == 7
dilate!(p1) # 12 vertices are added, so size is now 19
@assert size(p1) == 19

# make a copy of p1 and do some various resizing
p1′ = Parcel(p1)
dilate!(p1′)
erode!(p1′)
@assert isequal(p1, p1′)

# resize to an arbitrary size, say 500 vertices, by repeated dilation:
resize!(p1′, 500)

# or shrink (erode) it to 100 vertices:
resize!(p1′, 100)

# dilate it once more, but don't add more than 10 new vertices:
dilate!(p1′; limit = 10)
@assert size(p1′) <= 110

# if you want to see which vertices belong to a parcel:
vertices(p1′)

# remove all vertices from the parcel
clear!(p1′)

# grow p2 iteratively until there's only a small margin or interstice 
# separating it from p1:
while sum(interstices(p1, p2)) == 0
	n_new_verts = dilate!(p2)
	println("Added $n_new_verts vertices to p2 (total size: $(size(p2)))")
end

# they still don't quite overlap yet ...
@assert overlap(p1, p2) == 0
@assert complement(p1, p2) == size(p1)
@assert complement(p2, p1) == size(p2)

# but there's only a thin margin or interstice, 4 vertices long, between them:
margin_vertices = findall(interstices(p1, p2))
@assert length(margin_vertices) == 3

# now make an empty parcellation within the space of our left Hemisphere struct,
# using keys (parcel IDs) of type Int:
px = HemisphericParcellation{Int}(c[L])

# give it *copies* of the two parcels we were working with above
px[1] = Parcel(p1)
px[2] = Parcel(p2)

@assert size(px) == 2 # two parcels

# now combine parcels 1 and 2 from px; parcel 2 will be incorporated into 
# parcel 1, along with the 3 interstitial vertices in between, and then deleted
merge!(px, 1, 2)
@assert size(px) == 1 # just one parcel remains now
@assert size(px[1]) == size(p1) + size(p2) + length(margin_vertices)

# now reverse those operations and go back to the way it was a minute ago
setdiff!(px[1], p2)
setdiff!(px[1], margin_vertices)
@assert isequal(px[1], p1)

# add a copy of parcel #2 back in
px[2] = Parcel(p2)

# add just one vertex from the interstices so that the two parcels become connected
append!(px[1], margin_vertices[1])

# make a new parcel p3 just for demo purposes
p3 = Parcel(px[1])
union!(p3, p2) # combine p3 and p2
@assert size(p3) == 1 + size(p1) + size(p2)

# now p3 has all the vertices from p1, all the vertices from p2,
# plus one vertex linking those two regions; we can cut that vertex
# (an articulation point or cut vertex in graph theory terms) and then
# get back the two resulting connected components, i.e. recover the 
# two original parcels p1 and p2 though not necessarily in the same order:
orig_parcels = cut(p3)
@assert isequal(orig_parcels[1], p2)
@assert overlap(orig_parcels[2], p1) == size(p1) - 1

# load in a real parcellation form a CIFTI file:
parcel_file = joinpath(data_dir, "test_parcels.dtseries.nii")
temp = CIFTI.load(parcel_file)
px = HemisphericParcellation{Int}(c[L], temp[L]) # just use left hem for demo

# every time you show px, it will display properties of a few random parcels
px
px
px

# some miscellaneous functions:
keys(px)       # get the parcel IDs (of type T) from Parcellation{T} px
vec(px)        # turn a Parcellation{T} into a Vector{T}
unassigned(px) # get a BitVector representing the unassigned vertices in px
union(px)      # collapse all Parcels within px to a single BitVector
nnz(px)        # the number of vertices in px that have parcel membership
length(px)     # the length of px's vector space representation
density(px)    # proportion of assigned verices: nnz(px) / length(px)



