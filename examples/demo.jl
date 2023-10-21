
using CorticalSurfaces
using CorticalParcels
using CIFTI
using JLD

# first we need to set up a surface to use (refer to CorticalSurfaces.jl for details)
data_dir = joinpath(dirname(@__FILE__), "..", "data")
surf_file = joinpath(data_dir, "MSC01.jld")
temp = load(surf_file)
coords = temp["pointsets"]["midthickness"][L]
mw = temp["medial wall"][L]
triangle = temp["triangle"][L]
hem = Hemisphere(coords, mw; triangles = triangle)
hem[:neighbors] = make_adjacency_list(hem)
hem[:A] = make_adjacency_matrix(hem)

# now make a "parcel" with just a single vertex, 17344
p1 = Parcel(hem, [17344])

# make another parcel at vertex 8423 (which happens to be 30mm from the previous one)
p2 = Parcel(hem, [8423])

# grow the first parcel a couple of times ...
dilate!(p1) # 6 vertices are added
@assert size(p1) == 7
dilate!(p1) # 12 vertices are added
@assert size(p1) == 19

# the parcels are still distant enough that there's not yet a small margin
# or interstice separating them:
@assert sum(interstices(p1, p2)) == 0

# but we can grow one of the parcels until this is the case,
# i.e. until they almost touch:
while sum(interstices(p1, p2)) == 0
	dilate!(p2)
end

# they still don't overlap yet ...
@assert overlap(p1, p2) == 0
@assert complement(p1, p2) == size(p1)
@assert complement(p2, p1) == size(p2)

# but now there's just thin margin or interstice, 4 vertices long, between them:
margin_vertices = findall(interstices(p1, p2))
@assert length(margin_vertices) == 4

# now make a parcellation from *copies* of these two parcels
px = Parcellation{Int}(hem)
px[1] = Parcel(p1)
px[2] = Parcel(p2)

@assert size(px) == 2 # two parcels

# now combine both parcels from px; note that the latter will be deleted
merge!(px, 1, 2)
@assert size(px) == 1 # just one parcel remains now

# in the merge operation above, parcel 1 has consumed all the vertices of parcel 2, 
# plus the 4 margin/interstitial vertices that separated them
@assert size(px[1]) == size(p1) + size(p2) + length(margin_vertices) 

# now reverse those operations and go back to the way it was a minute ago ...
setdiff!(px[1], margin_vertices)
setdiff!(px[1], p2)
@assert overlap(px[1], p1) == size(px[1]) # px[1] is now identical with p1 again

# add parcel #2 back in
px[2] = Parcel(p2)

# add just one vertex from the interstices so that the two parcels become contiguous
append!(px[1], margin_vertices[1])

# make a new parcel p3 just for demo purposes
p3 = Parcel(px[1])
union!(p3, p2)
@assert size(p3) == 1 + size(p1) + size(p2)

# now p3 has all the vertices from p1, all the vertices from p2,
# plus one vertex linking those two regions; we can cut that vertex
# (an articulation point or cut vertex in graph theory terms) and then
# get back the two resulting connected components, i.e. recover the 
# two original parcels p1 and p2:
orig_parcels = cut(p3)
@assert isequal(orig_parcels[1], p2)
@assert overlap(orig_parcels[2], p1) == size(p1) - 1

# load in a real parcellation form a CIFTI file:
parcel_file = joinpath(data_dir, "test_parcels.dtseries.nii")
temp = CIFTI.load(parcel_file)
px = Parcellation{Int}(hem, temp[L]) # just use left hem for demo

# a few miscellaneous functions:
keys(px)       # get the parcel IDs (of type T) from Parcellation{T} px
vec(px)        # turn a Parcellation{T} into a Vector{T}
unassigned(px) # get a BitVector representing the unassigned vertices in px
union(px)      # collapse all Parcels within px to a single BitVector
nnz(px)        # the number of vertices in px that have parcel membership
length(px)     # the length of a px's vector space representation
density(px)    # proportion of assigned verices: nnz(px) / length(px)



