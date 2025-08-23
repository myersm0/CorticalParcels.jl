
module CorticalParcels

using CIFTI
using CorticalSurfaces
using NearestNeighbors
using SparseArrays
using StatsBase: sample

# import some type-aliasing constants for convenience
import CorticalSurfaces: AdjacencyList, AdjacencyMatrix, DistanceMatrix

include("types.jl")
export Parcel, AbstractParcellation, HemisphericParcellation, BilateralParcellation

include("constructors.jl")

include("accessors.jl")
export vertices, size, length, keys, haskey, values, getindex
export vec, union, unassigned, nnz, density

include("set_ops.jl")
export intersect, union, setdiff, intersect!, union!, setdiff!
export overlap, complement

include("morphology.jl")
export dilate!, erode!, close!, resize!
export dilate, erode, interstices, borders

include("editing.jl")
export setindex!, cut, split, clear!, delete!, append!, merge!, deepcopy

include("distances.jl")
export DistanceMethod, CentroidToCentroid, ClosestVertices, centroid, distance

include("show.jl")

end

