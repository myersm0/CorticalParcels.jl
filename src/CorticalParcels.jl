
module CorticalParcels

using CIFTI
using CorticalSurfaces
using Chain
using NearestNeighbors
using HDF5
using SparseArrays
using StatsBase: sample
using ThreadsX

# import some type-aliasing constants for convenience
import CorticalSurfaces: AdjacencyList, AdjacencyMatrix, DistanceMatrix

include("parcels.jl")
include("parcellations.jl")
include("morphology.jl")
include("editing.jl")
include("show.jl")

end

