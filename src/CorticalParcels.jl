
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

include("types.jl")
include("constructors.jl")
include("accessors_mutators.jl")
include("set_ops.jl")
include("morphology.jl")
include("editing.jl")
include("distances.jl")
include("show.jl")

end

