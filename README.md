# CorticalParcels
This Julia package supplies a set of tools for conveniently and efficiently working with parcels, or regions of interest, in the context of the surface-space representation of the cerebral cortex. It builds upon the `Hemisphere` and `CorticalSurface` types (with supertype `SurfaceSpace`) from [CorticalSurfaces.jl](https://github.com/myersm0/CorticalSurfaces.jl) and provides the foundation for my implementation of an important parcel-generation and -evaluation method in [WatershedParcellation.jl](https://github.com/myersm0/WatershedParcellation.jl). Although the abstractions and implementations here are new, at its core it's based on MATLAB code developed at Washington University by Tim Laumann and Evan Gordon from their 2016 paper ["Generation and Evaluation of a Cortical Area Parcellation from Resting-State Correlations."](https://pubmed.ncbi.nlm.nih.gov/25316338/).

A `Parcel` is a discrete region of interest on the cortical surface, and in this implementation is stored internally as a `SparseVector` of vertices where each element denotes membership (`true` or `false`). The total length of that vector constitutes the surface-space representation in which the parcel is defined. The size of a parcel `size(p::Parcel)` is given as the number of non-zero elements of that vector, i.e. the number of vertices belonging to that parcel. This implementation was chosen in order to enable flexibility and fast performance of common operations such as getting the vertex indices, getting the size, computing overlap, dilating and eroding, etc.

A `Parcellation` is a collection of `Parcel`s that all share the same space. It's typical, but not verified or required, that the parcels within it are non-overlapping. The struct contains two fields:
- a `SurfaceSpace` supplying details of the geometry (particularly, the size of the space) that all its component `Parcel`s must conform to
  - (if however the geometry is not of interest in your application, then a dummy surface can be created by, for example, `Hemisphere(32492)` where the only piece of information that's strictly required is the number of vertices, 32492 in this case)
- a `Dict{T, Parcel}` mapping keys of type `T` to parcels, where `T` can be any type that you want to use as keys for accessing and labeling individual parcels

## Installation
Within Julia:
```
using Pkg
Pkg.add(url = "http://github.com/myersm0/CorticalParcels.jl")
```

## Usage

### Constructors
The following are two basic ways in which to initialize a `Parcel`::
```
Parcel(32492)                # create an empty parcel within a space of 32492 vertices
Parcel([1, 2, 3]; n = 32492) # create a parcel with 3 vertices within a space of 32492
```

A `Parcellation` can be initialized in several ways, such as:
```
hem = Hemisphere(32492) # create a Hemisphere of 39492 vertices that will define the space
Parcellation{Int}(hem)  # create an empty parcellation within that space

# as above, but this time fill the space with 10 randomly assigned parcels
Parcellation{Int}(hem, rand(1:10, 32492))
```

The above examples use `Int` as the initialization parameter, and this defines the type of key that will be assigned to each parcel. Any type should be usable, however, provided that its `typemax` can represent the largest value you anticipate needing to represent. You could use `String` keys, for example, if you want to provide descriptive labels for your parcels and index them in that way.

### Accessors
Coming soon.

[![Build Status](https://github.com/myersm0/CorticalParcels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/myersm0/CorticalParcels.jl/actions/workflows/CI.yml?query=branch%3Amain)
