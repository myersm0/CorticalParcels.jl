# CorticalParcels
This Julia package supplies a set of tools for conveniently and efficiently working with parcels, or regions of interest, in the context of the surface-space representation of the cerebral cortex. It builds upon the `Hemisphere` and `CorticalSurface` types (with supertype `SurfaceSpace`) from [CorticalSurfaces.jl](https://github.com/myersm0/CorticalSurfaces.jl) and provides the foundation for my implementation of an important parcel-generation and -evaluation method in [WatershedParcellation.jl](https://github.com/myersm0/WatershedParcellation.jl).

Internally, a `Parcel` is stored as a `SparseVector` of vertices where each element denotes membership (`true` or `false`), and the total length of that vector constitutes the surface-space representation in which the parcel is defined. The size of a parcel `size(p::Parcel)` is given as the number of non-zero elements of that vector, i.e. the number of vertices belonging to that parcel.

A `Parcellation` is a collection of `Parcel`s that all share the same space. It's typical, but not verified and not necessarily required, that the parcels within it are non-overlapping. The struct contains two fields:
- a `SurfaceSpace` supplying details of the geometry (particularly, the size of the space) that all its component `Parcel`s must conform to
- a `Dict{T, Parcel}` mapping keys of type `T` to parcels, where `T` can be any type that you want to use as keys for accessing and labeling individual parcels

## Installation
Within Julia:
```
using Pkg
Pkg.add(url = "http://github.com/myersm0/CorticalParcels.jl")
```

## Usage

### Constructors
The following are two basic ways in which to initialize a `Parcel`, i.e. a discrete region of interest on the cortical surface:
```
Parcel(32492)                # create an empty parcel within a space of 32492 vertices
Parcel([1, 2, 3]; n = 32492) # create a parcel with 3 vertices within a space of 32492
```

A `Parcellation` is a collection of parcels, all within the same space, and can be initialized in several ways, such as:
```
hem = Hemisphere(32492) # create a Hemisphere of 39492 vertices that will define the space
Parcellation{Int}(hem)  # create an empty parcellation within that space

# as above, but this time fill the space with 10 randomly assigned parcels
Parcellation{Int}(hem, rand(1:10, 32492))
```

The above examples use `Int` as the initialization parameter, and this defines the type of key that will be assigned to each parcel. Any type should be usable, however, provided that its `typemax` can represent the largest value you anticipate needing to represent. You could use `String` keys, for example, if you want to provide descriptive labels for your parcels and index them in that way.

### Accessors
To get the number of vertices that constitute a parcel:
```
Parcel([1, 2, 3]; n = 32492) # create a parcel with 3 vertices within a space of 32492

[![Build Status](https://github.com/myersm0/CorticalParcels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/myersm0/CorticalParcels.jl/actions/workflows/CI.yml?query=branch%3Amain)
