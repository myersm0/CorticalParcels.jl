# CorticalParcels
This Julia package supplies a set of tools for conveniently and efficiently working with parcels, or regions of interest, in the context of the surface-space representation of the cerebral cortex. It builds upon the `Hemisphere` and `CorticalSurface` types (with supertype `SurfaceSpace`) from [CorticalSurfaces.jl](https://github.com/myersm0/CorticalSurfaces.jl) and provides the foundation for my implementation of an important parcel-generation and -evaluation method in [WatershedParcellation.jl](https://github.com/myersm0/WatershedParcellation.jl). Although the abstractions and implementations here are new, at its core it's based on MATLAB code developed at Washington University by Tim Laumann and Evan Gordon from their 2016 paper ["Generation and Evaluation of a Cortical Area Parcellation from Resting-State Correlations."](https://pubmed.ncbi.nlm.nih.gov/25316338/).

A `Parcel` is a discrete region of interest on the cortical surface, and in this implementation is stored internally as a `BitVector` of vertices where each element denotes membership (`true` or `false`). The total length of that vector constitutes the surface-space representation in which the parcel is defined. The size of a parcel `size(p::Parcel)` is given as the number of non-zero elements of that vector, i.e. the number of vertices belonging to that parcel. This implementation was chosen to enable very fast performance of common operations such as getting the size, computing overlap with other parcels, dilating and eroding, etc, by reducing them internally to simple bitwise operations.

A `Parcellation` is a collection of `Parcel`s that all share the same space. It's typically the case that the parcels within it are non-overlapping. The struct contains two fields:
- `surface`: a `SurfaceSpace` supplying details of the geometry (particularly, the size of the space) that all its component `Parcel`s must conform to
  - (if however the geometry is not of interest in your application, then a dummy surface can be created by, for example, `Hemisphere(32492)` where the only piece of information that's strictly required is the number of vertices, 32492 in this case)
- `parcels`: a `Dict{T, Parcel}` mapping keys of type `T` to parcels, where `T` can be any type (preferably one for which a `zero(::T)` method is defined) that you want to use as keys for accessing and labeling individual parcels

Rather than having to create the `Dict{T, Parcel}` yourself, I anticipate that a `Parcellation` will most often be initialized via its `Parcellation(surface::SurfaceSpace, x::Vector{T})` constructor, since the `Vector{T}` representation is probably the way you read in an existing parcellation from disk, e.g. from a [CIFTI](https://github.com/myersm0/CIFTI.jl") file. See the Usage section below.

A `Parcellation` can be mapped back to a vanilla `Vector{T}` representation if desired via `vec(px::Parcellation)`.

Some notation notes: in the following documentation and in demos, `p`, `p1`, `p2` will refer to individual parcels; and `px` will refer to a whole parcellation.

## Performance
The performance is going to depend on several factors. The benchmarks below are based on using a single-hemisphere parcellation of 185 parcels, in a space of 32492 vertices, and compares the current `BitVector`-based implementation to an alternative using `SparseVector`s as well as to a naive `Vector{Int}` representation (simply a list of vertex index numbers).
- *Adding or removing vertices to/from a `Parcel`*. This is where the current implementation shines most, via operations like `union!(p1::Parcel, p2::Parcel)` and analogous calls to `setdiff!` and `intersect!`.
- *Computing the amount of overlap of two `Parcel`s*. This is fast because it reduces to just taking the dot product of their respective membership vectors.
- *Checking the size of a `Parcel`.* This is the only case where the current implementation lags behind alternatives.
- *Checking a `Parcellation` for unassigned values*. This is relatively "slow" compared to `Parcel`-level operations supplied. But it should be infrequent enough that it doesn't matter much; and the present `BitVector` is still faster than alternatives.

|              |`intersect!(p1, p2)`|`overlap(p1, p2`|`size(p)`|`unassigned(px)`|
|:-------------|-------------------:|-------------------:|-------------------:|-------------------:|
|**`BitVector`**|**85 ns**|**108 ns**|104 ns|**22000** ns|
|`SparseVector`|3047 ns|812 ns|83 ns|39000 ns|
|`Vector`|7692 ns|49110 ns|**9 ns**|1024000 ns|

While the need to compute the size of a parcel is indeed a common operation and we'd like it to be as fast as possible, this implementation's considerable advantage in the other basic operations should still make it the clear frontrunner in most use cases.

If we assume for simplicity that the above operations occur equally often, the `SparseVector` implementation (used in this package version 0.1.0 only) achieves a 25x speedup relative to the naive case, and the present `BitVector` implementation (package version 0.2+) achieves a 48x speedup relative to the same. If we ignore the `unassigned(px)` call, the current implementation improves to a 191x speedup over baseline.

## Installation
Within Julia:
```
using Pkg
Pkg.add("CorticalParcels")
```

## Usage
A full demo of the basic functionality can be found in `examples/demo.jl`.

The package [CorticalSurfaces.jl](https://github.com/myersm0/CorticalSurfaces.jl) provides the definitions of `Hemisphere` and `CorticalSurface` types (and their supertype `SurfaceSpace`), on which many of the operations in this package depend. So first of all, load both packages and create a `Hemisphere` struct that will define the vertex space. At a minimum, you need to specify the number of vertices in that space, for example 32492; but see [CorticalSurfaces.jl](https://github.com/myersm0/CorticalSurfaces.jl) for further details.
```
using CorticalSurfaces
using CorticalParcels

hem = Hemisphere(32492) # create a Hemisphere of 32492 vertices that will define the space
```

### Constructors
Once a `Hemisphere` has been created (we'll call it `hem` here), the following are two basic ways in which to initialize a `Parcel`::
```
Parcel(hem)                  # create an empty parcel within the same space as `hem`
Parcel(hem, [5, 6, 7])       # create a parcel with 3 vertices within the same space as `hem`
```

A `Parcellation` can be initialized in several ways, such as:
```
hem = Hemisphere(32492) # create a Hemisphere of 39492 vertices that will define the space
Parcellation{Int}(hem)  # create an empty parcellation within that space

# as above, but this time fill the space with 10 randomly assigned parcels
Parcellation{Int}(hem, rand(1:10, 32492))
```

The above examples use `Int` as the initialization parameter, and this defines the type of key that will be assigned to each parcel. Any type should be usable, however, provided that its `typemax` can represent the largest value you anticipate needing to store. You could use `String` keys, for example, if you want to provide descriptive labels for your parcels and index them in that way.

### Accessors
`unassigned(px::Parcellation)` may be used to dynamically determine the elements in the vector space that are not assigned to any parcel.

`vec(px::Parcellation)` will reduce the parcellation to a single `Vector{T}`. If you constructed `px` from a `Vector{T}` (and have not changed any of its elements), this operation should return that same vector.

More details and a demo coming soon.

[![Build Status](https://github.com/myersm0/CorticalParcels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/myersm0/CorticalParcels.jl/actions/workflows/CI.yml?query=branch%3Amain)
