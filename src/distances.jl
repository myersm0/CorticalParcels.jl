
export DistanceMetric, CentroidToCentroid, ClosestVertices, centroid, distance

abstract type DistanceMetric end
struct CentroidToCentroid <: end
struct ClosestVertices <: end

"""
    centroid(p, distances)

Find the centroid of a parcel (the vertex that has the least summed distance
to all other vertices in the parcel). `distances` is expected to be
a square distance matrix of dimensions (length(p), length(p)).
"""
function centroid(p::Parcel, distances::AbstractMatrix)
	all(size(distances) .== length(p)) || error(DimensionMismatch)
	verts = vertices(p)
	summed_dists = sum(distances[verts, verts]; dims = 1)[:]
	return verts[argmin(summed_dists)]
end

"""
	 distance(p1, p2, distances; method = CentroidToCentroid())

Find the distance between `Parcel`s `p1` and `p2` according to distance matrix
`distances`, using `method` (one of `CentroidToCentroid()` or `ClosestVertices`)
"""
function distance(
		p1::Parcel, p2::Parcel, distances::AbstractMatrix; 
		method::DistanceMetric = CentroidToCentroid()
	)
	distance(method, p1, p2, distances)
end

"""
	 distance(p1, p2; method = CentroidToCentroid())

Find the distance between `Parcel`s `p1` and `p2` using `method` (one of 
`CentroidToCentroid()` or `ClosestVertices`). This method call will expect to find a
distance matrix `:distances` belonging to the first parcel's `SurfaceSpace` struct
"""
function distance(p1::Parcel, p2::Parcel, method::DistanceMetric = CentroidToCentroid())
	distance(method, p1, p2, p1.surface[:distances])
end

function distance(
		::CentroidToCentroid, p1::Parcel, p2::Parcel, distances::AbstractMatrix
	)
	c1 = centroid(p1, distances)
	c2 = centroid(p2, distances)
	return distances[c1, c2]
end

function distance(
		::ClosestVertices, p1::Parcel, p2::Parcel, distances::AbstractMatrix
	)
	return minimum(distances[vertices(p1), vertices(p2)])
end


