# NOTE: this remains WIP, which is why it is not yet in the documentation
# These are routines that allow importing *.stl triangular surfaces and employ them

using FileIO
# using GeometryBasics: TriangleP, Mesh, normals, PointMeta, coordinates
# I removed TriangleP and PointMeta, as they are not existing in GeometryBasics anymore (past 0.5.0)
# HD, 2024-06-10
using GeometryBasics: Mesh, coordinates
using LinearAlgebra

# Warning: the TriangleIntersect dependency does not seem to work on different machines, as the developer did not add a version number..
# That forces us to remove it here, and
#using TriangleIntersect

export load, coordinates
export isinside_closed_STL

"""
    inside = isinside_closed_STL(mesh::Mesh, Pt, eps=1e-3)

Determine whether a point `Pt` is inside a 3D closed triangular `*.stl` surface or not.

This implements the winding number method, following the python code:
https://github.com/marmakoide/inside-3d-mesh

This again is described in the following [paper](https://igl.ethz.ch/projects/winding-number/) by Alec Jacobson, Ladislav Kavan and Olga Sorkine-Hornung.
"""
function isinside_closed_STL(mesh::Mesh, Pt::Vector, eps = 1.0e-3)

    # Accumulate generalized winding number per triangle
    winding_number = 0.0
    for iT in 1:length(mesh)
        t = mesh[iT].points
        M = zeros(3, 3)
        for i in 1:3
            M[i, :] = t[i] - Pt[:]
        end
        M_norm = sqrt.(sum(M .^ 2, dims = 2))  # per vertex (row)

        A, B, C = M[1, :], M[2, :], M[3, :]
        a, b, c = M_norm[1], M_norm[2], M_norm[3]

        winding_number += atan(det(M), (a * b * c) + c * dot(A, B) + a * dot(B, C) + b * dot(C, A))
    end

    # Job done
    if winding_number >= 2pi - eps
        isinside = true
    else
        isinside = false
    end


    return isinside
end
