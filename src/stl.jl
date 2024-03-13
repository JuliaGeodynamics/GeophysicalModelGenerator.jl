# NOTE: this remains WIP, which is why it is not yet in the documentation
# These are routines that allow importing *.stl triangular surfaces and employ them

using FileIO
using GeometryBasics: TriangleP, Mesh, normals, PointMeta, coordinates
using LinearAlgebra

# Warning: the TriangleIntersect dependency does not seem to work on different machines, as the developer did not add a version number..
# That forces us to remove it here, and
#using TriangleIntersect

export Ray, Intersection, IntersectRayTriangle, load, TriangleNormal, Point, IntersectRayMesh, coordinates
export STLToSurface, isInsideClosedSTL

#=
# Conversion routines from GeometryBasics triangles to TriangleIntersect triangles:
Base.convert(::Type{Point}, p::PointMeta)       =   Point(p[1],p[2],p[3])
Base.convert(::Type{Triangle}, t::TriangleP)    =   Triangle(convert(Point,t[1]),convert(Point,t[2]),convert(Point,t[3]))
Triangle(t::TriangleP) = convert(Triangle,t)


# Define a few routines to allow computing the intersection of a ray with a triangle
#/(p::GeometryBasics.Point, n::Number) = GeometryBasics.Point(p.x/n, p.y/n, p.z/n)
#unitize(p::GeometryBasics.Point) = p/(p*p)

#=
# Intersection
struct Intersection
    ip::Point    # intersecting point
    id::Float64                 # intersecting distance
    is_intersection::Bool
end

# Ray
struct Ray
    origin::Point
    direction::Point
end


const no_intersection = Intersection(Point(0,0,0), 0.0, false)

cross(p1::Any, p2::Any) = Point(p1[2]*p2[3]-p1[3]*p2[2], -p1[1]*p2[3]+p1[3]*p2[1], p1[1]*p2[2]-p1[2]*p2[1])


function triangleNormal(t::TriangleP)    # normal of a triangle
    a,b,c   =   t[1], t[2], t[3];
    v1      =   b-a
    v2      =   c-a
    normal  =   cross(v1, v2)'
    normal  =   normal/sqrt(sum(normal.^2))    # normalize


    return normal
end


# Taken from TriangleIntersect.jl & adapted for the GeometryBasics TriangleP
function intersectRayTriangle(r::Ray, t::TriangleP, normal)
    v1      =   t[2]-t[1];
    v2      =   t[3]-t[1];
    v1v1    =   dot(v1,v1)
    v2v2    =   dot(v2,v2)
    v1v2    =   dot(v1,v2)
    t_denom =   v1v2*v1v2 - v1v1*v2v2

    denom = dot(normal,r.direction)
    denom == 0 && return no_intersection
    ri = normal*(t[1] - r.origin) / denom
    ri <= 0 && return no_intersection
    plane_intersection =  ri * r.direction + r.origin
    w = plane_intersection - t[1]
    wv1 = w'*v1
    wv2 = w'*v2
    s_intersection = (v1v2*wv2 - v2v2*wv1) / t_denom
    s_intersection <= 0 && return no_intersection
    s_intersection >= 1 && return no_intersection
    t_intersection = (v1v2*wv1 - v1v1*wv2) / t_denom
    t_intersection <= 0 && return no_intersection
    t_intersection >= 1 && return no_intersection
    s_intersection + t_intersection >= 1 && return no_intersection
    return Intersection(t[1] + s_intersection*v1+t_intersection*v2, ri, true)


end
=#

"""
    ind, pts, dist = IntersectRayMesh(r::Ray, mesh::GeometryBasics.Mesh)

Intersects a ray with a triangular `*.stl` mesh
"""
function increasentersectRayMesh(r::Ray, mesh::Mesh)

    N = normals(mesh);

    Intersecting_Ind    =   [];
    Intersecting_Points =   [];
    Intersecting_Dist   =   [];
    for i=1:length(mesh)
        is = intersect(r, Triangle(mesh[i]))
        if is.is_intersection
            Intersecting_Ind    = [Intersecting_Ind;    i];
            Intersecting_Points = [Intersecting_Points; is.ip];
            Intersecting_Dist   = [Intersecting_Dist;   is.id];
        end
    end

    return Intersecting_Ind, Intersecting_Points, Intersecting_Dist
end


function stlToSurface(name::String,  Xin, Yin, minZ)

    mesh =  load(name)
    if length(size(Xin))==3
        X    =  Xin[:,:,1];
        Y    =  Yin[:,:,1];
    else
        X,Y  =  Xin, Yin;
    end

    Z    =  ones(size(X))*NaN

    minM  =  minimum.(mesh)
    maxM  =  maximum.(mesh)
    max_x = [maxM[i][1] for i=1:length(maxM)]
    min_x = [minM[i][1] for i=1:length(maxM)]
    max_y = [maxM[i][2] for i=1:length(maxM)]
    min_y = [minM[i][2] for i=1:length(maxM)]

    for i in eachindex(X)

        r_up =  Ray(Point(X[i],Y[i],minZ),Point(0.0, 0.0, 1.0))

        ind = findall(  (X[i] .>= min_x)  .& (X[i] .<= max_x) .&
                        (Y[i] .>= min_y)  .& (Y[i] .<= max_y) )

        for iT in eachindex(ind)

            is = intersect(r_up, Triangle(mesh[ind[iT]]))
            if is.is_intersection
                Z[i] = is.ip.z
            end


        end
    end

    X_out = zeros( size(X,1), size(X,2), 1)
    Y_out = zeros( size(X,1), size(X,2), 1)
    Z_out = zeros( size(X,1), size(X,2), 1)
    X_out[:,:,1] = X;
    Y_out[:,:,1] = Y;
    Z_out[:,:,1] = Z;

    Data = ParaviewData(X_out,Y_out,Z_out,(depth=Z_out,))

    return Data

end
=#


"""
    inside = isInsideClosedSTL(mesh::Mesh, Pt, eps=1e-3)

Determine whether a point `Pt` is inside a 3D closed triangular `*.stl` surface or not.

This implements the winding number method, following the python code:
https://github.com/marmakoide/inside-3d-mesh

This again is described in the following [paper](https://igl.ethz.ch/projects/winding-number/) by Alec Jacobson, Ladislav Kavan and Olga Sorkine-Hornung.
"""
function isInsideClosedSTL(mesh::Mesh, Pt::Vector, eps=1e-3)

     # Compute triangle vertices and their norms relative to X
     M_vec  = [mesh.position[i]-Pt[:]   for i in eachindex(mesh.position)];
     M      = zeros(length(M_vec),3);
     for i=1:length(M_vec)
        M[i,:] = Float64.(M_vec[i][1:3]);
     end
     M_norm = sqrt.(sum(M.^2,dims=2))

     # Accumulate generalized winding number per triangle
     winding_number = 0.
     for iT=1:length(mesh)
            t = mesh[iT].points;
            M = zeros(3,3)
            for i=1:3
                M[i,:] = t[i]-Pt[:];
            end
            M_norm = sqrt.(sum(M.^2,dims=1))

            A,B,C = M[1,:], M[2,:], M[3,:]
            a,b,c = M_norm[1], M_norm[2], M_norm[3]

         winding_number += atan(det(M), (a * b * c) + c * dot(A, B) + a * dot(B, C) + b * dot(C, A))
      end

     # Job done
     if winding_number >= 2pi - eps
        isinside = true;
     else
        isinside = false
     end


     return isinside
end
