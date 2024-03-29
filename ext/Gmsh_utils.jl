# Support for Gmsh meshes

module Gmsh_utils

using GridapGmsh, StaticArrays

import GeophysicalModelGenerator: import_Gmsh, FEData


"""
    fe_data, tag_names = import_Gmsh(fname::String)

Reads a Gmsh file and returns a `FEData` object with info about the mesh. `tag_names` contains the names of the regions of the FE mesh
"""
function import_Gmsh(fname::String)

    mesh = GmshDiscreteModel(fname, renumber=false)

    # Extract vertices
    nverts   = length(mesh.grid.node_coordinates);
    dims     = length(mesh.grid.node_coordinates[1])
    vertices = [mesh.grid.node_coordinates[n][i] for i=1:dims,n=1:nverts]

    # write coords as 1D double array
    nvertices_cell = length(mesh.grid.cell_node_ids[1])
    connectivity   = [c[i] for i=1:nvertices_cell, c in mesh.grid.cell_node_ids]

    # extract tag of each of the tetrahedrons
    regions, tag_names = cell_tags_from_gmsh(mesh)     

    cellfields  = (regions=regions,)
    fields      = nothing

    return FEData(vertices,connectivity, fields, cellfields), tag_names
end


"""
    tags, tag_names = cell_tags_from_gmsh(mesh::GmshDiscreteModel)

Returns a list with integers that are the tags for each of the cells
"""
function cell_tags_from_gmsh(mesh)
    cell_entities   = mesh.face_labeling.d_to_dface_to_entity[4]    # volumetric entities
    cell_entities_unique = unique(cell_entities)
    tag_unique = zeros(Int64,size(cell_entities_unique))
    
    for i=1:length(cell_entities_unique)
        for (n,tag) in enumerate(mesh.face_labeling.tag_to_entities)
           if any(tag .== cell_entities_unique[i])
               tag_unique[i] = n
           end
        end
    end

    # create tags for cells
    tags = zeros(Int64,length(cell_entities))
    for (i,entity) in enumerate(cell_entities_unique)
        id = findall(cell_entities.==entity)
        tags[id] .= tag_unique[i]
    end
   
    return tags,  mesh.face_labeling.tag_to_name
end





"""
    inside = point_in_tetrahedron(p::_T, a::_T, b::_T, c::_T, d::_T, tol=1e-10)
Determines if a point `p` is inside a tetrahedron specified by `a`,`b`,`c`,`d` or not    
"""
function point_in_tetrahedron(p::_T, a::_T, b::_T, c::_T, d::_T, tol=1e-10) where _T<:Vector{Float64}

    # check bounding box
    xmin = min(a[1],b[1],c[1],d[1])
    xmax = max(a[1],b[1],c[1],d[1])
    ymin = min(a[2],b[2],c[2],d[2])
    ymax = max(a[2],b[2],c[2],d[2])
    zmin = min(a[3],b[3],c[3],d[3])
    zmax = max(a[3],b[3],c[3],d[3])
    
    inside = true
    if p[1] < xmin || p[1] > xmax 
        inside = false
    end
    if (p[2] < ymin || p[2] > ymax) && inside  
        inside = false
    end
    if (p[3] < zmin || p[3] > zmax) && inside  
        inside = false
    end
    
    if inside
        v0 = @SVector [d[i] - a[i] for i in 1:3]
        v1 = @SVector [b[i] - a[i] for i in 1:3]
        v2 = @SVector [c[i] - a[i] for i in 1:3]
        v3 = @SVector [p[i] - a[i] for i in 1:3]
    
        denom = dot(v0, cross(v1, v2))
    
        u = dot(v3, cross(v1, v2)) / denom
        v = dot(v0, cross(v3, v2)) / denom
        w = dot(v0, cross(v1, v3)) / denom
    
        inside =  (u >= -tol) && (v >= -tol) && (w >= -tol) && (u + v + w <= 1 + tol)
    end

    return inside
end


end