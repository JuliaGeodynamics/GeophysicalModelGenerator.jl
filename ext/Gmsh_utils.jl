# Support for Gmsh meshes

module Gmsh_utils

using GridapGmsh, StaticArrays

import GeophysicalModelGenerator: import_Gmsh, FEData, CartData


"""
    fe_data, tag_names = import_Gmsh(fname::String)

Reads a Gmsh file and returns a `FEData` object with info about the mesh. `tag_names` contains the names of the regions of the FE mesh
"""
function import_Gmsh(fname::String)

    mesh = GmshDiscreteModel(fname, renumber = false)

    # Extract vertices
    nverts = length(mesh.grid.node_coordinates)
    dims = length(mesh.grid.node_coordinates[1])
    vertices = [mesh.grid.node_coordinates[n][i] for i in 1:dims,n in 1:nverts]

    # write coords as 1D double array
    nvertices_cell = length(mesh.grid.cell_node_ids[1])
    connectivity = [c[i] for i in 1:nvertices_cell, c in mesh.grid.cell_node_ids]

    # extract tag of each of the tetrahedrons
    regions, tag_names = cell_tags_from_gmsh(mesh)

    cellfields = (regions = regions,)
    fields = nothing

    return FEData(vertices, connectivity, fields, cellfields), tag_names
end


"""
    tags, tag_names = cell_tags_from_gmsh(mesh::GmshDiscreteModel)

Returns a list with integers that are the tags for each of the cells
"""
function cell_tags_from_gmsh(mesh)
    cell_entities = mesh.face_labeling.d_to_dface_to_entity[4]    # volumetric entities
    cell_entities_unique = unique(cell_entities)
    tag_unique = zeros(Int64, size(cell_entities_unique))

    for i in 1:length(cell_entities_unique)
        for (n, tag) in enumerate(mesh.face_labeling.tag_to_entities)
            if any(tag .== cell_entities_unique[i])
                tag_unique[i] = n
            end
        end
    end

    # create tags for cells
    tags = zeros(Int64, length(cell_entities))
    for (i, entity) in enumerate(cell_entities_unique)
        id = findall(cell_entities .== entity)
        tags[id] .= tag_unique[i]
    end

    return tags, mesh.face_labeling.tag_to_name
end


end
