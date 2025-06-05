# This generates input model setups for pTatin (an open source parallel finite element code for geodynamic application)

# We follow the python routines given here https://github.com/hpc4geo/gmsh_to_point_cloud, written by Dave May

import Base: show, size
export write_pTatin_mesh, swap_yz_dims


"""
    tags = cell_tags_from_gmsh(mesh::GmshDiscreteModel)
Returns a list with integers that are the tags for each of the cells
"""
function cell_tags_from_gmsh(mesh)
    cell_entities = mesh.face_labeling.d_to_dface_to_entity[4]
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

    return tags
end


"""
    write_pTatin_mesh(fe_mesh::FEData; out_file="md.bin", connectivity_zero_based=true)

Write a binary file with the mesh information for pTatin
"""
function write_pTatin_mesh(fe_mesh::FEData; out_file = "md.bin", connectivity_zero_based = true)

    # Write mesh
    write_FEmesh_mesh(fe_mesh; out_file = out_file, connectivity_zero_based = connectivity_zero_based)

    # Write cell and vertex fields
    write_FEmesh_fields(fe_mesh)

    return nothing
end

"""
    write_pTatin_mesh(q1_mesh::Q1Data; out_file="md.bin", connectivity_zero_based=true)

Write a binary file with the mesh information for pTatin
"""
write_pTatin_mesh(q1_mesh::Q1Data; out_file = "md.bin", connectivity_zero_based = true) = write_pTatin_mesh(convert2FEData(q1_mesh), out_file = out_file, connectivity_zero_based = connectivity_zero_based)


"""
    write_FEmesh_mesh(vdata::FEData; out_file="md.bin", connectivity_zero_based=true)

Writes a binary file with the mesh information for pTatin
"""
function write_FEmesh_mesh(data::FEData; out_file = "md.bin", connectivity_zero_based = true)
    dims = size(data.vertices, 1)
    nverts = size(data.vertices, 2)
    arange = 0:(nverts - 1)
    @assert nverts > dims

    nvertices_cell = size(data.connectivity, 1)
    ncells = size(data.connectivity, 2)
    @assert ncells > nvertices_cell

    shift_connectivity = 0
    if connectivity_zero_based
        shift_connectivity = 1
    end

    npart = 1 # don't partition currenrtly

    f = open(out_file, "w")
    write(f, Int32(nverts))              # nverts
    write(f, Int32(dims))                # coordinate_dimension
    write(f, data.vertices[:])                # vertices

    write(f, Int32(ncells))                  # ncells
    write(f, Int32(nvertices_cell))          # points-per-cell
    write(f, Int32.(data.connectivity[:] .- shift_connectivity))   # cells

    write(f, Int32(npart))               # npartitions
    write(f, minimum(data.vertices, dims = 2))   # bbmin
    write(f, maximum(data.vertices, dims = 2))   # bbmax

    write(f, Int32(ncells))              # ncells
    write(f, Int32.(arange))             # arange

    close(f)
    println("Wrote pTatin mesh  : $(out_file)")

    return nothing
end


"""
    write_FEmesh_fields(data::FEData)

Writes cell and vertex fields to disk to be read by pTatin
"""
function write_FEmesh_fields(data::FEData)

    ncells = size(data.connectivity, 2)
    nverts = size(data.vertices, 2)

    # Mapping used internally
    dtype_map = Dict()
    dtype_map[:Int16] = 10
    dtype_map[:Int32] = 11
    dtype_map[:Int64] = 12
    dtype_map[:Float32] = 20
    dtype_map[:Float64] = 21

    # Write the cell fields
    if !isnothing(data.cellfields)
        names = keys(data.cellfields)
        for i in 1:length(data.cellfields)
            _write_field_file(String(names[i]) * "_cell.bin", data.cellfields[i], ncells, dtype_map)
        end
    end

    # Write the vertex fields
    if !isnothing(data.fields)
        names = keys(data.fields)
        for i in 1:length(data.fields)
            _write_field_file(String(names[i]) * "_vertex.bin", data.fields[i], nverts, dtype_map)
        end
    end

    return nothing
end

# Writes a field to disk
function _write_field_file(fname::String, field, len, dtype_map)
    @assert length(field) == len

    f = open(fname, "w")
    write(f, Int32(len))                     # length
    write(f, Int32(dtype_map[Symbol(eltype(field))])) # field type
    write(f, field[:])                        # field
    close(f)
    println("Wrote pTatin field : $fname")

    return nothing
end

"""
    fe_swap = swap_yz_dims(fe_data::FEData)

This swaps the `y` and `z` dimensions of the FEData object, which is useful for pTatin as it uses `y` for what is `z` in GMG.
"""
function swap_yz_dims(fe_data::FEData)
    vertices = copy(fe_data.vertices)
    return FEData(vertices[[1, 3, 2], :], fe_data.connectivity, fe_data.fields, fe_data.cellfields)
end
