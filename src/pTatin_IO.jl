# This generates input model setups for pTatin (an open source parallel finite element code for geodynamic application)

# We follow the python routines given here https://github.com/hpc4geo/gmsh_to_point_cloud, written by Dave May

import Base: show, size
export write_pTatin_mesh, swap_yz_dims


"""
    write_pTatin_mesh(fe_mesh::FEData; out_file="md.bin", connectivity_zero_based=true)

Write a binary file with the mesh information for pTatin, and one file per field
(`<name>_cell.bin`, `<name>_vertex.bin`) in the same directory as `out_file`.
"""
function write_pTatin_mesh(fe_mesh::FEData; out_file = "md.bin", connectivity_zero_based = true)

    # Write mesh
    write_FEmesh_mesh(fe_mesh; out_file, connectivity_zero_based)

    # Write cell and vertex fields
    write_FEmesh_fields(fe_mesh; dir = dirname(out_file))

    return nothing
end

"""
    write_pTatin_mesh(q1_mesh::Q1Data; out_file="md.bin", connectivity_zero_based=true)

Write a binary file with the mesh information for pTatin, and one file per field in the same
directory as `out_file`.
"""
write_pTatin_mesh(q1_mesh::Q1Data; out_file = "md.bin", connectivity_zero_based = true) = write_pTatin_mesh(convert2FEData(q1_mesh); out_file, connectivity_zero_based)


"""
    write_FEmesh_mesh(vdata::FEData; out_file="md.bin", connectivity_zero_based=true)

Writes a binary file with the mesh information for pTatin
"""
function write_FEmesh_mesh(data::FEData; out_file = "md.bin", connectivity_zero_based = true)
    dims = size(data.vertices, 1)
    nverts = size(data.vertices, 2)
    arange = 0:(nverts - 1)

    nvertices_cell = size(data.connectivity, 1)
    ncells = size(data.connectivity, 2)
    (nverts >= 1 && ncells >= 1) || error("the mesh has $nverts vertices and $ncells cells; it needs at least one of each")
    all(in(1:nverts), data.connectivity) || error("the connectivity refers to vertices outside 1:$nverts")

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
    write_FEmesh_fields(data::FEData; dir=".")

Writes cell and vertex fields to `<name>_cell.bin` and `<name>_vertex.bin` in directory `dir`,
to be read by pTatin
"""
function write_FEmesh_fields(data::FEData; dir = ".")

    ncells = size(data.connectivity, 2)
    nverts = size(data.vertices, 2)

    # Mapping used internally
    dtype_map = Dict()
    dtype_map[:Int16] = 10
    dtype_map[:Int32] = 11
    dtype_map[:Int64] = 12
    dtype_map[:Float32] = 20
    dtype_map[:Float64] = 21

    for (name, field) in pairs(data.cellfields)
        _write_field_file(joinpath(dir, "$(name)_cell.bin"), field, ncells, dtype_map)
    end
    for (name, field) in pairs(data.fields)
        _write_field_file(joinpath(dir, "$(name)_vertex.bin"), field, nverts, dtype_map)
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

Swapping two axes mirrors the mesh, so the vertices of each cell are reordered to keep the
cells' orientation. This is supported for 3D tetrahedra and for 3D hexahedra whose vertices
are numbered lexicographically (`x` fastest, then `y`, then `z`), as `convert2FEData` does.
"""
function swap_yz_dims(fe_data::FEData)
    nvertices_cell = size(fe_data.connectivity, 1)
    size(fe_data.vertices, 1) == 3 || error("swap_yz_dims needs 3D vertices, got $(size(fe_data.vertices, 1))D")
    # an odd permutation of a tetrahedron's vertices restores its orientation; for a hexahedron
    # this one swaps the local y and z directions
    perm = nvertices_cell == 4 ? [1, 3, 2, 4] :
        nvertices_cell == 8 ? [1, 2, 5, 6, 3, 4, 7, 8] :
        error("swap_yz_dims supports tetrahedra (4 vertices) and hexahedra (8 vertices) only, got $nvertices_cell vertices per cell")
    return FEData(fe_data.vertices[[1, 3, 2], :], fe_data.connectivity[perm, :], fe_data.fields, fe_data.cellfields)
end
