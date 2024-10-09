#helper functions to make GMG work with Chmy grids and fields
module Chmy_utils

using Chmy, Chmy.Grids, Chmy.Fields

import GeophysicalModelGenerator: create_CartGrid, CartGrid, CartData
import GeophysicalModelGenerator: add_box!, add_sphere!, add_ellipsoid!, add_cylinder!
import GeophysicalModelGenerator: add_layer!, add_polygon!, add_slab!, add_stripes!, add_volcano!
import GeophysicalModelGenerator: above_surface, below_surface

println("Loading Chmy-GMG tools")

"""
    CartGrid = create_CartGrid(grid::StructuredGrid; ylevel=0.0)

Creates a GMG `CartGrid` data structure from a `Chmy` grid object
"""
function create_CartGrid(grid::StructuredGrid; ylevel=0.0)

    coord1D     = Vector.(coords(grid, Vertex()))
    coord1D_cen = Vector.(coords(grid, Center()))
    N           = length.(coord1D)
    L           = extent(grid, Vertex())
    X₁          = origin(grid, Vertex())  
    Δ           = spacing(grid)
    ConstantΔ   = false;
    if isa(grid, UniformGrid)
        ConstantΔ = true
    end
    if ndims(grid)==2
        # we need a special treatment of this, as all GMG routines work with 3D coordinates
        X₁  = (X₁[1], ylevel, X₁[2]) 
        L   = (L[1], 0.0, L[2]) 
        Δ   = (Δ[1], 0.0, Δ[2])
        N   = (N[1],1,N[2])
        coord1D = (coord1D[1], [0.0], coord1D[2])
        coord1D_cen = (coord1D_cen[1], [0.0], coord1D_cen[2])
    end
    Xₙ          = X₁ .+ L
    

    return CartGrid(ConstantΔ,N,Δ,L,X₁,Xₙ,coord1D, coord1D_cen)
end

# all functions to be ported
function_names = (:add_box!, :add_sphere!, :add_ellipsoid!, :add_cylinder!, :add_layer!, :add_polygon!, :add_slab!, :add_stripes!, :add_volcano!)

for fn in function_names

    @eval begin 
        """
            $($fn)( Phase::Field, 
                    Temp::Field,
                    Grid::StructuredGrid;      # required input
                    kwargs...) 

        Sets `$($fn)` function for `Chmy` fields and grids.           
        """
        function $fn(  Phase::Field, 
                            Temp::Field,
                            Grid::StructuredGrid;      # required input
                            kwargs...)  

            CartGrid = create_CartGrid(Grid)

            cell = false
            if all(location(Phase).==Center()) 
                cell = true
            end
            
            return ($fn)(Phase, Temp, CartGrid; cell=cell, kwargs...)
        end
    end
    
end


# all functions to be ported
function_names = (:above_surface, :below_surface)

for fn in function_names

    @eval begin 
        """
            $($fn)( Grid::StructuredGrid, field::Field, DataSurface_Cart::CartData; kwargs...) 

        Sets `$($fn)` function for `Chmy` grids and the field `field` which can be either on vertices or centers           
        """
        function $fn( Grid::StructuredGrid,
                      field::Field, 
                      DataSurface_Cart::CartData;     
                      kwargs...)  

            CartGrid = create_CartGrid(Grid)
            
            cell = false
            if all(location(field).==Center()) 
                cell = true
            end

            return ($fn)(CartGrid, DataSurface_Cart; cell=cell, kwargs...)
        end
    end
    
end


end # end of module
