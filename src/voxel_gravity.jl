# This is voxel_gravity.jl
#
# This function allows the computation of bouguer anomalies over a regular grid
#
# Author: Arne Spang 06/2021

using Printf            # do we need this?
using Statistics        # do we need that?

export voxGrav

"""
    voxGrav(X::Array{Float64, 3}, Y::Array{Float64, 3}, Z::Array{Float64, 3}, RHO::Array{Float64, 3};
    refMod="AVG", lengthUnit="m", rhoTol=1e-9, Topo=[], outName="Bouguer", printing=true)

    Computes Bouguer anomalies and gradients

    Required arguments:
    X,Y,Z:       3D matrices with the coordinates of the grid (X should vary in the first dimension, Y in the second, Z (vertical) in the third)
    RHO:         3D matrix with the density at each grid point [kg/m^3]

    Optional arguments:
    refMod:      1D vector with the reference density for each depth
                 Alternatively, the strings "NE", "SE", "SW", "NW", "AVG" can be used.
                 In that case, one of the corners of `RHO` is used as reference model.
                 In case of "AVG" the reference model is the average of each depth slice.
    lengthUnit:  The unit of the coordinates and topography file. Either "m" or "km"
    rhoTol:      density differences smaller than rhoTol will be ignored [kg/m^3]
    Topo:        2D matrix with the topography of the surface (only relevant for the paraview output)
    outName:     name of the paraview output (do not include file type)
    printing:    activate printing of additional information [true or false]
"""
function voxGrav(X::Array{Float64, 3}, Y::Array{Float64, 3}, Z::Array{Float64, 3}, RHO::Array{Float64, 3};
                 refMod="AVG", lengthUnit="m", rhoTol=1e-9, Topo=[], outName="Bouguer", printing=true)

    ## check input
    X, Y, Z, RHO, RefMod, rhoTol, Topo, orient = checkInput(X, Y, Z, RHO, refMod, lengthUnit, rhoTol, Topo, outName, printing)

    ################ precompute things ################
    # define constants
    G      = 6.67408e-11

    # get coordinate vectors
    x_vec  = X[:,1,1]
    y_vec  = Y[1,:,1]
    z_vec  = Z[1,1,:]

    # cut everything above sea level
    ind    = findall(x->x<=0, z_vec)
    X      = X[:,:,ind]
    Y      = Y[:,:,ind]
    Z      = Z[:,:,ind]
    RHO    = RHO[:,:,ind]
    RefMod = RefMod[ind]

    # check dimensions
    nx     = size(X,1);
    ny     = size(X,2);
    nz     = size(X,3);

    # subtract reference model
    for i = 1 : nz
        RHO[:,:,i] .= RHO[:,:,i] .- RefMod[i]
    end

    # interpolate density grid to cell centers
    DRHO = RHO[1:end-1,1:end-1,1:end-1] .+ RHO[2:end,1:end-1,1:end-1] .+ RHO[2:end,2:end,1:end-1] .+ RHO[1:end-1,2:end,1:end-1] +
            RHO[1:end-1,1:end-1,2:end]   .+ RHO[2:end,1:end-1,2:end]   .+ RHO[2:end,2:end,2:end]   .+ RHO[1:end-1,2:end,2:end]
    DRHO = DRHO ./ 8

    # voxel volume
    dx     = X[2,1,1] - X[1,1,1]
    dy     = Y[1,2,1] - Y[1,1,1]
    dz     = Z[1,1,2] - Z[1,1,1]
    dV     = abs(dx*dy*dz)

    # dV * G is a constants
    VG_vox = dV * G

    # coordinate vector of cell center grid
    xCells = x_vec[1:end-1] .+ dx/2
    yCells = y_vec[1:end-1] .+ dy/2
    zCells = z_vec[1:end-1] .+ dz/2

    # precompute distances
    if printing
        @printf "Precomputing Distances:"
        @time d_cube = precompDist(x_vec, y_vec, xCells, yCells, zCells, nx, ny, nz)
        @printf "\n"
    else
        d_cube = precompDist(x_vec, y_vec, xCells, yCells, zCells, nx, ny, nz)
    end
    ###################################################

    ############# compute bouguer anomaly #############
    if printing
        @printf "Computing Bouguer anomaly:"
        @time dg = computeBoug(nx,ny,nz,DRHO,d_cube,VG_vox,zCells,rhoTol)
        @printf "\n"
    else
        dg = computeBoug(nx,ny,nz,DRHO,d_cube,VG_vox,zCells,rhoTol)
    end
    ###################################################

    ############ compute bouguer gradients ############
    if printing
        @printf "Computing Bouguer gradients:"
        @time gradX, gradY = computeBougGrads(nx,ny,dg)
        @printf "\n"
    else
        gradX, gradY = computeBougGrads(nx,ny,dg)
    end
    ###################################################

    ############## write info and output ##############
    numRel = length(findall(x->abs(x)>rhoTol, DRHO))
    frac   = numRel/length(DRHO)
    coords = cat(X[:,:,1],Y[:,:,1],Topo,dims=4)
    coords = permutedims(coords,[4,1,2,3])
    vtkfile = vtk_grid(outName,coords)

    if printing
        @printf "%.3f %% of the domain contained anomalous densities. If this is more than expected, adjust rhoTol for faster computation.\n" 100*frac
        @printf "Writing output...\n"
        vtkfile["Bouguer Anomaly [mGal]"]   = dg
        vtkfile["Boug_Gradient_X [mGal/m]"] = gradX
        vtkfile["Boug_Gradient_Y [mGal/m]"] = gradY
        outfiles = vtk_save(vtkfile)
        @printf "Wrote output to: %s\n\n" outfiles
    else
        vtkfile["Bouguer Anomaly [mGal]"]   = dg
        vtkfile["Boug_Gradient_X [mGal/m]"] = gradX
        vtkfile["Boug_Gradient_Y [mGal/m]"] = gradY
        outfiles = vtk_save(vtkfile)
    end
    ###################################################

    if orient == 1
        return dg, gradX, gradY
    else
        return permutedims(dg, [2,1]), permutedims(gradX, [2,1]), permutedims(gradY, [2,1])
    end
end





function checkInput(X, Y, Z, RHO, refMod, lengthUnit, rhoTol, Topo, outName, printing)
    # orientation
    if X[1,1,1] ≠ X[2,1,1] && X[1,1,1] == X[1,2,1] && X[1,1,1] == X[1,1,2]
        orientation = 1
    elseif X[1,1,1] == X[2,1,1] && X[1,1,1] ≠ X[1,2,1] && X[1,1,1] == X[1,1,2]
        orientation = 2
        X   = permutedims(X,   [2,1,3])
        Y   = permutedims(Y,   [2,1,3])
        Z   = permutedims(Z,   [2,1,3])
        RHO = permutedims(RHO, [2,1,3])
    else
        error("Coordinate orientation looks wrong!")
    end

    # dimensions
    if !(size(X) == size(Y) && size(X) == size(Z) && size(X) == size(RHO))
        error("X, Y, Z, RHO must be 3D matrices of the same size.")
    end
    nz = size(X,3)

    # check if grid is regular
    dx = diff(X,dims=1)[:]
    dy = diff(Y,dims=2)[:]
    dz = diff(Z,dims=3)[:]
    tol = 1e-12
    if !(all(a->a<dx[1]+tol && a > dx[1]-tol,dx) && all(a->a<dy[1]+tol && a > dy[1]-tol,dy) && all(a->a<dz[1]+tol && a > dz[1]-tol,dz))
        error("Non-regular grids are not supported yet")
    end

    # Topo
    if !isempty(Topo)
        if !(size(Topo,1) == size(X,1) && size(Topo,2) == size(X,2))
            if printing
                @printf "Topo input dimensions do not fit X, Y, Z, RHO. Using a flat topography.\n"
            end
            Topo = zeros(size(X,1),size(X,2))
        end
    else
        Topo = zeros(size(X,1),size(X,2))
    end

    # tolerance
    if printing
        @printf "Using density tolerance: %.3e\n" rhoTol
    end

    # length unit
    if lengthUnit == "m" || lengthUnit == "M"
        if printing
            @printf "Assuming coordinates to be in meters.\n\n"
        end
    elseif lengthUnit == "km" || lengthUnit == "KM"
        if printing
            @printf "Converting coordinates to meters.\n\n"
        end
        X     = X    .* 1000;
        Y     = Y    .* 1000;
        Z     = Z    .* 1000;
        Topo  = Topo .* 1000;
    else
        error("lengthUnit should be \"m\" or \"km\".")
    end

    # reference model
    if !(typeof(refMod) == String)
        if length(refMod) == nz
            RefMod = refMod
        else
            error("RefMod must have the same length as the third dimension of X, Y, Z, RHO.")
        end
    else
        if refMod == "NE"
            RefMod = RHO[end,end,:]
        elseif refMod == "SE"
            RefMod = RHO[end,1,:]
        elseif refMod == "SW"
            RefMod = RHO[1,1,:]
        elseif refMod == "NW"
            RefMod = RHO[1,end,:]
        elseif refMod == "AVG"
            RefMod = !mean([1.,1.,1.],RHO)
        else
            error("RefMod should be NE, SE, SW, NW, AVG or a vector with one value for each depth.")
        end
    end

    # outName
    if !(typeof(outName) == String)
        error("outName must be a string.")
    end

    return X, Y, Z, RHO, RefMod, rhoTol, Topo, orientation

end

function precompDist(x_vec, y_vec, xCells, yCells, zCells, nx, ny, nz)
    d_square = zeros(nx-1,ny-1,nz-1)

    for iX = 1 : nx - 1
        for iY = 1 : ny - 1
            for iZ = 1 : nz - 1
                d_square[iX,iY,iZ] = (xCells[iX]-x_vec[1])^2 + (yCells[iY]-y_vec[1])^2 + zCells[iZ]^2
            end
        end
    end

    d_cube = d_square .^1.5

    return d_cube
end

function computeBoug(nx,ny,nz,DRHO,d_cube,VG_vox,zCells,rhoTol)
    mGal   = 1e5

    dg     = zeros(nx,ny)

    for jX = 1 : nx-1
        for jY = 1 : ny-1
            for jZ = 1 : nz-1
                if (DRHO[jX,jY,jZ] > rhoTol || DRHO[jX,jY,jZ] < -rhoTol)
                    for iX = 1 : nx
                        for iY = 1 : ny
                            d_indX = jX-iX
                            if d_indX < 0; d_indX = abs(d_indX); else d_indX = d_indX + 1; end
                            d_indY = jY-iY;
                            if d_indY < 0; d_indY = abs(d_indY); else d_indY = d_indY + 1; end
                            d_c       = d_cube[d_indX,d_indY,jZ];
                            dg[iX,iY] = dg[iX,iY] + VG_vox * DRHO[jX,jY,jZ] * -zCells[jZ]/d_c;
                        end
                    end
                end
            end
        end
    end

    dg = dg * mGal

    return dg
end

function computeBougGrads(nx,ny,dg)
    gradX  = zeros(nx,ny)
    gradY  = zeros(nx,ny)

    for iY = 1 : ny
        itp  = Interpolations.interpolate(dg[:,iY], BSpline(Quadratic(Reflect(OnCell()))))
        for iX = 1 : nx
            grad         = Interpolations.gradient(itp,iX)
            gradX[iX,iY] = grad[1]
        end
    end

    for iX = 1 : nx
        itp  = Interpolations.interpolate(dg[iX,:], BSpline(Quadratic(Reflect(OnCell()))))
        for iY = 1 : ny
            grad         = Interpolations.gradient(itp,iY)
            gradY[iX,iY] = grad[1]
        end
    end

    return gradX, gradY
end
