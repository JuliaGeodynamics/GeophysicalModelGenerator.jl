using Test
using GeophysicalModelGenerator
using VectorizedRoutines
include("../src/voxel_gravity.jl")

function main()
############# Input #############
# survey
x        = [-20.,   100.]
y        = [100., 200.]
z        = [-100., 10.]
nx       = 111
ny       = 101
nz       = 91

# spheres
centerX  = [  50.,   10.,   80.]
centerY  = [ 150.,  150.,  180.]
centerZ  = [ -30.,  -10.,  -50.]
radius   = [  10.,    3.,   10.]
rho      = [2650., 2600., 2750.]

# background density model
background = [0. 2700.]
#################################

######### Set things up #########
# constants
G        = 6.67408e-11
mGal     = 1e5

# coordinates
x_vec    = LinRange(x[1],x[2],nx)
y_vec    = LinRange(y[1],y[2],ny)
z_vec    = LinRange(z[1],z[2],nz)
X,Y,Z    = Matlab.meshgrid(x_vec,y_vec,z_vec)

# reference model
RefMod   = zeros(nz)
RefMod[findall(x->x>=0, z_vec)] .= background[1]
RefMod[findall(x->x<0,  z_vec)] .= background[2]

# densities
RHO      = zeros(ny,nx,nz)
RHO[findall(x->x>=0, Z)] .= background[1]
RHO[findall(x->x<0,  Z)] .= background[2]

# put in spheres
if !(length(centerX) == length(centerY) && length(centerX) == length(centerZ) && length(centerX) == length(radius))
    error("Sphere definition is wrong")
end

numSpheres = length(rho)

for i = 1 : numSpheres
    d = (X .- centerX[i]).^2 .+ (Y .- centerY[i]).^2 .+ (Z .- centerZ[i]).^2
    RHO[findall(x->x<radius[i]^2, d)] .= rho[i]
end
#################################

############ Compute ############
# Analytical Solution
ana = zeros(ny,nx)
for iS = 1 : numSpheres
    for iX = 1 : nx
        for iY = 1 : ny
            d          = ((x_vec[iX] - centerX[iS])^2 + (y_vec[iY] - centerY[iS])^2 + (centerZ[iS])^2)^0.5
            depth      = -centerZ[iS]
            ana[iY,iX] = ana[iY,iX] + (4 * π * G * (radius[iS]^3) * (rho[iS] - background[2]) * depth) / (3 * (d^3)) * mGal
        end
    end
end

# Voxel Code
dg1, gradX1, gradY1 = voxGrav(X, Y, Z, RHO, refMod = RefMod, outName = "Benchmark1", printing = false)

# Other options (same result)
dg2, gradX2, gradY2 = voxGrav(X, Y, Z, RHO, rhoTol = 1, refMod = "SW", outName = "Benchmark2", printing = false)

# Other options (different result)
dg3, gradX3, gradY3 = voxGrav(X, Y, Z, RHO, rhoTol = 70, refMod = "NW", outName = "Benchmark3", printing = false)
#################################

############# Test ##############
maxVal = maximum(broadcast(abs, ana))
maxErr = maxVal / 20

# check differences
diff1 = broadcast(abs, ana-dg1)
diff2 = broadcast(abs, ana-dg2)
diff3 = broadcast(abs, ana-dg3)

# check
@test maximum(diff1) < maxErr
@test maximum(diff2) < maxErr
@test maximum(diff3) > maxErr
#################################
end

main()