using Test, GeophysicalModelGenerator

# 3D arrays
c_1D = CartData(XYZGrid(1:4,0,0))
c_2D = CartData(XYZGrid(1:4,1:5,2))
c_3D = CartData(XYZGrid(1:4,1:5,2:5))

# points
X_pt, Y_pt, Z_pt = XYZGrid(1:.05:5,0:.07:8,1:.4:5)

# 1D test
id_1D = nearest_point_indices(NumValue(c_1D.x),X_pt[:])
@test sum(id_1D) == 166336

# 2D test
id_2D = nearest_point_indices(NumValue(c_2D.x),NumValue(c_2D.y), X_pt[:], Y_pt[:])
@test sum(id_2D) == 992141

# 3D test
id_3D = nearest_point_indices(NumValue(c_2D.x),NumValue(c_2D.y),NumValue(c_2D.z), X_pt[:],Y_pt[:],Z_pt[:])
@test sum(id_3D) == 442556